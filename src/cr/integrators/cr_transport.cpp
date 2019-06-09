//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file rad_transport.cpp
//  \brief implementation of radiation integrators
//======================================================================================


// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"
#include "../cr.hpp"
#include "../../coordinates/coordinates.hpp" //
#include <algorithm>   // min,max

// class header
#include "cr_integrators.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif


void CRIntegrator::CalculateFluxes(MeshBlock *pmb,
      AthenaArray<Real> &w, AthenaArray<Real> &bcc, 
      AthenaArray<Real> &u_cr, int reconstruct_order)
{
  CosmicRay *pcr=pmy_cr;
  Coordinates *pco = pmb->pcoord;
  
  Real invlim = 1.0/pcr->vmax;
  
  AthenaArray<Real> &x1flux=pcr->flux[X1DIR];
  AthenaArray<Real> &x2flux=pcr->flux[X2DIR];
  AthenaArray<Real> &x3flux=pcr->flux[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;

  // First, calculate com-moving flux
  int n1 = pmb->block_size.nx1 + 2*NGHOST;
  int n2 = pmb->block_size.nx2;
  int n3 = pmb->block_size.nx3;
  if(pmb->block_size.nx2 > 1) n2 += 2*NGHOST;
  if(pmb->block_size.nx3 > 1) n3 += 2*NGHOST;  
  

  int tid=0;
  int nthreads = pmb->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) private(tid) num_threads(nthreads)
{
#ifdef OPENMP_PARALLEL
    tid=omp_get_thread_num();
#endif



    AthenaArray<Real> flx, vel_l, vel_r, wl, wr, vdiff_l, vdiff_r, eddl, eddr;
    AthenaArray<Real> cwidth;
    flx.InitWithShallowSlice(flx_,3,tid,1);
    vel_l.InitWithShallowSlice(vel_l_,2,tid,1);
    vel_r.InitWithShallowSlice(vel_r_,2,tid,1);
    vdiff_l.InitWithShallowSlice(vdiff_l_,2,tid,1);
    vdiff_r.InitWithShallowSlice(vdiff_r_,2,tid,1);
    eddl.InitWithShallowSlice(eddl_,3,tid,1);
    eddr.InitWithShallowSlice(eddr_,3,tid,1);
    wl.InitWithShallowSlice(wl_,3,tid,1);
    wr.InitWithShallowSlice(wr_,3,tid,1);
    cwidth.InitWithShallowSlice(cwidth_,2,tid,1);

    // The area functions needed to calculate Grad Pc
    AthenaArray<Real> x1area, x2area, x2area_p1, x3area, x3area_p1, vol, dflx;
    x1area.InitWithShallowSlice(x1face_area_,2,tid,1);
    x2area.InitWithShallowSlice(x2face_area_,2,tid,1);
    x2area_p1.InitWithShallowSlice(x2face_area_p1_,2,tid,1);
    x3area.InitWithShallowSlice(x3face_area_,2,tid,1);
    x3area_p1.InitWithShallowSlice(x3face_area_p1_,2,tid,1);
    vol.InitWithShallowSlice(cell_volume_,2,tid,1);

//--------------------------------------------------------------------------------------
    // First, calculate the diffusion velocity along three coordinate system
    for (int k=0; k<n3; ++k){
      for(int j=0; j<n2; ++j){

        // diffusion velocity along the direction of sigma vector
        // We first assume B is along x coordinate
        // Then rotate according to B direction to the actual acooridnate

        for(int i=0; i<n1; ++i){
          Real eddxx=pcr->prtensor_cr(PC11,k,j,i);
          Real totsigma = pcr->sigma_diff(0,k,j,i);
          if(pcr->stream_flag)
            totsigma = 1.0/(1.0/pcr->sigma_diff(0,k,j,i) 
                          + 1.0/pcr->sigma_adv(0,k,j,i));
          Real taux = taufact_ * totsigma * pco->dx1f(i);
          taux = taux * taux/(2.0 * eddxx);
          Real diffv = sqrt((1.0 - exp(-taux)) / taux);

          if(taux < 1.e-3)
            diffv = sqrt((1.0 - 0.5* taux));

          pcr->v_diff(0,k,j,i) = pcr->vmax * sqrt(eddxx) * diffv;
        }

        // y direction
        pco->CenterWidth2(k,j,0,n1-1,cwidth);
          // get the optical depth across the cell
        for(int i=0; i<n1; ++i){
          Real eddyy=pcr->prtensor_cr(PC22,k,j,i);
          Real totsigma = pcr->sigma_diff(1,k,j,i);
          if(pcr->stream_flag)
            totsigma = 1.0/(1.0/pcr->sigma_diff(1,k,j,i) 
                          + 1.0/pcr->sigma_adv(1,k,j,i));
          Real tauy = taufact_ * totsigma * cwidth(i);
          tauy = tauy * tauy/(2.0 * eddyy);
          Real diffv = sqrt((1.0 - exp(-tauy)) / tauy);

          if(tauy < 1.e-3)
            diffv = sqrt((1.0 - 0.5* tauy));

          pcr->v_diff(1,k,j,i) = pcr->vmax * sqrt(eddyy) * diffv;            
        }// end i
        // z direction
        pco->CenterWidth3(k,j,0,n1-1,cwidth);
          // get the optical depth across the cell
        for(int i=0; i<n1; ++i){
          Real eddzz=pcr->prtensor_cr(PC33,k,j,i);
          Real totsigma = pcr->sigma_diff(2,k,j,i);
          if(pcr->stream_flag)
            totsigma = 1.0/(1.0/pcr->sigma_diff(2,k,j,i) 
                          + 1.0/pcr->sigma_adv(2,k,j,i));
          Real tauz = taufact_ * totsigma * cwidth(i);
          tauz = tauz * tauz/(2.0 * eddzz);
          Real diffv = sqrt((1.0 - exp(-tauz)) / tauz);

          if(tauz < 1.e-3)
            diffv = sqrt((1.0 - 0.5* tauz));

          pcr->v_diff(2,k,j,i) = pcr->vmax * sqrt(eddzz) * diffv;            
        }

        //rotate the v_diff vector to the local coordinate
        if(MAGNETIC_FIELDS_ENABLED){
          for(int i=0; i<n1; ++i){
            

            InvRotateVec(pcr->b_angle(0,k,j,i),pcr->b_angle(1,k,j,i),
                        pcr->b_angle(2,k,j,i),pcr->b_angle(3,k,j,i), 
                          pcr->v_diff(0,k,j,i),pcr->v_diff(1,k,j,i),
                                              pcr->v_diff(2,k,j,i));
            // take the absolute value
            // Also add the Alfven velocity for the streaming flux
            pcr->v_diff(0,k,j,i) = fabs(pcr->v_diff(0,k,j,i));

            pcr->v_diff(1,k,j,i) = fabs(pcr->v_diff(1,k,j,i));
                                   
            pcr->v_diff(2,k,j,i) = fabs(pcr->v_diff(2,k,j,i));

          }

        }
        // need to add additional sound speed for stability
        for(int i=0; i<n1; ++i){
           Real cr_sound_x = vel_flx_flag_ * sqrt((4.0/3.0) * pcr->prtensor_cr(PC11,k,j,i) 
                                    * u_cr(k,j,i)/w(IDN,k,j,i)); 

           pcr->v_diff(0,k,j,i) += cr_sound_x;

           Real cr_sound_y = vel_flx_flag_ * sqrt((4.0/3.0) * pcr->prtensor_cr(PC22,k,j,i) 
                                    * u_cr(k,j,i)/w(IDN,k,j,i));

           pcr->v_diff(1,k,j,i) += cr_sound_y;

           Real cr_sound_z = vel_flx_flag_ * sqrt((4.0/3.0) * pcr->prtensor_cr(PC33,k,j,i) 
                                    * u_cr(k,j,i)/w(IDN,k,j,i)); 
           
           pcr->v_diff(2,k,j,i) += cr_sound_z;

        }


      }
    }
     
//--------------------------------------------------------------------------------------
// i-direction
    // set the loop limits
    jl=js, ju=je, kl=ks, ku=ke;
    for (int k=kl; k<=ku; ++k){
#pragma omp for schedule(static)
      for (int j=jl; j<=ju; ++j){
        // First, need to do reconstruction
        // to reconstruct Ec, Fc, vel, v_a and 
        // return Ec,Fc and signal speed at left and right state
        if(reconstruct_order == 1){
          DonorCellX1(k,j,is,ie+1,u_cr,w,pcr->prtensor_cr,wl,wr,
                                        vel_l,vel_r, eddl, eddr);
        }else {
          PieceWiseLinear(k,j,is,ie+1,u_cr,w,pcr->prtensor_cr,wl,wr,
                                     vel_l, vel_r, eddl, eddr, 0);
        }
        // get the optical depth across the cell
#pragma omp simd
        for(int i=is; i<=ie+1; ++i){
          vdiff_l(i) = pcr->v_diff(0,k,j,i-1);
          vdiff_r(i) = pcr->v_diff(0,k,j,i);
        }

        // calculate the flux
        CRFlux(CRF1, k, j, is, ie+1, wl, wr, vel_l, vel_r, eddl, eddr,  
                                               vdiff_l, vdiff_r, flx);
        // store the flux
        for(int n=0; n<NCR; ++n){
#pragma omp simd
          for(int i=is; i<=ie+1; ++i){
            x1flux(n,k,j,i) = flx(n,i);
          }
        }

      }
    }
    
// j-direction
    if(pmb->block_size.nx2 > 1){
      il=is; iu=ie; kl=ks; ku=ke;
      for (int k=kl; k<=ku; ++k){
#pragma omp for schedule(static)
        for (int j=js; j<=je+1; ++j){

          if(reconstruct_order == 1){
            DonorCellX2(k,j,il,iu,u_cr,w,pcr->prtensor_cr,wl,wr,
                                     vel_l,vel_r, eddl, eddr);
          }else {
            PieceWiseLinear(k,j,il,iu,u_cr,w,pcr->prtensor_cr,wl,wr,
                                     vel_l, vel_r, eddl, eddr, 1);
          }

          // get the optical depth across the cell
#pragma omp simd
          for(int i=il; i<=iu; ++i){
            vdiff_l(i) = pcr->v_diff(1,k,j-1,i);
            vdiff_r(i) = pcr->v_diff(1,k,j,i);
          }

         // calculate the flux
          CRFlux(CRF2, k, j, il, iu, wl, wr, vel_l, vel_r, eddl, eddr, 
                                                vdiff_l, vdiff_r, flx);
        
        // store the flux
          for(int n=0; n<NCR; ++n){
#pragma omp simd
            for(int i=is; i<=ie; ++i){
              x2flux(n,k,j,i) = flx(n,i);
            }
          }

        }
      }
    }// finish j direction

//  k-direction
    if(pmb->block_size.nx3 > 1){
      il=is; iu=ie; jl=js; ju=je;
      for (int k=ks; k<=ke+1; ++k){
#pragma omp for schedule(static)
        for (int j=jl; j<=ju; ++j){

          if(reconstruct_order == 1){
            DonorCellX3(k,j,il,iu,u_cr,w,pcr->prtensor_cr,wl,wr,
                                        vel_l,vel_r, eddl, eddr);
          }else {
            PieceWiseLinear(k,j,il,iu,u_cr,w,pcr->prtensor_cr,wl,wr,
                                       vel_l, vel_r, eddl, eddr, 2);
          }

          // get the optical depth across the cell
#pragma omp simd
          for(int i=il; i<=iu; ++i){
            vdiff_l(i) = pcr->v_diff(2,k-1,j,i);
            vdiff_r(i) = pcr->v_diff(2,k,j,i);
          }
         // calculate the flux
          CRFlux(CRF3, k, j, il, iu, wl, wr, vel_l, vel_r, eddl, eddr,
                                               vdiff_l, vdiff_r, flx);
        
        // store the flux
          for(int n=0; n<NCR; ++n){
#pragma omp simd
            for(int i=is; i<=ie; ++i){
              x3flux(n,k,j,i) = flx(n,i);
            }
          }

        }
      }
    }// finish k direction

  
// Now calculate Grad Pc and the associated heating term
#pragma omp for schedule(static)
    for (int k=ks; k<=ke; ++k) { 
      for (int j=js; j<=je; ++j) {
        pmb->pcoord->Face1Area(k,j,is,ie+1,x1area);
        pmb->pcoord->CellVolume(k,j,is,ie,vol);
        // x1 direction
        for(int n=0; n<3; ++n){
#pragma omp simd
          for(int i=is; i<=ie; ++i){
            grad_pc_(n,k,j,i) = (x1area(i+1)*x1flux(CRF1+n,k,j,i+1) 
                               - x1area(i)  *x1flux(CRF1+n,k,j,i))/vol(i);
          }
        } 

        if(pmb->block_size.nx2 > 1){
          pmb->pcoord->Face2Area(k,j  ,is,ie,x2area   );
          pmb->pcoord->Face2Area(k,j+1,is,ie,x2area_p1);
          for(int n=0; n<3; ++n){
#pragma omp simd
            for(int i=is; i<=ie; ++i){
              grad_pc_(n,k,j,i) += (x2area_p1(i)*x2flux(CRF1+n,k,j+1,i) 
                                 -  x2area(i)  *x2flux(CRF1+n,k,j,i))/vol(i);
            }// end i
          }  
        }// end nx2

        if(pmb->block_size.nx3 > 1){
          pmb->pcoord->Face3Area(k  ,j,is,ie,x3area   );
          pmb->pcoord->Face3Area(k+1,j,is,ie,x3area_p1);
          for(int n=0; n<3; ++n){
#pragma omp simd
            for(int i=is; i<=ie; ++i){
              grad_pc_(n,k,j,i) += (x3area_p1(i) *x3flux(CRF1+n,k+1,j,i) 
                                  - x3area(i)*x3flux(CRF1+n,k,j,i))/vol(i);
            } 
          } 
        }// end nx3


        for(int n=0; n<3; ++n){
#pragma omp simd
          for(int i=is; i<=ie; ++i){
            grad_pc_(n,k,j,i) *= invlim;
          }
        } 



       // calculate streaming velocity with magnetic field
        for(int i=is; i<=ie; ++i){
            Real vtotx = w(IVX,k,j,i) + pcr->v_adv(0,k,j,i);
            Real vtoty = w(IVY,k,j,i) + pcr->v_adv(1,k,j,i);
            Real vtotz = w(IVZ,k,j,i) + pcr->v_adv(2,k,j,i);
            Real v_dot_gradpc = vtotx * grad_pc_(0,k,j,i) 
                              + vtoty * grad_pc_(1,k,j,i) 
                              + vtotz * grad_pc_(2,k,j,i);


            if(pcr->stream_flag){
              if(MAGNETIC_FIELDS_ENABLED){

                Real inv_sqrt_rho = 1.0/sqrt(w(IDN,k,j,i));

                Real pb= bcc(IB1,k,j,i)*bcc(IB1,k,j,i)
                        +bcc(IB2,k,j,i)*bcc(IB2,k,j,i)
                        +bcc(IB3,k,j,i)*bcc(IB3,k,j,i);

                Real b_grad_pc = bcc(IB1,k,j,i) * grad_pc_(0,k,j,i) 
                               + bcc(IB2,k,j,i) * grad_pc_(1,k,j,i) 
                               + bcc(IB3,k,j,i) * grad_pc_(2,k,j,i);

                Real va1 = bcc(IB1,k,j,i) * inv_sqrt_rho;
                Real va2 = bcc(IB2,k,j,i) * inv_sqrt_rho;
                Real va3 = bcc(IB3,k,j,i) * inv_sqrt_rho;

                Real va = sqrt(pb) * inv_sqrt_rho;
                Real dpc_sign = 0.0;

                if(b_grad_pc > TINY_NUMBER) dpc_sign = 1.0;
                else if(-b_grad_pc > TINY_NUMBER) dpc_sign = -1.0;

                pcr->v_adv(0,k,j,i) = -va1 * dpc_sign;
                pcr->v_adv(1,k,j,i) = -va2 * dpc_sign;
                pcr->v_adv(2,k,j,i) = -va3 * dpc_sign;

                if(va > TINY_NUMBER){
                  pcr->sigma_adv(0,k,j,i) = fabs(b_grad_pc)/(sqrt(pb) * va * 
                                 (1.0 + pcr->prtensor_cr(PC11,k,j,i)) 
                                            * invlim * u_cr(CRE,k,j,i));
                  pcr->sigma_adv(1,k,j,i) = pcr->max_opacity;
                  pcr->sigma_adv(2,k,j,i) = pcr->max_opacity;
                }

                Real sigma_x = 1.0/(1.0/pcr->sigma_diff(0,k,j,i) + 
                               1.0/pcr->sigma_adv(0,k,j,i));

                vtotx = w(IVX,k,j,i) + pcr->v_adv(0,k,j,i);
                vtoty = w(IVY,k,j,i) + pcr->v_adv(1,k,j,i);
                vtotz = w(IVZ,k,j,i) + pcr->v_adv(2,k,j,i);

                Real v1 = w(IVX,k,j,i);
                Real v2 = w(IVY,k,j,i);
                Real v3 = w(IVZ,k,j,i);

                Real dpcdx = grad_pc_(0,k,j,i);
                Real dpcdy = grad_pc_(1,k,j,i);
                Real dpcdz = grad_pc_(2,k,j,i);

               // explicity method needs to use CR flux from previous step
                Real fr1 = u_cr(CRF1,k,j,i);
                Real fr2 = u_cr(CRF2,k,j,i);
                Real fr3 = u_cr(CRF3,k,j,i);

               // perform rotation
                RotateVec(pcr->b_angle(0,k,j,i),pcr->b_angle(1,k,j,i),
                         pcr->b_angle(2,k,j,i),pcr->b_angle(3,k,j,i),vtotx,vtoty,vtotz);

                RotateVec(pcr->b_angle(0,k,j,i),pcr->b_angle(1,k,j,i),
                         pcr->b_angle(2,k,j,i),pcr->b_angle(3,k,j,i),dpcdx,dpcdy,dpcdz);

                RotateVec(pcr->b_angle(0,k,j,i),pcr->b_angle(1,k,j,i),
                         pcr->b_angle(2,k,j,i),pcr->b_angle(3,k,j,i),v1,v2,v3);


                RotateVec(pcr->b_angle(0,k,j,i),pcr->b_angle(1,k,j,i),
                         pcr->b_angle(2,k,j,i),pcr->b_angle(3,k,j,i),fr1,fr2,fr3);


              // only calculate v_dot_gradpc perpendicular to B
              // perpendicular direction only has flow velocity, no streaming velocity
                v_dot_gradpc = v2 * dpcdy + v3 * dpcdz;

                Real fxx = pcr->prtensor_cr(PC11,k,j,i);
                Real fxy = pcr->prtensor_cr(PC12,k,j,i);
                Real fxz = pcr->prtensor_cr(PC13,k,j,i);
              
                Real fr_cm1 = fr1 - (v1 * (1.0 + fxx) + v2 * fxy + v3 * fxz) 
                              * u_cr(CRE,k,j,i) * invlim;

                ec_source_(k,j,i) = (v_dot_gradpc - vtotx * sigma_x * fr_cm1);

             }// end mhd
             else{


                Real v1 = w(IVX,k,j,i);
                Real v2 = w(IVY,k,j,i);
                Real v3 = w(IVZ,k,j,i);

                Real dpcdx = grad_pc_(0,k,j,i);
                Real dpcdy = grad_pc_(1,k,j,i);
                Real dpcdz = grad_pc_(2,k,j,i);

                Real fr1 = u_cr(CRF1,k,j,i);
                Real fr2 = u_cr(CRF2,k,j,i);
                Real fr3 = u_cr(CRF3,k,j,i);


                Real fxx = pcr->prtensor_cr(PC11,k,j,i);
                Real fxy = pcr->prtensor_cr(PC12,k,j,i);
                Real fxz = pcr->prtensor_cr(PC13,k,j,i);
                Real fyy = pcr->prtensor_cr(PC22,k,j,i);
                Real fyz = pcr->prtensor_cr(PC23,k,j,i);
                Real fzz = pcr->prtensor_cr(PC33,k,j,i);
              
                Real fr_cm1 = u_cr(CRF1,k,j,i) - (v1 * (1.0 + fxx) 
                            + v2 * fxy + v3 * fxz) * u_cr(CRE,k,j,i) * invlim; 
                Real fr_cm2 = u_cr(CRF2,k,j,i) - (v2 * (1.0 + fyy) 
                            + v1 * fxy + v3 * fyz) * u_cr(CRE,k,j,i) * invlim;                
                Real fr_cm3 = u_cr(CRF3,k,j,i) - (v3 * (1.0 + fzz) 
                            + v2 * fyz + v1 * fxz) * u_cr(CRE,k,j,i) * invlim; 

                Real sigma_x = 1.0/(1.0/pcr->sigma_diff(0,k,j,i) + 
                               1.0/pcr->sigma_adv(0,k,j,i));

                Real sigma_y = 1.0/(1.0/pcr->sigma_diff(1,k,j,i) + 
                               1.0/pcr->sigma_adv(1,k,j,i));

                Real sigma_z = 1.0/(1.0/pcr->sigma_diff(2,k,j,i) + 
                               1.0/pcr->sigma_adv(2,k,j,i));

                ec_source_(k,j,i) = -(vtotx * sigma_x * fr_cm1 + vtoty * sigma_y * fr_cm2
                                      + vtotz * sigma_z * fr_cm3);                

             }// no magnetic field
           }// if stream flag is turned on. Otherwise, use default values

         }// end i

      }// end k
    }// end j


}// end of omp parallel region
  

  
}


void CRIntegrator::FluxDivergence(MeshBlock *pmb, AthenaArray<Real> &u_cr1,
                      AthenaArray<Real> &u_cr2, const IntegratorWeight wght,
                      AthenaArray<Real> &u_out, AthenaArray<Real> &u, 
                      AthenaArray<Real> &w, AthenaArray<Real> &bcc)
{
  CosmicRay *pcr=pmb->pcr;
  Coordinates *pco = pmb->pcoord;

  AthenaArray<Real> &x1flux=pcr->flux[X1DIR];
  AthenaArray<Real> &x2flux=pcr->flux[X2DIR];
  AthenaArray<Real> &x3flux=pcr->flux[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  
  Real dt = pmb->pmy_mesh->dt;

  Real invlim = 1.0/pcr->vmax;

    
  int tid=0;
  int nthreads = pmb->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) private(tid) num_threads(nthreads)
{
#ifdef OPENMP_PARALLEL
  tid=omp_get_thread_num();
#endif
  AthenaArray<Real> x1area, x2area, x2area_p1, x3area, x3area_p1, vol, dflx;
  x1area.InitWithShallowSlice(x1face_area_,2,tid,1);
  x2area.InitWithShallowSlice(x2face_area_,2,tid,1);
  x2area_p1.InitWithShallowSlice(x2face_area_p1_,2,tid,1);
  x3area.InitWithShallowSlice(x3face_area_,2,tid,1);
  x3area_p1.InitWithShallowSlice(x3face_area_p1_,2,tid,1);
  vol.InitWithShallowSlice(cell_volume_,2,tid,1);
  dflx.InitWithShallowSlice(flx_,3,tid,1);
    

#pragma omp for schedule(static)
  for (int k=ks; k<=ke; ++k) { 
    for (int j=js; j<=je; ++j) {

      // calculate x1-flux divergence 
      pmb->pcoord->Face1Area(k,j,is,ie+1,x1area);
      for(int n=0; n<NCR; ++n){
#pragma omp simd
        for(int i=is; i<=ie; ++i){
          dflx(n,i) = (x1area(i+1) *x1flux(n,k,j,i+1) - x1area(i)*x1flux(n,k,j,i));
        }// end i
      }// End n


     // calculate x2-flux
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Face2Area(k,j  ,is,ie,x2area   );
        pmb->pcoord->Face2Area(k,j+1,is,ie,x2area_p1);
      for(int n=0; n<NCR; ++n){
#pragma omp simd
        for(int i=is; i<=ie; ++i){
            dflx(n,i) += (x2area_p1(i)*x2flux(n,k,j+1,i) - x2area(i)*x2flux(n,k,j,i));
          }
        }
      }// end nx2

      // calculate x3-flux divergence
      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Face3Area(k  ,j,is,ie,x3area   );
        pmb->pcoord->Face3Area(k+1,j,is,ie,x3area_p1);
        for(int n=0; n<NCR; ++n){
#pragma omp simd
          for(int i=is; i<=ie; ++i){

            dflx(n,i) += (x3area_p1(i)*x3flux(n,k+1,j,i) - x3area(i)*x3flux(n,k,j,i));
          }
        }
      }// end nx3
      // update variable with flux divergence
      pmb->pcoord->CellVolume(k,j,is,ie,vol);
      for(int n=0; n<NCR; ++n){
#pragma omp simd
        for(int i=is; i<=ie; ++i){
          u_out(n,k,j,i) = wght.a*u_cr1(n,k,j,i) + wght.b*u_cr2(n,k,j,i)
                       - wght.c*dt*dflx(n,i)/vol(i);
        }
      }
      // check cosmic ray energy density is always positive
      for(int i=is; i<=ie; ++i){
        if(u_out(CRE,k,j,i) < TINY_NUMBER)
          u_out(CRE,k,j,i) = TINY_NUMBER;
      }    

    }// end j
  }// End k
  


}// end omp parallel region

}
