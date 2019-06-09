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
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../parameter_input.hpp"
#include "../../../mesh/mesh.hpp"
#include "../tc.hpp"
#include "../../../coordinates/coordinates.hpp" //
#include "../../hydro.hpp"
#include <algorithm>   // min,max

// class header
#include "tc_integrators.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif


void NTCIntegrator::CalculateFluxes(MeshBlock *pmb,
      AthenaArray<Real> &w, AthenaArray<Real> &bcc, 
      AthenaArray<Real> &u_tc, int reconstruct_order)
{
  NewThermalConduction *ptc=pmy_tc;
  Coordinates *pco = pmb->pcoord;
  
  Real invlim = 1.0/ptc->vmax;
  
  AthenaArray<Real> &x1flux=ptc->flux[X1DIR];
  AthenaArray<Real> &x2flux=ptc->flux[X2DIR];
  AthenaArray<Real> &x3flux=ptc->flux[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;


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



    AthenaArray<Real> flx, wl, wr, vdiff_l, vdiff_r;
    AthenaArray<Real> rhol, rhor, tl, tr;
    AthenaArray<Real> cwidth;

    flx.InitWithShallowSlice(flx_,3,tid,1);

    vdiff_l.InitWithShallowSlice(vdiff_l_,2,tid,1);
    vdiff_r.InitWithShallowSlice(vdiff_r_,2,tid,1);


    wl.InitWithShallowSlice(wl_,3,tid,1);
    wr.InitWithShallowSlice(wr_,3,tid,1);

    rhol.InitWithShallowSlice(rho_l_,2,tid,1);
    rhor.InitWithShallowSlice(rho_r_,2,tid,1);
    tl.InitWithShallowSlice(tgas_l_,2,tid,1);
    tr.InitWithShallowSlice(tgas_r_,2,tid,1);

    cwidth.InitWithShallowSlice(cwidth_,2,tid,1);



//--------------------------------------------------------------------------------------
    // First, calculate the diffusion velocity along three coordinate system
    for (int k=0; k<n3; ++k){
      for(int j=0; j<n2; ++j){

        // diffusion velocity along the direction of sigma vector
        // We first assume B is along x coordinate
        // Then rotate according to B direction to the actual acooridnate

        for(int i=0; i<n1; ++i){
          Real kappa = ptc->rho(k,j,i)/ptc->ntc_kappa(0,k,j,i);

          Real taux = taufact_ * ptc->vmax * kappa * pco->dx1f(i);
          taux = taux * taux/2.0;
          Real diffv = sqrt((1.0 - exp(-taux)) / taux);

          if(taux < 1.e-3)
            diffv = sqrt((1.0 - 0.5* taux));

          vdiff_(0,k,j,i) = ptc->vmax * diffv;
        }

        // y direction
        pco->CenterWidth2(k,j,0,n1-1,cwidth);
          // get the optical depth across the cell
        for(int i=0; i<n1; ++i){
          Real kappa = ptc->rho(k,j,i)/ptc->ntc_kappa(1,k,j,i);
          Real tauy = taufact_ * ptc->vmax * kappa * cwidth(i);          
          tauy = tauy * tauy/2.0;
          Real diffv = sqrt((1.0 - exp(-tauy)) / tauy);

          if(tauy < 1.e-3)
            diffv = sqrt((1.0 - 0.5* tauy));

          vdiff_(1,k,j,i) = ptc->vmax * diffv;            
        }// end i
        // z direction
        pco->CenterWidth3(k,j,0,n1-1,cwidth);
          // get the optical depth across the cell
        for(int i=0; i<n1; ++i){
          Real kappa = ptc->rho(k,j,i)/ptc->ntc_kappa(2,k,j,i);
          Real tauz = taufact_ * ptc->vmax * kappa * cwidth(i);  
          tauz = tauz * tauz/2.0;
          Real diffv = sqrt((1.0 - exp(-tauz)) / tauz);

          if(tauz < 1.e-3)
            diffv = sqrt((1.0 - 0.5* tauz));

          vdiff_(2,k,j,i) = ptc->vmax * diffv;            
        }

        //rotate the v_diff vector to the local coordinate
        if(MAGNETIC_FIELDS_ENABLED){
          for(int i=0; i<n1; ++i){
            

            InvRotateVec(ptc->b_angle(0,k,j,i),ptc->b_angle(1,k,j,i),
                        ptc->b_angle(2,k,j,i),ptc->b_angle(3,k,j,i), 
                      vdiff_(0,k,j,i),vdiff_(1,k,j,i),vdiff_(2,k,j,i));
            // take the absolute value
            // Also add the Alfven velocity for the streaming flux
            vdiff_(0,k,j,i) = fabs(vdiff_(0,k,j,i));

            vdiff_(1,k,j,i) = fabs(vdiff_(1,k,j,i));
                                   
            vdiff_(2,k,j,i) = fabs(vdiff_(2,k,j,i));

          }

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
        // to reconstruct Ec, Fc, vel, and 
        // return Ec,Fc and signal speed at left and right state
        if(reconstruct_order == 1){
          DonorCellX1(k,j,is,ie+1,u_tc,ptc->rho,ptc->tgas,
                                    rhol,rhor,tl,tr,wl,wr);
        }else {
          PieceWiseLinear(k,j,is,ie+1,u_tc,ptc->rho,ptc->tgas,
                                  rhol,rhor,tl,tr,wl,wr,0);
        }
        // get the optical depth across the cell
#pragma omp simd
        for(int i=is; i<=ie+1; ++i){
          vdiff_l(i) = -vdiff_(0,k,j,i-1);
          vdiff_r(i) = vdiff_(0,k,j,i);
        }

        // calculate the flux
        NTCFlux(1, is, ie+1, tl, tr, rhol,rhor, wl, wr, 
                                  vdiff_l, vdiff_r, flx);
        // store the flux
        for(int n=0; n<4; ++n){
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
            DonorCellX2(k,j,il,iu,u_tc,ptc->rho,ptc->tgas, 
                                     rhol,rhor,tl,tr,wl,wr);
          }else {
            PieceWiseLinear(k,j,il,iu,u_tc,ptc->rho,ptc->tgas,
                                    rhol,rhor,tl,tr,wl,wr,1);
          }

          // get the optical depth across the cell
#pragma omp simd
          for(int i=il; i<=iu; ++i){
            vdiff_l(i) = -vdiff_(1,k,j-1,i);
            vdiff_r(i) = vdiff_(1,k,j,i);
          }

         // calculate the flux
          NTCFlux(2, il, iu, tl, tr, rhol, rhor, wl,wr, 
                                  vdiff_l, vdiff_r, flx);
        
        // store the flux
          for(int n=0; n<4; ++n){
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
            DonorCellX3(k,j,il,iu,u_tc,ptc->rho,ptc->tgas,
                                      rhol,rhor,tl,tr,wl,wr);
          }else {
            PieceWiseLinear(k,j,il,iu,u_tc,ptc->rho,ptc->tgas,
                                   rhol,rhor,tl,tr,wl,wr,2);
          }

          // get the optical depth across the cell
#pragma omp simd
          for(int i=il; i<=iu; ++i){
            vdiff_l(i) = -vdiff_(2,k-1,j,i); 
            vdiff_r(i) = vdiff_(2,k,j,i);
          }
         // calculate the flux
          NTCFlux(3, il, iu, tl, tr, rhol, rhor, wl, wr, 
                                    vdiff_l, vdiff_r, flx);
        
        // store the flux
          for(int n=0; n<4; ++n){
#pragma omp simd
            for(int i=is; i<=ie; ++i){
              x3flux(n,k,j,i) = flx(n,i);
            }
          }

        }
      }
    }// finish k direction

}// end of omp parallel region
  
  
}


void NTCIntegrator::FluxDivergence(MeshBlock *pmb, AthenaArray<Real> &w,
     AthenaArray<Real> &u_tc1, AthenaArray<Real> &u_tc2, 
     const IntegratorWeight wght, AthenaArray<Real> &u_out)
{
  NewThermalConduction *ptc=pmy_tc;
  Coordinates *pco = pmb->pcoord;

  AthenaArray<Real> &x1flux=ptc->flux[X1DIR];
  AthenaArray<Real> &x2flux=ptc->flux[X2DIR];
  AthenaArray<Real> &x3flux=ptc->flux[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  
  Real dt = pmb->pmy_mesh->dt;

    
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
      for(int n=0; n<4; ++n){
#pragma omp simd
        for(int i=is; i<=ie; ++i){
          dflx(n,i) = (x1area(i+1) *x1flux(n,k,j,i+1) - x1area(i)*x1flux(n,k,j,i));
        }// end i
      }// End n


     // calculate x2-flux
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Face2Area(k,j  ,is,ie,x2area   );
        pmb->pcoord->Face2Area(k,j+1,is,ie,x2area_p1);
      for(int n=0; n<4; ++n){
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
        for(int n=0; n<4; ++n){
#pragma omp simd
          for(int i=is; i<=ie; ++i){

            dflx(n,i) += (x3area_p1(i)*x3flux(n,k+1,j,i) - x3area(i)*x3flux(n,k,j,i));
          }
        }
      }// end nx3
      // update variable with flux divergence
      pmb->pcoord->CellVolume(k,j,is,ie,vol);
      for(int n=1; n<4; ++n){
#pragma omp simd
        for(int i=is; i<=ie; ++i){
          u_out(n,k,j,i) = wght.a*u_tc1(n,k,j,i) + wght.b*u_tc2(n,k,j,i)
                       - wght.c*dt*ptc->rho(k,j,i)*dflx(n,i)/vol(i);
        }
      }

      //get the energy source term
#pragma omp simd
      for(int i=is; i<=ie; ++i){
        tc_esource_(k,j,i) = -wght.c*dt*dflx(0,i)/vol(i);
      }      

    }// end j
  }// End k
  


}// end omp parallel region

}
