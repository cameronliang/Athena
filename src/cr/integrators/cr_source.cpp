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
//! \file rad_source.cpp
//  \brief Add radiation source terms to both radiation and gas
//======================================================================================


// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"
#include "../cr.hpp"
#include "../../coordinates/coordinates.hpp" //
#include "../../hydro/hydro.hpp"
#include "../../field/field.hpp"
#include "../../eos/eos.hpp"

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

//add the source terms implicitly



void CRIntegrator::AddSourceTerms(MeshBlock *pmb, const Real dt, AthenaArray<Real> &u,
        AthenaArray<Real> &w, AthenaArray<Real> &bcc,
        AthenaArray<Real> &u_cr, const int step)
{
  CosmicRay *pcr=pmb->pcr;
  Coordinates *pco = pmb->pcoord;

  Real vlim = pcr->vmax;
  Real invlim = 1.0/vlim;

// The information stored in the array
// b_angle is
// b_angle[0]=sin_theta_b
// b_angle[1]=cos_theta_b
// b_angle[2]=sin_phi_b
// b_angle[3]=cos_phi_b
  
    
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
 
  for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){

         Real *vel1 = &(w(IVX,k,j,0));
         Real *vel2 = &(w(IVY,k,j,0));
         Real *vel3 = &(w(IVZ,k,j,0));

         Real *pc11 = &(pcr->prtensor_cr(PC11,k,j,0));
         Real *pc22 = &(pcr->prtensor_cr(PC22,k,j,0));
         Real *pc33 = &(pcr->prtensor_cr(PC33,k,j,0));
         Real *pc12 = &(pcr->prtensor_cr(PC12,k,j,0));
         Real *pc13 = &(pcr->prtensor_cr(PC13,k,j,0));
         Real *pc23 = &(pcr->prtensor_cr(PC23,k,j,0));

         Real *ec = &(u_cr(CRE,k,j,0));
         Real *fc1 = &(u_cr(CRF1,k,j,0));
         Real *fc2 = &(u_cr(CRF2,k,j,0));
         Real *fc3 = &(u_cr(CRF3,k,j,0));

         // The angle of B
         Real *sint_b = &(pcr->b_angle(0,k,j,0));
         Real *cost_b = &(pcr->b_angle(1,k,j,0));
         Real *sinp_b = &(pcr->b_angle(2,k,j,0));
         Real *cosp_b = &(pcr->b_angle(3,k,j,0));

 // adv1 is dPc/dx, adv2 is dPc/dy, adv3 is dPc/dz

      for(int i=is; i<=ie; ++i){

         Real v1 = vel1[i];
         Real v2 = vel2[i];
         Real v3 = vel3[i];

         Real fr1 = fc1[i];
         Real fr2 = fc2[i];
         Real fr3 = fc3[i];



         Real fxx=pc11[i], fyy=pc22[i], fzz=pc33[i], fxy=pc12[i], 
              fyz=pc23[i], fxz=pc13[i];

         // in the case with magnetic field
        // rotate the vectors to oriante to the B direction
         if(MAGNETIC_FIELDS_ENABLED){
           // Apply rotation of the vectors
           RotateVec(sint_b[i],cost_b[i],sinp_b[i],cosp_b[i],v1,v2,v3);

           RotateVec(sint_b[i],cost_b[i],sinp_b[i],cosp_b[i],fr1,fr2,fr3);

         }

         Real sigma_x = pcr->sigma_diff(0,k,j,i);
         Real sigma_y = pcr->sigma_diff(1,k,j,i);
         Real sigma_z = pcr->sigma_diff(2,k,j,i);

         if(pcr->stream_flag){
           sigma_x = 1.0/(1.0/pcr->sigma_diff(0,k,j,i) + 
                             1.0/pcr->sigma_adv(0,k,j,i));

           sigma_y = 1.0/(1.0/pcr->sigma_diff(1,k,j,i) +
                             1.0/pcr->sigma_adv(1,k,j,i));

           sigma_z = 1.0/(1.0/pcr->sigma_diff(2,k,j,i) + 
                             1.0/pcr->sigma_adv(2,k,j,i));
         }



         // Now update the momentum equation
         //\partial F/\partial t=-V_m\sigma (F-v(E+Pc_/v_m)) 


         Real dtsigma1 = dt * sigma_x * vlim;
         Real dtsigma2 = dt * sigma_y * vlim;
         Real dtsigma3 = dt * sigma_z * vlim;

         Real rhs1 = (v1 * (1.0 + fxx) + v2 * fxy + v3 * fxz)
                            * invlim * ec[i];
         Real rhs2 = (v2 * (1.0 + fyy) + v1 * fxy + v3 * fyz)
                            * invlim * ec[i];
         Real rhs3 = (v3 * (1.0 + fzz) + v1 * fxz + v2 * fyz)
                            * invlim * ec[i];

         Real newfr1 = (fr1 - rhs1) / (1.0 + dtsigma1) + rhs1;
         Real newfr2 = (fr2 - rhs2) / (1.0 + dtsigma2) + rhs2;
         Real newfr3 = (fr3 - rhs3) / (1.0 + dtsigma3) + rhs3;
   

        // Now apply the invert rotation
         if(MAGNETIC_FIELDS_ENABLED){
           // Apply rotation of the vectors
           InvRotateVec(sint_b[i],cost_b[i],sinp_b[i],cosp_b[i], 
                                           newfr1,newfr2,newfr3);
         }        

         // Add the energy source term
         Real new_ec = ec[i] + dt * ec_source_(k,j,i);
         Real new_eg = u(IEN,k,j,i) - dt * ec_source_(k,j,i);

         if(new_ec > 0.0 && new_eg > 0.0){

            new_sol_(CRE,i) = new_ec;
            u(IEN,k,j,i) = new_eg;
         }

         new_sol_(CRF1,i) = newfr1;  
         new_sol_(CRF2,i) = newfr2;
         new_sol_(CRF3,i) = newfr3;


         u(IM1,k,j,i) += (-(newfr1 - fc1[i]) * invlim);
         u(IM2,k,j,i) += (-(newfr2 - fc2[i]) * invlim);
         u(IM3,k,j,i) += (-(newfr3 - fc3[i]) * invlim);

         
      }// end i
        

      // update Cosmic Ray quantities
      for(int n=0; n<NCR; ++n){
        for(int i=is; i<=ie; ++i){
           u_cr(n,k,j,i) = new_sol_(n,i);
        }
      }



    }// end j
  }// end k
      
}


