//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

// C++ headers
#include <iostream>   // endl
#include <fstream>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>      // sqrt
#include <algorithm>  // min

// Athena++ headers
#include "../globals.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../bvals/bvals.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"
#include "../cr/cr.hpp"
#include "../cr/integrators/cr_integrators.hpp"
#include <stdio.h>  // fopen and fwrite

// data for the cooling function

static AthenaArray<Real> cool_t;
static AthenaArray<Real> cool_coef;
static AthenaArray<Real> cool_index;
static AthenaArray<Real> cool_tef;
static int nline=68;

static Real tunit;
static Real rhounit;
static Real time_unit;
static Real gamma_idx;
//======================================================================================
/*! \file beam.cpp
 *  \brief Beam test for the radiative transfer module
 *
 *====================================================================================*/


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================

static Real sigma=1.e8;

void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr, 
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc, Real dt);

void ExactCooling(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

void FixCRsource(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, Real time, Real dt, 
     int is, int ie, int js, int je, int ks, int ke);

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // isothermal sound speed for 1.e4 K is 1.16061e6 cm/s
  // length unit is pc=3.086*10^18 cm
  rhounit = pin->GetOrAddReal("problem", "rhounit", 1.25e-25);
  tunit = pin->GetOrAddReal("problem", "Tunit", 1.e4);
  time_unit = pin->GetOrAddReal("problem", "Timeunit", 2.65895e12);

  gamma_idx = pin->GetOrAddReal("hydro", "gamma", 5.0/3.0);

  FILE *fp;

  if(CR_ENABLED){
    cool_t.NewAthenaArray(nline);
    cool_coef.NewAthenaArray(nline);
    cool_index.NewAthenaArray(nline);
    cool_tef.NewAthenaArray(nline);

    if ((fp = fopen("cool_func.dat","r")) == NULL) {
        std::cout << "### ERROR to read in the cooling function data" << std::endl
              << "Cannot open cool_func.dat" << std::endl;
      return;
    }

    // temperature has nline numbers
    // coef and index only has nline-1 numbers

    for(int i=0; i<nline; ++i){
      fscanf(fp,"%lf",&(cool_t(i)));
      cool_t(i) /= tunit;
    }
    for(int i=0; i<nline; ++i){
      fscanf(fp,"%lf",&(cool_coef(i)));
    }
    for(int i=0; i<nline; ++i){
      fscanf(fp,"%lf",&(cool_index(i)));
    }
    
    // use X=0.7, Z=0.02
    Real miu =0.617284;
    Real miue=1.17647;
    Real miuh=1.42857;
    Real kb=1.3807e-16;
    Real mp=1.6605e-24;


    // Scale the unit
    for(int i=0; i<nline-1; ++i){
      cool_coef(i) *= (pow(tunit,cool_index(i)) * (gamma_idx-1.0) * miu 
                       * rhounit * time_unit/(tunit * kb * mp * miue * miuh));
    }
    // After scale the unit, the equation for cooling is just:
    // dT/dt = -coef T^alpha in code unit
    // The solution is (T_ref/Lambda_ref)(Y(T^n) - Y(T^n+1))=-Delta t

    // The TEF is just Lambda(T_ref)/T_{ref} \int _T ^ref dT
    // Starting from Npoint, d
    cool_tef(nline-1) = 0.0;
    Real lambda_t_n=cool_coef(nline-2)*pow(cool_t(nline-1),cool_index(nline-2))
                                                            /cool_t(nline-1);

    for(int i=nline-2; i>=0; i--){
      Real slope = cool_index(i);
      Real coef = cool_coef(i);
      if(fabs(slope-1.0) < TINY_NUMBER){
        cool_tef(i) = cool_tef(i+1) + lambda_t_n*log(cool_t(i+1)/cool_t(i))/coef;
      }else{
        cool_tef(i) = cool_tef(i+1) + lambda_t_n*(pow(cool_t(i+1),1.0-slope) - 
                                       pow(cool_t(i),1.0-slope))/(coef*(1.0-slope));        
      }

    }
 
    fclose(fp);
    

  }          
     
 
  EnrollUserExplicitSourceFunction(ExactCooling);
  EnrollUserBoundaryFunction(INNER_X1, FixCRsource);


}

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{ 

  if(CR_ENABLED){
     cool_t.DeleteAthenaArray();
     cool_coef.DeleteAthenaArray();
     cool_index.DeleteAthenaArray();
     cool_tef.DeleteAthenaArray();
  }
    
}


void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{

  if(CR_ENABLED){
    pcr->EnrollDiffFunction(Diffusion);
  }

}



void MeshBlock::ProblemGenerator(ParameterInput *pin)
{


  Real vx=0.0;
  Real rho_c = 1.0;
  Real rho_h = 0.001;
  Real delta_z = 25.0;
  Real z_back = 200.0;
  Real z_front = 200.0;
  Real pgas=1.0;
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {

        Real x1 = pcoord->x1v(i);
        Real density = rho_h + (rho_c - rho_h) * 
                               (1.0 + 1.0*tanh((x1-z_front)/delta_z))
                              *(1.0 + 1.0*tanh((z_back-x1)/delta_z));
      
        phydro->u(IDN,k,j,i) = density;
        phydro->u(IM1,k,j,i) = vx;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS){
          phydro->u(IEN,k,j,i) = 0.5*vx*vx+pgas/(gamma_idx-1.0);
        }
        
        if(CR_ENABLED){
            pcr->u_cr(CRE,k,j,i) = 1.e-6;
            pcr->u_cr(CRF1,k,j,i) = 0.0;
            pcr->u_cr(CRF2,k,j,i) = 0.0;
            pcr->u_cr(CRF3,k,j,i) = 0.0;
        }
      }// end i
    }
  }
  //Need to set opactiy sigma in the ghost zones
  if(CR_ENABLED){

  // Default values are 1/3
    int nz1 = block_size.nx1 + 2*(NGHOST);
    int nz2 = block_size.nx2;
    if(nz2 > 1) nz2 += 2*(NGHOST);
    int nz3 = block_size.nx3;
    if(nz3 > 1) nz3 += 2*(NGHOST);
    for(int k=0; k<nz3; ++k){
      for(int j=0; j<nz2; ++j){
        for(int i=0; i<nz1; ++i){
          pcr->sigma_diff(0,k,j,i) = sigma;
          pcr->sigma_diff(1,k,j,i) = sigma;
          pcr->sigma_diff(2,k,j,i) = sigma;
        }
      }
    }// end k,j,i

  }// End CR

    // Add horizontal magnetic field lines, to show streaming and diffusion 
  // along magnetic field ines
  if(MAGNETIC_FIELDS_ENABLED){

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          pfield->b.x1f(k,j,i) = 1.0;
        }
      }
    }

    if(block_size.nx2 > 1){

      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je+1; ++j) {
          for (int i=is; i<=ie; ++i) {
            pfield->b.x2f(k,j,i) = 0.0;
          }
        }
      }

    }

    if(block_size.nx3 > 1){

      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
    }// end nx3

    // set cell centerd magnetic field
    // Add magnetic energy density to the total energy
    pfield->CalculateCellCenteredField(pfield->b,pfield->bcc,pcoord,is,ie,js,je,ks,ke);

    for(int k=ks; k<=ke; ++k){
      for(int j=js; j<=je; ++j){
        for(int i=is; i<=ie; ++i){
          phydro->u(IEN,k,j,i) +=
            0.5*(SQR((pfield->bcc(IB1,k,j,i)))
               + SQR((pfield->bcc(IB2,k,j,i)))
               + SQR((pfield->bcc(IB3,k,j,i))));
      
        }
      }
    }

  }// end MHD
  
  
  return;
}



void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr, 
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc, Real dt)
{ 


  // set the default opacity to be a large value in the default hydro case
  CosmicRay *pcr=pmb->pcr;
  int kl=pmb->ks, ku=pmb->ke;
  int jl=pmb->js, ju=pmb->je;
  int il=pmb->is-1, iu=pmb->ie+1;
  if(pmb->block_size.nx2 > 1){
    jl -= 1;
    ju += 1;
  }
  if(pmb->block_size.nx3 > 1){
    kl -= 1;
    ku += 1;
  }

  for(int k=kl; k<=ku; ++k){
    for(int j=jl; j<=ju; ++j){
#pragma omp simd
      for(int i=il; i<=iu; ++i){

        pcr->sigma_diff(0,k,j,i) = sigma;
        pcr->sigma_diff(1,k,j,i) = sigma;
        pcr->sigma_diff(2,k,j,i) = sigma;  

      }
    }
  }

  Real invlim=1.0/pcr->vmax;

  // The information stored in the array
  // b_angle is
  // b_angle[0]=sin_theta_b
  // b_angle[1]=cos_theta_b
  // b_angle[2]=sin_phi_b
  // b_angle[3]=cos_phi_b




  if(MAGNETIC_FIELDS_ENABLED){
    //First, calculate B_dot_grad_Pc
    for(int k=kl; k<=ku; ++k){
      for(int j=jl; j<=ju; ++j){
    // x component
        pmb->pcoord->CenterWidth1(k,j,il-1,iu+1,pcr->cwidth);
        for(int i=il; i<=iu; ++i){
          Real distance = 0.5*(pcr->cwidth(i-1) + pcr->cwidth(i+1))
                         + pcr->cwidth(i);
          Real dprdx=(pcr->prtensor_cr(PC11,k,j,i+1) * u_cr(CRE,k,j,i+1)
                       - pcr->prtensor_cr(PC11,k,j,i-1) * u_cr(CRE,k,j,i-1));
          dprdx /= distance;
          pcr->b_grad_pc(k,j,i) = bcc(IB1,k,j,i) * dprdx;
        }
    //y component
        pmb->pcoord->CenterWidth2(k,j-1,il,iu,pcr->cwidth1);       
        pmb->pcoord->CenterWidth2(k,j,il,iu,pcr->cwidth);
        pmb->pcoord->CenterWidth2(k,j+1,il,iu,pcr->cwidth2);

        for(int i=il; i<=iu; ++i){
          Real distance = 0.5*(pcr->cwidth1(i) + pcr->cwidth2(i))
                         + pcr->cwidth(i);
          Real dprdy=(pcr->prtensor_cr(PC22,k,j+1,i) * u_cr(CRE,k,j+1,i)
                           - pcr->prtensor_cr(PC22,k,j-1,i) * u_cr(CRE,k,j-1,i));
          dprdy /= distance;
          pcr->b_grad_pc(k,j,i) += bcc(IB2,k,j,i) * dprdy;

        }
    // z component
        pmb->pcoord->CenterWidth3(k-1,j,il,iu,pcr->cwidth1);       
        pmb->pcoord->CenterWidth3(k,j,il,iu,pcr->cwidth);
        pmb->pcoord->CenterWidth3(k+1,j,il,iu,pcr->cwidth2);

        for(int i=il; i<=iu; ++i){
          Real distance = 0.5*(pcr->cwidth1(i) + pcr->cwidth2(i))
                          + pcr->cwidth(i);
          Real dprdz=(pcr->prtensor_cr(PC33,k+1,j,i) * u_cr(CRE,k+1,j,i)
                           - pcr->prtensor_cr(PC33,k-1,j,i) * u_cr(CRE,k-1,j,i));
          dprdz /= distance;
          pcr->b_grad_pc(k,j,i) += bcc(IB3,k,j,i) * dprdz;

          // now only get the sign
//          if(pcr->b_grad_pc(k,j,i) > TINY_NUMBER) pcr->b_grad_pc(k,j,i) = 1.0;
//          else if(-pcr->b_grad_pc(k,j,i) > TINY_NUMBER) pcr->b_grad_pc(k,j,i) 
//            = -1.0;
//          else pcr->b_grad_pc(k,j,i) = 0.0;
        }

      // now calculate the streaming velocity
      // streaming velocity is calculated with respect to the current coordinate 
      //  system
      // diffusion coefficient is calculated with respect to B direction
        for(int i=il; i<=iu; ++i){
          Real pb= bcc(IB1,k,j,i)*bcc(IB1,k,j,i)
                  +bcc(IB2,k,j,i)*bcc(IB2,k,j,i)
                  +bcc(IB3,k,j,i)*bcc(IB3,k,j,i);
          Real inv_sqrt_rho = 1.0/sqrt(prim(IDN,k,j,i));
          Real va1 = bcc(IB1,k,j,i)*inv_sqrt_rho;
          Real va2 = bcc(IB2,k,j,i)*inv_sqrt_rho;
          Real va3 = bcc(IB3,k,j,i)*inv_sqrt_rho;

          Real va = sqrt(pb/prim(IDN,k,j,i));

          Real dpc_sign = 0.0;
          if(pcr->b_grad_pc(k,j,i) > TINY_NUMBER) dpc_sign = 1.0;
          else if(-pcr->b_grad_pc(k,j,i) > TINY_NUMBER) dpc_sign = -1.0;
          
          pcr->v_adv(0,k,j,i) = -va1 * dpc_sign;
          pcr->v_adv(1,k,j,i) = -va2 * dpc_sign;
          pcr->v_adv(2,k,j,i) = -va3 * dpc_sign;

          // now the diffusion coefficient

          if(va < TINY_NUMBER){
            pcr->sigma_adv(0,k,j,i) = pcr->max_opacity;
          }else{
            pcr->sigma_adv(0,k,j,i) = fabs(pcr->b_grad_pc(k,j,i))
                                   /(va * (1.0 + pcr->prtensor_cr(PC11,k,j,i)) 
                                    * invlim * u_cr(CRE,k,j,i)); 
          }

          pcr->sigma_adv(1,k,j,i) = pcr->max_opacity;
          pcr->sigma_adv(2,k,j,i) = pcr->max_opacity;  

          // Now calculate the angles of B
          Real bxby = sqrt(bcc(IB1,k,j,i)*bcc(IB1,k,j,i) +
                           bcc(IB2,k,j,i)*bcc(IB2,k,j,i));
          Real btot = sqrt(pb);
          if(btot > TINY_NUMBER){
            pcr->b_angle(0,k,j,i) = bxby/btot;
            pcr->b_angle(1,k,j,i) = bcc(IB3,k,j,i)/btot;
          }else{
            pcr->b_angle(0,k,j,i) = 1.0;
            pcr->b_angle(1,k,j,i) = 0.0;
          }
          if(bxby > TINY_NUMBER){
            pcr->b_angle(2,k,j,i) = bcc(IB2,k,j,i)/bxby;
            pcr->b_angle(3,k,j,i) = bcc(IB1,k,j,i)/bxby;
          }else{
            pcr->b_angle(2,k,j,i) = 0.0;
            pcr->b_angle(3,k,j,i) = 1.0;            
          }

        }//        

      }// end j
    }// end k

  }// End MHD  
  else{



  for(int k=kl; k<=ku; ++k){
    for(int j=jl; j<=ju; ++j){
  // x component
      pmb->pcoord->CenterWidth1(k,j,il-1,iu+1,pcr->cwidth);
      for(int i=il; i<=iu; ++i){
         Real distance = 0.5*(pcr->cwidth(i-1) + pcr->cwidth(i+1))
                        + pcr->cwidth(i);
         Real grad_pr=(pcr->prtensor_cr(PC11,k,j,i+1) * u_cr(CRE,k,j,i+1)
                     - pcr->prtensor_cr(PC11,k,j,i-1) * u_cr(CRE,k,j,i-1));
         grad_pr /= distance;

         Real va = 1.0;

         if(va < TINY_NUMBER){
           pcr->sigma_adv(0,k,j,i) = sigma;
           pcr->v_adv(0,k,j,i) = 0.0;
         }else{
           Real sigma2 = fabs(grad_pr)/(va * (1.0 + pcr->prtensor_cr(PC11,k,j,i)) 
                             * invlim * u_cr(CRE,k,j,i)); 
           if(fabs(grad_pr) < TINY_NUMBER){
             pcr->sigma_adv(0,k,j,i) = 0.0;
             pcr->v_adv(0,k,j,i) = 0.0;
           }else{
             pcr->sigma_adv(0,k,j,i) = sigma2;
             pcr->v_adv(0,k,j,i) = -va * grad_pr/fabs(grad_pr);     
           }
        }

        pcr->sigma_adv(1,k,j,i) = pcr->max_opacity;
        pcr->sigma_adv(2,k,j,i) = pcr->max_opacity;
       
        pcr->v_adv(1,k,j,i) = 0.0;
        pcr->v_adv(2,k,j,i) = 0.0;




      }

    }
  }

  }
}

void ExactCooling(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{

    // scale the temperature unit and cooling time
    // so the cooling will be just 
    // dT/dt = coef() Lambda(T)


  int kl=pmb->ks, ku=pmb->ke;
  int jl=pmb->js, ju=pmb->je;
  int il=pmb->is, iu=pmb->ie;


  Real lambda_t_n=cool_coef(nline-2)*pow(cool_t(nline-1),cool_index(nline-2))
                                                            /cool_t(nline-1);



  for(int k=kl; k<=ku; ++k){
    for(int j=jl; j<=ju; ++j){
#pragma omp simd
      for(int i=il; i<=iu; ++i){
        // Need to take density and temperature at time step n from 
        // cons, not from prim
        // Because we do not need intermediate step to calculate the cooling 
        // function
        Real rho = cons(IDN,k,j,i);
        Real eint = cons(IEN,k,j,i) 
                    - 0.5 *(cons(IM1,k,j,i)*cons(IM1,k,j,i)
                      + cons(IM2,k,j,i)*cons(IM2,k,j,i)
                      + cons(IM3,k,j,i)*cons(IM3,k,j,i))/rho;
        if(MAGNETIC_FIELDS_ENABLED){
             eint -= 0.5 * (bcc(IB1,k,j,i) * bcc(IB1,k,j,i)
                      + bcc(IB2,k,j,i) * bcc(IB2,k,j,i)
                      + bcc(IB3,k,j,i) * bcc(IB3,k,j,i));
        }
        Real t_i = eint *(gamma_idx - 1.0)/rho;

        if(t_i > cool_t(0)){

          if(t_i > cool_t(nline-1)) t_i = cool_t(nline-1);

          int t_loc=0;
          while((t_loc < nline-2) && (cool_t(t_loc+1) < t_i) ){
            ++t_loc;
          }

          Real slope = cool_index(t_loc);
          Real coef = cool_coef(t_loc);

          Real tef = cool_tef(t_loc+1);
          if(fabs(slope-1.0) < TINY_NUMBER){
            tef += lambda_t_n*log(cool_t(t_loc+1)/t_i)/coef;
          }else{
            tef += lambda_t_n*(pow(cool_t(t_loc+1),1.0-slope) - 
                            pow(t_i,1.0-slope))/(coef*(1.0-slope));        
          }

          Real new_tef = tef + rho * dt * lambda_t_n;
          // Now invert TEF to get the current temperature
          // new_tef > tef
          int tef_loc=t_loc+1;
          while((tef_loc > 0) && (new_tef > cool_tef(tef_loc))){
            --tef_loc;
          }

          Real diff_tef = (new_tef - cool_tef(tef_loc+1))/lambda_t_n;
          slope = cool_index(tef_loc);
          coef = cool_coef(tef_loc);

          Real tnew = t_i;
          if(fabs(slope-1.0) < TINY_NUMBER){
            tnew = exp(log(cool_t(tef_loc+1))-(coef*diff_tef));
          }else{
            tnew = pow(cool_t(tef_loc+1),1.0-slope) 
                               - (1.0-slope) * coef * diff_tef;
            tnew = pow(tnew,1.0/(1.0-slope));
          }


          cons(IEN,k,j,i) += (tnew - t_i) * rho/(gamma_idx - 1.0);  

        }
      }
    }
  }




}



// reflecting boundary condition for hydro variables
// Fix U_cr, reflecting CR flux
void FixCRsource(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, Real time, Real dt, 
     int is, int ie, int js, int je, int ks, int ke)
{
  Real fix_u = 3.0;

  // copy hydro variables into ghost zones, reflecting v1
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVX)) {
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          prim(IVX,k,j,is-i) = -prim(IVX,k,j,(is+i-1));  // reflect 1-velocity
        }
      }}
    } else {
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          prim(n,k,j,is-i) = prim(n,k,j,(is+i-1));
        }
      }}
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b1
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) { 
    for (int j=js; j<=je; ++j) { 
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) { 
        b.x1f(k,j,(is-i)) = b.x1f(k,j,(is+i  ));  // reflect 1-field
      } 
    }}
  
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x2f(k,j,(is-i)) =  b.x2f(k,j,(is+i-1));
      }
    }}  
        
    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x3f(k,j,(is-i)) =  b.x3f(k,j,(is+i-1));
      }
    }}
  }


  if(CR_ENABLED){
    for (int n=0; n<(NCR); ++n) {
      if (n==(CRF1)) {
        for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(NGHOST); ++i) {
            u_cr(CRF1,k,j,is-i) = -u_cr(CRF1,k,j,(is+i-1));  // reflect 1-velocity
          }
        }}
      }else if (n==(CRE)){
        for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(NGHOST); ++i) {
            u_cr(CRE,k,j,is-i) = fix_u;  // reflect 1-velocity
          }
        }}        
      } 

      else {
        for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(NGHOST); ++i) {
            u_cr(n,k,j,is-i) = u_cr(n,k,j,(is+i-1));
          }
        }}
      }
    }
  }


}




