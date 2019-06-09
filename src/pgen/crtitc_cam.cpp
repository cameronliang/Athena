//======================================================================================
// the code runs now. last check it agrees down to temperature of tnew in python version. 
// i can continue to cool uniform gas to see what happens. 
// double check if it agrees everything in python version. 
// check if it agrees with exactcooling2..(general behavior)

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


static int ninputline=100001;

//temperature unit at 2e5K, cool_0=0.1414
//kappa_0=0.0414
//
static Real unit_n,unit_T,unit_rho, unit_length;
static Real unit_pressure,unit_velocity,unit_time;
//static Real heat_coef;
//static Real rhounit;
//static Real time_unit;
static Real gamma_idx;
static Real totP=1.5;

static AthenaArray<Real> cool_t;
static AthenaArray<Real> cool_coef;
static AthenaArray<Real> cool_index;
static AthenaArray<Real> cool_tef;


static Real const_pb=0.5;

static AthenaArray<Real> inputx;
static AthenaArray<Real> inputrho;
static AthenaArray<Real> inputtg;
static AthenaArray<Real> inputec;
static AthenaArray<Real> tbot;
static AthenaArray<Real> rhobot;
static AthenaArray<Real> ecbot;
static AthenaArray<Real> ttop;
static AthenaArray<Real> rhotop;
static AthenaArray<Real> ectop;

static Real cool_coef1=0.1414;
static Real heat_coef=0.0;
static Real cool_alpha1=-1.0;
static Real BOLTZ = 1.380648e-16; // erg/K
static Real mH    = 1.6726e-24; // g
static int nline=11;
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

void TCKappa(MeshBlock *pmb, 
	         AthenaArray<Real> &prim, AthenaArray<Real> &bcc);

void ExactCooling(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

void ExactCooling3(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

void FixCRsourceLeft(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, Real time, Real dt, 
     int is, int ie, int js, int je, int ks, int ke);
void FixCRsourceRight(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, Real time, Real dt, 
     int is, int ie, int js, int je, int ks, int ke);

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // isothermal sound speed for 1.e4 K is 1.16061e6 cm/s
  // length unit is pc=3.086*10^18 cm
  unit_n      = pin->GetOrAddReal("problem", "n0", 1.0e-1);
  unit_T      = pin->GetOrAddReal("problem", "T0", 1.0e5);
  unit_length = pin->GetOrAddReal("problem", "domainsize", 3.086e18); // 1 pc

  unit_rho      = unit_n * mH;
  unit_pressure = BOLTZ * unit_n * unit_T;
  unit_velocity = sqrt(unit_pressure / unit_rho);
  unit_time     = unit_length / unit_velocity;

  //heat_coef = pin->GetOrAddReal("problem", "heat_coef", 0.1);

 // time_unit = pin->GetOrAddReal("problem", "Timeunit", 3.3497e11);

  gamma_idx = pin->GetOrAddReal("hydro", "gamma", 5.0/3.0);



  //const_heatrate = pin->GetOrAddReal("problem", "heat_rate", 0.3);
 
  EnrollUserExplicitSourceFunction(ExactCooling3);
//  EnrollUserBoundaryFunction(INNER_X1, FixCRsourceLeft);
//  EnrollUserBoundaryFunction(OUTER_X1, FixCRsourceRight);


  // the picesewise power law cooling

  cool_t.NewAthenaArray(nline);

  cool_coef.NewAthenaArray(nline);
  cool_index.NewAthenaArray(nline);
  cool_tef.NewAthenaArray(nline);

  cool_t(0) = 10.000000 / unit_T;
  cool_t(1) = 53.366992 / unit_T;  
  cool_t(2) = 284.803587 / unit_T;
  cool_t(3) = 1519.911083 / unit_T;
  cool_t(4) = 8111.308308 / unit_T;
  cool_t(5) = 43287.612811 / unit_T;
  cool_t(6) = 231012.970008 / unit_T;
  cool_t(7) = 1232846.739442 / unit_T;
  cool_t(8) = 6579332.246576 / unit_T;
  cool_t(9) = 35111917.342151 / unit_T;
  cool_t(10) = 187381742.286039 / unit_T;
  

  cool_coef(0) = 1.000000e-15;
  cool_coef(1) = 1.000000e-15;
  cool_coef(2) = 1.000000e-15;
  cool_coef(3) = 1.000000e-15;
  cool_coef(4) = 3.995462e-01;
  cool_coef(5) = 1.319280e+01;
  cool_coef(6) = 2.491086e+01;
  cool_coef(7) = 5.856870e+00;
  cool_coef(8) = 1.629738e+00;
  cool_coef(9) = 1.000864e+00;
  cool_coef(10) = 2.014851e+00;


  cool_index(0) = 2.390287e-15;
  cool_index(1) = 4.243041e-15;
  cool_index(2) = 1.246475e-15;
  cool_index(3) = 9.077155e+00;
  cool_index(4) = 2.088309e+00;
  cool_index(5) = 3.795713e-01;
  cool_index(6) = -8.644943e-01;
  cool_index(7) = -7.638781e-01;
  cool_index(8) = -2.911465e-01;
  cool_index(9) = 4.178183e-01;
  cool_index(10) = 8.006204e-02;
  
  //Real lambda_t_n=cool_coef(nline-2)*pow(cool_t(nline-1),cool_index(nline-2))
  //                                                          /cool_t(nline-1);
  Real lambda_t_n=cool_coef(nline-2); // last Lambda. units?
  cool_tef(nline-1) = 0.0; // Last Y_N = 0

  // Calculate All Yk. (Eq. A6)
  // loop backward using Y_N = 0asdas
  for(int i=nline-2; i>=0; i--){
    Real slope = cool_index(i);
    Real coef = cool_coef(i);
    if(fabs(slope-1.0) < TINY_NUMBER){
      cool_tef(i) = cool_tef(i+1) - (lambda_t_n/coef)*(cool_t(i)/cool_t(nline-1))*log(cool_t(i+1)/cool_t(i));

    }else{
      Real term = (lambda_t_n/coef)*(cool_t(i)/cool_t(nline-1)) / (1-slope);

      Real a = (lambda_t_n/coef); 
      Real b = (cool_t(i)/cool_t(nline-1));
      Real c = 1-slope; 
      printf("(cool_t(%d) = %e\n", i, cool_t(i));
      //printf("a = %e,b = %e,c = %e\n", a,b,c);

      printf("1 = %d, %e\n", i, term);
      term *= 1-pow(cool_t(i)/cool_t(i+1),slope-1);
      printf("2 = %d, %e\n", i, term);
      cool_tef(i) = cool_tef(i+1) - term;
      
      //std::cout << "cool_tef,i = " << cool_tef(i) << ", " <<i << std::endl;
      
    }
    //printf("%d, %e\n", i, cool_tef(i)); 
  }
  
}


void MeshBlock::UserWorkInLoop(void)
{
/*

    int il=is, iu=ie, jl=js, ju=je, kl=ks, ku=ke;
    il -= NGHOST;
    iu += NGHOST;
    if(ju>jl){
       jl -= NGHOST;
       ju += NGHOST;
    }
    if(ku>kl){
      kl -= NGHOST;
      ku += NGHOST;
    }
    Real gm1 = peos->GetGamma() - 1.0;

    Real wi[(NWAVE)];
    Real cfmax;
    
    
    for (int k=kl; k<=ku; ++k){
      for (int j=jl; j<=ju; ++j){
       for (int i=il; i<=iu; ++i){

          Real& vx=phydro->w(IVX,k,j,i);
          Real& vy=phydro->w(IVY,k,j,i);
          Real& vz=phydro->w(IVZ,k,j,i);

          Real& rho=phydro->w(IDN,k,j,i);

          Real vel = sqrt(vx*vx+vy*vy+vz*vz);
          Real ke = 0.5 * rho * vel * vel;

          Real pb=0.0;

          int flag=0;

          // case 1, check superlum velocity
          if(vel > pcr->vlim * pcr->vmax){
            Real ratio = pcr->vmax * pcr->vlim / vel;
            vx *= ratio;
            vy *= ratio;
            vz *= ratio;

            phydro->u(IDN,k,j,i) = rho;
            phydro->u(IM1,k,j,i) = rho*vx;
            phydro->u(IM2,k,j,i) = rho*vy;
            phydro->u(IM3,k,j,i) = rho*vz;

            ke = 0.5 * rho * (vx*vx+vy*vy+vz*vz);
           
            flag=1;
          }
         
          // case 2, check fast speed too large
          if(MAGNETIC_FIELDS_ENABLED){
         
            pb = 0.5*(SQR(pfield->bcc(IB1,k,j,i))+SQR(pfield->bcc(IB2,k,j,i))
                 +SQR(pfield->bcc(IB3,k,j,i)));
          }

          Real tgas = phydro->w(IEN,k,j,i)/rho;
          if(tgas < 0.1) {tgas = 0.1; flag=1;}
          if(tgas > 10.0) {tgas = 10.0; flag=1;}
          if(flag > 0){
             // Only do this for bad cells
             user_out_var(2,k,j,i) = phydro->u(IEN,k,j,i);
             Real  eint = rho * tgas/gm1;
             phydro->u(IEN,k,j,i) = eint + ke + pb;
             phydro->w(IEN,k,j,i) = rho * tgas;
             user_out_var(2,k,j,i) = phydro->u(IEN,k,j,i)
                                   - user_out_var(2,k,j,i);

          }

      }}}
    
*/
  return;

}



void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{ 

  cool_t.DeleteAthenaArray();
  cool_coef.DeleteAthenaArray();
  cool_index.DeleteAthenaArray();
  cool_tef.DeleteAthenaArray();  

}


void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{


  AllocateUserOutputVariables(4);

  if(CR_ENABLED)
    pcr->EnrollDiffFunction(Diffusion);

  if(NEW_TH_CON_ENABLED){

      phydro->ptc->EnrollKappaFunction(TCKappa);

  }

}



void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

  std::srand(gid);
  Real amp = 0.1;

  // Initialize hydro variable
  for(int i=is; i<=ie; ++i) {
    Real &x1 = pcoord->x1v(i);
    
    
    Real rho = 1.0;

    Real tem = 1.0;

    Real pc = 0.5;

    Real va = sqrt(2.0*const_pb/rho);
    Real vm = 100.0;
    if(CR_ENABLED) vm = pcr->vmax;
    Real fc = va * (pc*3.0 + pc)/vm;


    for(int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {

        phydro->u(IDN,k,j,i) = rho * (1.0 + amp 
        	                       * ((double)rand()/(double)RAND_MAX-0.5));
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;

        if (NON_BAROTROPIC_EOS){

          phydro->u(IEN,k,j,i) = tem * rho/(gamma_idx-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }
        
        if(CR_ENABLED){
            pcr->u_cr(CRE,k,j,i) = pc*3.0 * (1.0 + amp 
        	                       * ((double)rand()/(double)RAND_MAX-0.5));
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
          pfield->b.x1f(k,j,i) = sqrt(const_pb);
        }
      }
    }

    if(block_size.nx2 > 1){

      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je+1; ++j) {
          for (int i=is; i<=ie; ++i) {
            pfield->b.x2f(k,j,i) = sqrt(const_pb);
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



void TCKappa(MeshBlock *pmb, AthenaArray<Real> &prim, AthenaArray<Real> &bcc)
{ 
  // set the default opacity to be a large value in the default hydro case
  NewThermalConduction *ptc=pmb->phydro->ptc;
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
        Real tgas = prim(IEN,k,j,i)/prim(IDN,k,j,i);
        Real kappa_parallel = std::min(0.005*tgas*tgas*sqrt(tgas),0.005);
        Real kappa_pernd = 0.00001;
        kappa_parallel = std::max(kappa_parallel, 10.0*kappa_pernd);

        ptc->ntc_kappa(0,k,j,i) = kappa_parallel;
        ptc->ntc_kappa(1,k,j,i) = kappa_pernd;
        ptc->ntc_kappa(2,k,j,i) = kappa_pernd;

      }
    }
  }

  if(MAGNETIC_FIELDS_ENABLED){
    //First, calculate B_dot_grad_Pc
    for(int k=kl; k<=ku; ++k){
      for(int j=jl; j<=ju; ++j){
      // diffusion coefficient is calculated with respect to B direction
      // Use a simple estimate of Grad Pc
        for(int i=il; i<=iu; ++i){
          // Now calculate the angles of B
          Real bxby = sqrt(bcc(IB1,k,j,i)*bcc(IB1,k,j,i) +
                           bcc(IB2,k,j,i)*bcc(IB2,k,j,i));
          Real btot = sqrt(bcc(IB1,k,j,i)*bcc(IB1,k,j,i) +
                           bcc(IB2,k,j,i)*bcc(IB2,k,j,i) +
                           bcc(IB3,k,j,i)*bcc(IB3,k,j,i));

          if(btot > TINY_NUMBER){
            ptc->b_angle(0,k,j,i) = bxby/btot;
            ptc->b_angle(1,k,j,i) = bcc(IB3,k,j,i)/btot;
          }else{
            ptc->b_angle(0,k,j,i) = 1.0;
            ptc->b_angle(1,k,j,i) = 0.0;
          }
          if(bxby > TINY_NUMBER){
            ptc->b_angle(2,k,j,i) = bcc(IB2,k,j,i)/bxby;
            ptc->b_angle(3,k,j,i) = bcc(IB1,k,j,i)/bxby;
          }else{
            ptc->b_angle(2,k,j,i) = 0.0;
            ptc->b_angle(3,k,j,i) = 1.0;
          }


        }//end i        

      }// end j
    }// end k
  }// end MHD



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
    // dT/dt =  -cool_coef *rho * T^alpha
  // for the current unit, cool_coef=1.4142


  int kl=pmb->ks, ku=pmb->ke;
  int jl=pmb->js, ju=pmb->je;
  int il=pmb->is, iu=pmb->ie;



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
          Real newt= pow(t_i,1.0-cool_alpha1)-(1.0-cool_alpha1)*rho*cool_coef1*dt;
          newt = pow(newt,1.0/(1.0-cool_alpha1));

//        Real newt2 = t_i - cool_coef * rho * pow(t_i,cool_alpha) * dt;
//        newt = std::max(newt2,newt);

          cons(IEN,k,j,i) += (newt - t_i) * rho/(gamma_idx - 1.0); 

          pmb->user_out_var(0,k,j,i) = (newt - t_i) * rho/(gamma_idx - 1.0);
 
 //       pmb->user_out_var(0,k,j,i) = ((newt - t_i) * rho/(gamma_idx - 1.0))/dt;
 //       pmb->user_out_var(1,k,j,i) = -cool_coef1 * rho * pow(t_i, cool_alpha1) * rho/(gamma_idx - 1.0);  
        // add a constant heating rate
          cons(IEN,k,j,i) += dt * heat_coef/(gamma_idx-1.0);     
          pmb->user_out_var(1,k,j,i) = dt * heat_coef/(gamma_idx-1.0);
      }     
    }
  }

}


void ExactCooling2(MeshBlock *pmb, const Real time, const Real dt,
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

          int t_loc=0; // location of current temperature bin. 
          while((t_loc < nline-2) && (cool_t(t_loc+1) < t_i) ){
            ++t_loc;
          }

          Real slope = cool_index(t_loc);
          Real coef = cool_coef(t_loc);

          Real tef = cool_tef(t_loc+1); // 
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

          pmb->user_out_var(0,k,j,i) = (tnew - t_i) * rho/(gamma_idx - 1.0);

          //add a constant heating rate
          if(tnew < 100.0){
            cons(IEN,k,j,i) += dt * heat_coef/(gamma_idx-1.0);
            //printf("Here..\n");

            pmb->user_out_var(1,k,j,i) = dt * heat_coef/(gamma_idx-1.0);
          }else{
            pmb->user_out_var(1,k,j,i) = 0.0;          	
          }

        }
      }
    }
  }

}

// ==================================================================
// ==================================================================
// ==================================================================
// ==================================================================
void ExactCooling3(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  int kl=pmb->ks, ku=pmb->ke;
  int jl=pmb->js, ju=pmb->je;
  int il=pmb->is, iu=pmb->ie;

  Real radE_total = 0; 
  int V = (ju-jl+1) * (iu-il+1);

  Real lambda_t_n=cool_coef(nline-2); // last Lambda?


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
        // Get current temperature 
        Real t_i = eint *(gamma_idx - 1.0)/rho;


        if(t_i > cool_t(0)){

          if(t_i > cool_t(nline-1)) t_i = cool_t(nline-1);

          int t_loc=0; // location of current temperature bin. 
          while((t_loc < nline-2) && (cool_t(t_loc+1) < t_i) ){
            ++t_loc;
          }

          Real slope = cool_index(t_loc);
          Real coef = cool_coef(t_loc);

          // Compute current Y for current temperature (Eq. A5)
          // since each temperature is not exactly at bin k. --> add a little bit Y to it. 
          Real tef = cool_tef(t_loc+1); 



          Real term1 = (lambda_t_n/coef)*(cool_t(t_loc)/cool_t(nline-1)) / (1-slope);
          term1 *= (1-pow(cool_t(t_loc)/t_i,slope-1));
          tef = cool_tef(t_loc) + term1;
          //if(fabs(slope-1.0) < TINY_NUMBER){
          //  tef = cool_tef(t_loc) + (lambda_t_n/coef)*(t_i/cool_t(nline-1))*log(cool_t(t_loc)/cool_t(t_loc));
          //}else{
          //  Real term = (lambda_t_n/coef)*(t_i/cool_t(nline-1)) / (1-slope);
          //  term *= (1-pow(t_i/cool_t(t_loc),slope-1));
          //  tef = cool_tef(t_loc) + term;
            //tef += lambda_t_n*(pow(cool_t(t_loc+1),1.0-slope) - 
            //                pow(t_i,1.0-slope))/(coef*(1.0-slope));        
          //}


          Real constant_mu = 0.62; 
          if (t_i*unit_T < 1.0e4){
            constant_mu = 1.1;
          }

          // compute new Y (TEF) in the new timestep (Eq. 26)
          Real tcool_i_cgs = (gamma_idx-1.0)*BOLTZ*(t_i*unit_T) / ((rho*unit_rho/constant_mu/mH)*(coef*1.0e-23));
          Real tcool_i = tcool_i_cgs / unit_time;
          Real new_tef = tef + (t_i / cool_t(nline-1))*(lambda_t_n/coef) * (dt / tcool_i);
          // Now invert TEF to get the current temperature
          // new_tef > tef


          int tef_loc=t_loc+1;
          while((tef_loc > 0) && (new_tef > cool_tef(tef_loc))){
            --tef_loc;
          }


          Real diff_tef = (new_tef - cool_tef(tef_loc));
          slope = cool_index(tef_loc);
          coef = cool_coef(tef_loc);

          // compute new temperature by inv of Y (Eq. A7)
          Real tnew = t_i;
          if(fabs(slope-1.0) < TINY_NUMBER){
            tnew = (coef/lambda_t_n)*(cool_t(nline-1)/cool_t(tef_loc))*diff_tef;
            tnew = cool_t(tef_loc)*exp(-tnew);
          }else{
            tnew = (coef/lambda_t_n)*(cool_t(nline-1)/cool_t(tef_loc))*diff_tef;
            tnew = cool_t(tef_loc)*pow(1.0 - (1.0-slope)*tnew, 1.0/(1.0-slope));
            //tnew = pow(1 - tnew, 1.0/(1.0-slope));
            //tnew = cool_t(tef_loc) * tnew;
            //tnew = 1 - (1.0-slope) * (coef/lambda_t_n) *(cool_t(nline-1)/cool_t(tef_loc)) * diff_tef;
            //tnew = pow(tnew,1.0/(1.0-slope));
          }
          cons(IEN,k,j,i) += (tnew - t_i) * rho/(gamma_idx - 1.0);  
          radE_total += (tnew - t_i) * rho/(gamma_idx - 1.0);  
          pmb->user_out_var(0,k,j,i) = (tnew - t_i) * rho/(gamma_idx - 1.0);

          //add a constant heating rate
          if(tnew < 100.0){
            cons(IEN,k,j,i) += dt * heat_coef/(gamma_idx-1.0);
            //printf("Here..\n");

            pmb->user_out_var(1,k,j,i) = dt * heat_coef/(gamma_idx-1.0);
          }else{
            pmb->user_out_var(1,k,j,i) = 0.0;          	
          }

        }
      }
    }
  }

  for(int k=kl; k<=ku; ++k){
    for(int j=jl; j<=ju; ++j){
#pragma omp simd
      for(int i=il; i<=iu; ++i){
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
        // Get current temperature 
        //Real t_i = eint *(gamma_idx - 1.0)/rho;
        cons(IEN,k,j,i) -= radE_total /  (float) V; 
        //if(t_i*unit_T < 1.e7){
        //  cons(IEN,k,j,i) -= radE_total /  (float) V; 
        //}
        //cons(IEN,k,j,i) -= pmb->user_out_var(0,k,j,i) /  (float) V; 
        
      }
    }
  }


}
// ==================================================================
// ==================================================================
// ==================================================================
// ==================================================================



// reflecting boundary condition for hydro variables
// Fix U_cr, reflecting CR flux
void FixCRsourceLeft(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, Real time, Real dt, 
     int is, int ie, int js, int je, int ks, int ke)
{

  // copy hydro variables into ghost zones, reflecting v1

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        Real ecratio=ecbot(i-1)/u_cr(CRE,k,j,is);
        prim(IDN,k,j,is-i) = prim(IDN,k,j,is) * pow(ecratio,1.5);
//        prim(IDN,k,j,is-i) = rhobot(i-1);
        prim(IVX,k,j,is-i) = prim(IVX,k,j,is);  // reflect 1-velocity
        prim(IVY,k,j,is-i) = prim(IVY,k,j,is);
        prim(IVZ,k,j,is-i) = prim(IVZ,k,j,is);
        Real sum_press=prim(IEN,k,j,is) + u_cr(CRE,k,j,is)/3.0;
        prim(IEN,k,j,is-i) = std::max(sum_press - ecbot(i-1)/3.0,1.e-7);
      }
    }}

  

  // copy face-centered magnetic fields into ghost zones, reflecting b1
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) { 
    for (int j=js; j<=je; ++j) { 
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) { 
//        b.x1f(k,j,(is-i)) = sqrt(2.0*const_pb);  // reflect 1-field
          b.x1f(k,j,(is-i)) =  b.x1f(k,j,is);
      } 
    }}
    if(je > js){ 
     for (int k=ks; k<=ke; ++k) {
     for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x2f(k,j,(is-i)) =  b.x2f(k,j,is);
      }
     }}  
    }
    if(ke > ks){        
     for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
       for (int i=1; i<=(NGHOST); ++i) {
         b.x3f(k,j,(is-i)) =  b.x3f(k,j,is);
       }
      }}
    }
  }


  if(CR_ENABLED){
    for (int n=0; n<(NCR); ++n) {
      if (n==(CRF1)) {
        for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(NGHOST); ++i) {
            Real pb = (b.x1f(k,j,(is-i)) * b.x1f(k,j,(is-i)) 
                    + b.x2f(k,j,(is-i)) * b.x2f(k,j,(is-i))
                    + b.x3f(k,j,(is-i)) * b.x3f(k,j,(is-i)));
            Real va = sqrt(pb/prim(IDN,k,j,is-i));
            Real vm = pmb->pcr->vmax;
            Real fc = (va +  prim(IVX,k,j,is-i)) * 
                      (ecbot(i-1) + ecbot(i-1)/3.0)/vm;
            u_cr(CRF1,k,j,is-i) = fc;  // reflect 1-velocity

          }
        }}
      }else if (n==(CRE)){
        for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(NGHOST); ++i) {
            u_cr(CRE,k,j,is-i) = ecbot(i-1);  // reflect 1-velocity
          }
        }}        
      } 

      else {
        for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(NGHOST); ++i) {
            u_cr(n,k,j,is-i) = u_cr(n,k,j,is);
          }
        }}
      }
    }
  }
}


// reflecting boundary condition for hydro variables
// Fix U_cr, reflecting CR flux
void FixCRsourceRight(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, Real time, Real dt, 
     int is, int ie, int js, int je, int ks, int ke)
{

  // copy hydro variables into ghost zones, reflecting v1

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
//        Real ecratio=ectop(i-1)/ecbot(0);
        Real ecratio = ectop(i-1)/u_cr(CRE,k,j,ie);
//        prim(IDN,k,j,ie+i) = rhobot(0) * pow(ecratio,1.5);
 //       prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie);
//        prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie) * pow(ecratio,1.5);
//        prim(IDN,k,j,ie+i) = rhotop(i-1);
        prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie-i+1);        
        prim(IVX,k,j,ie+i) = -prim(IVX,k,j,ie-i+1);  // reflect 1-velocity
        prim(IVY,k,j,ie+i) = prim(IVY,k,j,ie-i+1);
        prim(IVZ,k,j,ie+i) = prim(IVZ,k,j,ie-i+1);
        Real sum_press = prim(IEN,k,j,ie) + u_cr(CRE,k,j,ie)/3.0;
        prim(IEN,k,j,ie+i) = std::max(sum_press-ectop(i-1)/3.0,1.e-7);
      }
    }}

  

  // copy face-centered magnetic fields into ghost zones, reflecting b1
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) { 
    for (int j=js; j<=je; ++j) { 
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) { 
        b.x1f(k,j,(ie+i+1)) = b.x1f(k,j,(ie+1));  // reflect 1-field
      } 
    }}
    if(je > js){  
     for (int k=ks; k<=ke; ++k) {
     for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x2f(k,j,(ie+i)) =  b.x2f(k,j,ie);
      }
     }}  
    }
    if(ke > ks){  
     for (int k=ks; k<=ke+1; ++k) {
     for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x3f(k,j,(ie+i)) =  b.x3f(k,j,ie);
      }
     }}
    }
  }


  if(CR_ENABLED){
/*    for (int n=0; n<(NCR); ++n) {
      if (n==(CRF1)) {
        for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(NGHOST); ++i) {

            Real va = sqrt(2.0*const_pb/prim(IDN,k,j,ie+i));
            Real vm = pmb->pcr->vmax;
            Real fc = (va + prim(IVX,k,j,ie+i)) * 
                      (ectop(i-1) + ectop(i-1)/3.0)/vm;

            u_cr(CRF1,k,j,ie+i) = fc;  // reflect 1-velocity
          }
        }}
      }else if (n==(CRE)){
        for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(NGHOST); ++i) {
            u_cr(CRE,k,j,ie+i) = ectop(i-1);  // reflect 1-velocity

          }
        }}        
      } else {
        for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(NGHOST); ++i) {
            u_cr(n,k,j,ie+i) = u_cr(n,k,j,(ie-i+1));
          }
        }}
      }
    }
*/
        for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(NGHOST); ++i) {
            u_cr(CRE,k,j,ie+i) = ectop(i-1);
            Real pb = (b.x1f(k,j,(ie+i)) * b.x1f(k,j,(ie+i)) 
                    + b.x2f(k,j,(ie+i)) * b.x2f(k,j,(ie+i))
                    + b.x3f(k,j,(ie+i)) * b.x3f(k,j,(ie+i)));

            Real va = sqrt(pb/prim(IDN,k,j,ie+i));
            Real vm = pmb->pcr->vmax;
            Real fc = (va + prim(IVX,k,j,ie+i)) * 
                      (u_cr(CRE,k,j,ie+i) * 4.0/3.0)/vm;

            u_cr(CRF1,k,j,ie+i) = fc;             
            u_cr(CRF2,k,j,ie+i) = u_cr(CRF2,k,j,ie-i+1);
            u_cr(CRF3,k,j,ie+i) = u_cr(CRF3,k,j,ie-i+1);

          }
        }}


  }


}



