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
#include "../radiation/radiation.hpp"

// File scope variables
static int ang;
static int octnum;

//======================================================================================
/*! \file beam.cpp
 *  \brief Beam test for the radiative transfer module
 *
 *====================================================================================*/

void TwoBeams(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  ang = pin->GetOrAddInteger("problem","ang",0);
  octnum = pin->GetOrAddInteger("problem","octnum",0);
  
    // Enroll boundary functions
  if(RADIATION_ENABLED)
  EnrollUserRadBoundaryFunction(INNER_X2, TwoBeams);

  return;
}



//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  
  Real gamma = peos->GetGamma();
  
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phydro->u(IDN,k,j,i) = 1.0;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS){

          phydro->u(IEN,k,j,i) = 1.0/(gamma-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }
      }
    }
  }
  
  //Now initialize opacity and specific intensity
  if(RADIATION_ENABLED){
    int nfreq = prad->nfreq;
    int nang = prad->nang;
    for(int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ifr=0; ifr < nfreq; ++ifr){
            prad->sigma_s(k,j,i,ifr) = 0.0;
            prad->sigma_a(k,j,i,ifr) = 0.0;
            prad->sigma_ae(k,j,i,ifr) = 0.0;
          }
          for(int n=0; n<prad->n_fre_ang; ++n){
              prad->ir(k,j,i,n) = 0.0;
          }
        }
      }
    }
  }
  
  return;
}




void TwoBeams(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
         Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  Radiation *prad=pmb->prad;
  int nang=prad->nang;
  int noct=prad->noct;
  int nfreq=prad->nfreq;
  int ang_oct=nang/noct;
  
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      for (int i=is; i<=ie; ++i) {
        Real &x1 = pco->x1v(i);
        Real &x2 = pco->x2v(j);
        for(int ifr=0; ifr<nfreq; ++ifr){
        for(int l=0; l<noct; ++l){
        for(int n=0; n<ang_oct; ++n){
          int n_ang=l*ang_oct + n;
          Real slope1=prad->mu(1,k,j,i,0)/prad->mu(0,k,j,i,0);
          Real slope2=prad->mu(1,k,j,i,n_ang)/prad->mu(0,k,j,i,n_ang);
          Real dis1=fabs(slope1*(x1-0.1)+(x2+2.0));
          Real dis2=fabs(slope2*(x1+0.1)+(x2+2.0));
          if(((l==0)&&(n==0)&&(dis1<pco->dx1v(i))) ||
               ((l==1)&&(n==0)&&(dis2<pco->dx1v(i)))){
            a(k,js-j,i,n_ang+ifr*nang) = 10.0;
          }else{
            a(k,js-j,i,n_ang+ifr*nang) = 0.0;
          }
        }
        }
        }

    }}
  }

  return;
}



