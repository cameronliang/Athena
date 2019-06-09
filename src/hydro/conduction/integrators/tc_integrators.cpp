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
//! \file rad_integrators.cpp
//  \brief implementation of radiation integrators
//======================================================================================


// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../parameter_input.hpp"
#include "../../../mesh/mesh.hpp"
#include "../tc.hpp"
#include "tc_integrators.hpp"
#include "../../hydro.hpp"


NTCIntegrator::NTCIntegrator(NewThermalConduction *ptc, ParameterInput *pin)
{

  pmy_tc = ptc;
  

  taufact_ = pin->GetOrAddReal("tc","taucell",5.0);
  int nthreads = ptc->pmy_hydro->pmy_block->pmy_mesh->GetNumMeshThreads();
  int ncells1 = ptc->pmy_hydro->pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1;
  int ncells3 = 1;

  flx_.NewAthenaArray(nthreads,4,ncells1);
  vdiff_l_.NewAthenaArray(nthreads,ncells1);
  vdiff_r_.NewAthenaArray(nthreads,ncells1);
  rho_l_.NewAthenaArray(nthreads,ncells1);
  rho_r_.NewAthenaArray(nthreads,ncells1);
  tgas_l_.NewAthenaArray(nthreads,ncells1);
  tgas_r_.NewAthenaArray(nthreads,ncells1);
  wl_.NewAthenaArray(nthreads,4,ncells1);
  wr_.NewAthenaArray(nthreads,4,ncells1);

    
  cwidth_.NewAthenaArray(nthreads,ncells1);


  
  x1face_area_.NewAthenaArray(nthreads,ncells1+1);
  if(ptc->pmy_hydro->pmy_block->block_size.nx2 > 1) {
    x2face_area_.NewAthenaArray(nthreads,ncells1);
    x2face_area_p1_.NewAthenaArray(nthreads,ncells1);
    ncells2 = ptc->pmy_hydro->pmy_block->block_size.nx2 + 2*(NGHOST);

  }
  if(ptc->pmy_hydro->pmy_block->block_size.nx3 > 1) {
    x3face_area_.NewAthenaArray(nthreads,ncells1);
    x3face_area_p1_.NewAthenaArray(nthreads,ncells1);
    ncells3 = ptc->pmy_hydro->pmy_block->block_size.nx3 + 2*(NGHOST);

  }
  cell_volume_.NewAthenaArray(nthreads,ncells1);


  tc_esource_.NewAthenaArray(ncells3,ncells2,ncells1);
  vdiff_.NewAthenaArray(3,ncells3,ncells2,ncells1);

}

// destructor

NTCIntegrator::~NTCIntegrator()
{
  flx_.DeleteAthenaArray();
  vdiff_l_.DeleteAthenaArray();
  vdiff_r_.DeleteAthenaArray();

  wl_.DeleteAthenaArray();
  wr_.DeleteAthenaArray();
  rho_l_.DeleteAthenaArray();
  rho_r_.DeleteAthenaArray();
  tgas_l_.DeleteAthenaArray();
  tgas_r_.DeleteAthenaArray();
  tc_esource_.DeleteAthenaArray();
  vdiff_.DeleteAthenaArray();


  cwidth_.DeleteAthenaArray();

  
  x1face_area_.DeleteAthenaArray();
  if(pmy_tc->pmy_hydro->pmy_block->block_size.nx2 > 1) {
    x2face_area_.DeleteAthenaArray();
    x2face_area_p1_.DeleteAthenaArray();
  }
  if(pmy_tc->pmy_hydro->pmy_block->block_size.nx3 > 1) {
    x3face_area_.DeleteAthenaArray();
    x3face_area_p1_.DeleteAthenaArray();
  }
  cell_volume_.DeleteAthenaArray();


}




void NTCIntegrator::RotateVec(const Real sint, const Real cost, 
                 const Real sinp, const Real cosp, 
                    Real &v1, Real &v2, Real &v3)
{
  // vel1, vel2, vel3 are input
  // v1, v2, v3 are output
  // The two rotation matrix 
  //R_1=
  //[cos_p  sin_p 0]
  //[-sin_p cos_p 0]
  //[0       0    1]

  //R_2=
  //[sin_t  0 cos_t]
  //[0      1    0]
  //[-cos_t 0 sin_t]


  // First apply R1, then apply R2
  Real newv1 =  cosp * v1 + sinp * v2;
  v2 = -sinp * v1 + cosp * v2;

  // now apply R2
  v1 =  sint * newv1 + cost * v3;
  Real newv3 = -cost * newv1 + sint * v3;
  v3 = newv3;

}

void NTCIntegrator::InvRotateVec(const Real sint, const Real cost, 
                 const Real sinp, const Real cosp, 
                 Real &v1, Real &v2, Real &v3)
{
  // vel1, vel2, vel3 are input
  // v1, v2, v3 are output
  // The two rotation matrix 
  //R_1^-1=
  //[cos_p  -sin_p 0]
  //[sin_p cos_p 0]
  //[0       0    1]

  //R_2^-1=
  //[sin_t  0 -cos_t]
  //[0      1    0]
  //[cos_t 0 sin_t]


  // First apply R2^-1, then apply R1^-1
  Real newv1 = sint * v1 - cost * v3;
  v3 = cost * v1 + sint * v3;

  // now apply R1^-1
  v1 = cosp * newv1 - sinp * v2;
  Real newv2 = sinp * newv1 + cosp * v2;
  v2 = newv2;

}

