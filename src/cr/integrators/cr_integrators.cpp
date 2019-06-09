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
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"
#include "../cr.hpp"
#include "cr_integrators.hpp"



CRIntegrator::CRIntegrator(CosmicRay *pcr, ParameterInput *pin)
{

  pmy_cr = pcr;
  

  taufact_ = pin->GetOrAddReal("cr","taucell",1.0);
  vel_flx_flag_ = pin->GetOrAddInteger("cr","vflx",1);
  int nthreads = pcr->pmy_block->pmy_mesh->GetNumMeshThreads();
  int ncells1 = pcr->pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1;
  int ncells3 = 1;

  flx_.NewAthenaArray(nthreads,NCR,ncells1);
  vel_l_.NewAthenaArray(nthreads,ncells1);
  vel_r_.NewAthenaArray(nthreads,ncells1);
  vdiff_l_.NewAthenaArray(nthreads,ncells1);
  vdiff_r_.NewAthenaArray(nthreads,ncells1);
  eddl_.NewAthenaArray(nthreads,6,ncells1);
  eddr_.NewAthenaArray(nthreads,6,ncells1);
  wl_.NewAthenaArray(nthreads,(NCR),ncells1);
  wr_.NewAthenaArray(nthreads,(NCR),ncells1);
  new_sol_.NewAthenaArray(NCR,ncells1);



    
  cwidth_.NewAthenaArray(nthreads,ncells1);
  
  x1face_area_.NewAthenaArray(nthreads,ncells1+1);
  if(pcr->pmy_block->block_size.nx2 > 1) {
    x2face_area_.NewAthenaArray(nthreads,ncells1);
    x2face_area_p1_.NewAthenaArray(nthreads,ncells1);
    ncells2 = pcr->pmy_block->block_size.nx2 + 2*(NGHOST);
  }
  if(pcr->pmy_block->block_size.nx3 > 1) {
    x3face_area_.NewAthenaArray(nthreads,ncells1);
    x3face_area_p1_.NewAthenaArray(nthreads,ncells1);
    ncells3 = pcr->pmy_block->block_size.nx3 + 2*(NGHOST);
  }
  cell_volume_.NewAthenaArray(nthreads,ncells1);


  grad_pc_.NewAthenaArray(3,ncells3,ncells2,ncells1);
  ec_source_.NewAthenaArray(ncells3,ncells2,ncells1);

}

// destructor

CRIntegrator::~CRIntegrator()
{
  flx_.DeleteAthenaArray();
  vel_l_.DeleteAthenaArray();
  vel_r_.DeleteAthenaArray();
  vdiff_l_.DeleteAthenaArray();
  vdiff_r_.DeleteAthenaArray();
  eddl_.DeleteAthenaArray();
  eddr_.DeleteAthenaArray();
  wl_.DeleteAthenaArray();
  wr_.DeleteAthenaArray();
  new_sol_.DeleteAthenaArray();



  cwidth_.DeleteAthenaArray();
  
  x1face_area_.DeleteAthenaArray();
  if(pmy_cr->pmy_block->block_size.nx2 > 1) {
    x2face_area_.DeleteAthenaArray();
    x2face_area_p1_.DeleteAthenaArray();
  }
  if(pmy_cr->pmy_block->block_size.nx3 > 1) {
    x3face_area_.DeleteAthenaArray();
    x3face_area_p1_.DeleteAthenaArray();
  }
  cell_volume_.DeleteAthenaArray();

  grad_pc_.DeleteAthenaArray();
  ec_source_.DeleteAthenaArray();


}




void CRIntegrator::RotateVec(const Real sint, const Real cost, 
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

void CRIntegrator::InvRotateVec(const Real sint, const Real cost, 
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

