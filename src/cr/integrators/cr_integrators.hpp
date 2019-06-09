#ifndef CRINTEGRATORS_HPP
#define CRINTEGRATORS_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file radiation.hpp
//  \brief definitions for Radiation class
//======================================================================================

// Athena++ classes headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../cr.hpp" // radiation
#include "../../task_list/task_list.hpp"

class MeshBlock;
class ParameterInput;
class CosmicRay;

//! \class RadIntegrator
//  \brief integrate algorithm for radiative transfer


class CRIntegrator {
  friend class CosmicRay;
public:
  CRIntegrator(CosmicRay *pcr, ParameterInput *pin);
  ~CRIntegrator();
  
  CosmicRay *pmy_cr;
    

  void AddSourceTerms(MeshBlock *pmb, const Real dt, AthenaArray<Real> &u,
        AthenaArray<Real> &w, AthenaArray<Real> &bcc, 
        AthenaArray<Real> &u_cr, const int step);

  void CalculateFluxes(MeshBlock *pmb,
      AthenaArray<Real> &w, AthenaArray<Real> &bcc, AthenaArray<Real> &u_cr, 
      int reconstruct_order);




  void FluxDivergence(MeshBlock *pmb, AthenaArray<Real> &u_cr1,
                      AthenaArray<Real> &u_cr2, const IntegratorWeight wght,
                      AthenaArray<Real> &u_out, AthenaArray<Real> &u, 
                      AthenaArray<Real> &w, AthenaArray<Real> &bcc);

  void CRFlux(int dir, int k, int j, int il, int iu, AthenaArray<Real> &w_l, 
      AthenaArray<Real> &w_r, AthenaArray<Real> &vel_l, AthenaArray<Real> &vel_r,  
      AthenaArray<Real> &eddl, AthenaArray<Real> &eddr,      
      AthenaArray<Real> &vdiff_l, AthenaArray<Real> &vdiff_r, AthenaArray<Real> &flx);


  void DonorCellX1(const int k, const int j,
      const int il, const int iu, const AthenaArray<Real> &u_cr,
      const AthenaArray<Real> &prim, const AthenaArray<Real> &edd, 
      AthenaArray<Real> &w_l, AthenaArray<Real> &w_r, 
      AthenaArray<Real> &v_l, AthenaArray<Real> &v_r,
      AthenaArray<Real> &eddl, AthenaArray<Real> &eddr);   
                                  
  void DonorCellX2(const int k, const int j,
      const int il, const int iu, const AthenaArray<Real> &u_cr,
      const AthenaArray<Real> &prim, const AthenaArray<Real> &edd, 
      AthenaArray<Real> &w_l, AthenaArray<Real> &w_r, 
      AthenaArray<Real> &v_l, AthenaArray<Real> &v_r,
      AthenaArray<Real> &eddl, AthenaArray<Real> &eddr);   

  void DonorCellX3(const int k, const int j,
      const int il, const int iu, const AthenaArray<Real> &u_cr,
      const AthenaArray<Real> &prim, const AthenaArray<Real> &edd,
      AthenaArray<Real> &w_l, AthenaArray<Real> &w_r, 
      AthenaArray<Real> &v_l, AthenaArray<Real> &v_r,
      AthenaArray<Real> &eddl, AthenaArray<Real> &eddr);  


  void PieceWiseLinear(const int k, const int j,
      const int il, const int iu, AthenaArray<Real> &u_cr,
      AthenaArray<Real> &prim, AthenaArray<Real> &edd,
      AthenaArray<Real> &w_l, AthenaArray<Real> &w_r, 
      AthenaArray<Real> &v_l, AthenaArray<Real> &v_r, 
      AthenaArray<Real> &eddl, AthenaArray<Real> &eddr, int dir);

  void GetOneVariableX1(const int k, const int j, 
      const int il, const int iu, const AthenaArray<Real> &q, 
      AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void GetOneVariableX2(const int k, const int j, 
      const int il, const int iu, const AthenaArray<Real> &q, 
      AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void GetOneVariableX3(const int k, const int j, 
      const int il, const int iu, const AthenaArray<Real> &q, 
      AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void RotateVec(const Real sint, const Real cost, 
                 const Real sinp, const Real cosp, 
                     Real &v1, Real &v2, Real &v3);

  void InvRotateVec(const Real sint, const Real cost, 
                 const Real sinp, const Real cosp, 
                     Real &v1, Real &v2, Real &v3);


private:
  AthenaArray<Real> new_sol_;
  AthenaArray<Real> flx_;
  AthenaArray<Real> vel_l_,vel_r_,wl_,wr_,vdiff_l_,vdiff_r_;
  AthenaArray<Real> eddl_,eddr_;
  AthenaArray<Real> grad_pc_, ec_source_;

    // temporary array to store the flux
  Real taufact_;
  int vel_flx_flag_;
  AthenaArray<Real> x1face_area_, x2face_area_, x3face_area_;
  AthenaArray<Real> x2face_area_p1_, x3face_area_p1_;
  AthenaArray<Real> cell_volume_, cwidth_;

};

#endif // CRINTEGRATORS_HPP
