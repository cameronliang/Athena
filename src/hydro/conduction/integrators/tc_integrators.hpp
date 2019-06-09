#ifndef NTCINTEGRATORS_HPP
#define NTCINTEGRATORS_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file radiation.hpp
//  \brief definitions for Radiation class
//======================================================================================

// Athena++ classes headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../task_list/task_list.hpp"

class Hydro;
class ParameterInput;
class NewThermalConduction;

//! \class RadIntegrator
//  \brief integrate algorithm for radiative transfer


class NTCIntegrator {
  friend class NewThermalConduction;
public:
  NTCIntegrator(NewThermalConduction *ptc, ParameterInput *pin);
  ~NTCIntegrator();
  
  NewThermalConduction *pmy_tc;

    

  void AddSourceTerms(MeshBlock *pmb, const Real dt, AthenaArray<Real> &u,
                                  AthenaArray<Real> &u_tc, const int step);

  void CalculateFluxes(MeshBlock *pmb,
      AthenaArray<Real> &w, AthenaArray<Real> &bcc, AthenaArray<Real> &u_tc, 
      int reconstruct_order);


  void FluxDivergence(MeshBlock *pmb, AthenaArray<Real> &w,
                      AthenaArray<Real> &u_tc1, AthenaArray<Real> &u_tc2,
                      const IntegratorWeight wght, AthenaArray<Real> &u_out);

  void NTCFlux(int fdir, int il, int iu, 
      AthenaArray<Real> &t_l, AthenaArray<Real> &t_r,
      AthenaArray<Real> &rho_l, AthenaArray<Real> &rho_r,
      AthenaArray<Real> &w_l, AthenaArray<Real> &w_r,  
      AthenaArray<Real> &vdiff_l, AthenaArray<Real> &vdiff_r, 
                                      AthenaArray<Real> &flx);


  void DonorCellX1(const int k, const int j,
      const int il, const int iu, const AthenaArray<Real> &u_tc,
      AthenaArray<Real> &rho, AthenaArray<Real> &tgas,
      AthenaArray<Real> &rho_l, AthenaArray<Real> &rho_r,
      AthenaArray<Real> &t_l, AthenaArray<Real> &t_r,
      AthenaArray<Real> &w_l, AthenaArray<Real> &w_r);   
                                  
  void DonorCellX2(const int k, const int j,
      const int il, const int iu, const AthenaArray<Real> &u_tc,
      AthenaArray<Real> &rho, AthenaArray<Real> &tgas,
      AthenaArray<Real> &rho_l, AthenaArray<Real> &rho_r,
      AthenaArray<Real> &t_l, AthenaArray<Real> &t_r,
      AthenaArray<Real> &w_l, AthenaArray<Real> &w_r);


  void DonorCellX3(const int k, const int j,
      const int il, const int iu, const AthenaArray<Real> &u_tc,
      AthenaArray<Real> &rho, AthenaArray<Real> &tgas,
      AthenaArray<Real> &rho_l, AthenaArray<Real> &rho_r,
      AthenaArray<Real> &t_l, AthenaArray<Real> &t_r,
      AthenaArray<Real> &w_l, AthenaArray<Real> &w_r);


  void PieceWiseLinear(const int k, const int j,
      const int il, const int iu, AthenaArray<Real> &u_tc,
      AthenaArray<Real> &rho, AthenaArray<Real> &tgas,
      AthenaArray<Real> &rho_l, AthenaArray<Real> &rho_r,
      AthenaArray<Real> &t_l, AthenaArray<Real> &t_r,
      AthenaArray<Real> &w_l, AthenaArray<Real> &w_r, int fdir);

        
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
  AthenaArray<Real> flx_;
  AthenaArray<Real> wl_,wr_,vdiff_,vdiff_l_,vdiff_r_;
  AthenaArray<Real> rho_l_, rho_r_, tgas_l_, tgas_r_;
   //The final change of internal energy due to thermal conduction
  AthenaArray<Real> tc_esource_;

    // temporary array to store the flux
  Real taufact_;
  AthenaArray<Real> x1face_area_, x2face_area_, x3face_area_;
  AthenaArray<Real> x2face_area_p1_, x3face_area_p1_;
  AthenaArray<Real> cell_volume_, cwidth_;

};

#endif // NTCINTEGRATORS_HPP
