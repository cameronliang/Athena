//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file history.cpp
//  \brief writes history output data, volume-averaged quantities that are output
//         frequently in time to trace their history.

// C/C++ headers
#include <sstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../globals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../mesh/mesh.hpp"
#include "../radiation/radiation.hpp"
#include "../cr/cr.hpp"
#include "outputs.hpp"
//14 radiation variables, 
// if no RADIATION_ENABLED, they are always 0
#if (RADIATION_ENABLED > 0)
  #define NHISTORY_VARS ((NHYDRO)+(NFIELD)+3+14+(NCR))
  #define NRAD (14)
#else
  #define NHISTORY_VARS ((NHYDRO)+(NFIELD)+3+(NCR))
  #define NRAD 0
#endif
//----------------------------------------------------------------------------------------
// HistoryOutput constructor
// destructor - not needed for this derived class

HistoryOutput::HistoryOutput(OutputParameters oparams)
  : OutputType(oparams)
{
}

//----------------------------------------------------------------------------------------
//! \fn void OutputType::HistoryFile()
//  \brief Writes a history file

void HistoryOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag)
{
  MeshBlock *pmb=pm->pblock;
  AthenaArray<Real> vol;

  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  vol.NewAthenaArray(ncells1);
  int nhistory_output=NHISTORY_VARS+pm->nuser_history_output_;

  Real *data_sum = new Real[nhistory_output];
  for (int n=0; n<nhistory_output; ++n) data_sum[n]=0.0;

  // Loop over MeshBlocks
  while (pmb != NULL) {
    Hydro *phyd = pmb->phydro;
    Field *pfld = pmb->pfield;
    Radiation *prad = pmb->prad;
    CosmicRay *pcr=pmb->pcr;

    // Sum history variables over cells.  Note ghost cells are never included in sums
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      pmb->pcoord->CellVolume(k,j,pmb->is,pmb->ie,vol);
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real& u_d  = phyd->u(IDN,k,j,i);
        Real& u_mx = phyd->u(IM1,k,j,i);
        Real& u_my = phyd->u(IM2,k,j,i);
        Real& u_mz = phyd->u(IM3,k,j,i);

        data_sum[0] += vol(i)*u_d;
        data_sum[1] += vol(i)*u_mx;
        data_sum[2] += vol(i)*u_my;
        data_sum[3] += vol(i)*u_mz;
        data_sum[4] += vol(i)*0.5*SQR(u_mx)/u_d;
        data_sum[5] += vol(i)*0.5*SQR(u_my)/u_d;
        data_sum[6] += vol(i)*0.5*SQR(u_mz)/u_d;

        if (NON_BAROTROPIC_EOS) {
          Real& u_e = phyd->u(IEN,k,j,i);;
          data_sum[7] += vol(i)*u_e;
        }
        if (MAGNETIC_FIELDS_ENABLED) {
          Real& bcc1 = pfld->bcc(IB1,k,j,i);
          Real& bcc2 = pfld->bcc(IB2,k,j,i);
          Real& bcc3 = pfld->bcc(IB3,k,j,i);
          data_sum[NHYDRO + 3] += vol(i)*0.5*bcc1*bcc1;
          data_sum[NHYDRO + 4] += vol(i)*0.5*bcc2*bcc2;
          data_sum[NHYDRO + 5] += vol(i)*0.5*bcc3*bcc3;
        }// End MHD
        if (RADIATION_ENABLED){
          data_sum[NHYDRO + NFIELD + 3] += vol(i)*prad->rad_mom(IER,k,j,i);
          data_sum[NHYDRO + NFIELD + 4] += vol(i)*prad->rad_mom(IFR1,k,j,i);
          data_sum[NHYDRO + NFIELD + 5] += vol(i)*prad->rad_mom(IFR2,k,j,i);
          data_sum[NHYDRO + NFIELD + 6] += vol(i)*prad->rad_mom(IFR3,k,j,i);
          data_sum[NHYDRO + NFIELD + 7] += vol(i)*prad->rad_mom_cm(IER,k,j,i);
          data_sum[NHYDRO + NFIELD + 8] += vol(i)*prad->rad_mom_cm(IFR1,k,j,i);
          data_sum[NHYDRO + NFIELD + 9] += vol(i)*prad->rad_mom_cm(IFR2,k,j,i);
          data_sum[NHYDRO + NFIELD + 10] += vol(i)*prad->rad_mom_cm(IFR3,k,j,i);
          data_sum[NHYDRO + NFIELD + 11] += vol(i)*prad->rad_mom(IPR11,k,j,i);
          data_sum[NHYDRO + NFIELD + 12] += vol(i)*prad->rad_mom(IPR12,k,j,i);
          data_sum[NHYDRO + NFIELD + 13] += vol(i)*prad->rad_mom(IPR13,k,j,i);
          data_sum[NHYDRO + NFIELD + 14] += vol(i)*prad->rad_mom(IPR22,k,j,i);
          data_sum[NHYDRO + NFIELD + 15] += vol(i)*prad->rad_mom(IPR23,k,j,i);
          data_sum[NHYDRO + NFIELD + 16] += vol(i)*prad->rad_mom(IPR33,k,j,i);

        }
        if (CR_ENABLED){
          data_sum[NHYDRO + NFIELD + 3 + NRAD] += vol(i)*pcr->u_cr(CRE,k,j,i);
          data_sum[NHYDRO + NFIELD + 4 + NRAD] += vol(i)*pcr->u_cr(CRF1,k,j,i);
          data_sum[NHYDRO + NFIELD + 5 + NRAD] += vol(i)*pcr->u_cr(CRF2,k,j,i);
          data_sum[NHYDRO + NFIELD + 6 + NRAD] += vol(i)*pcr->u_cr(CRF3,k,j,i);
        }
      }
    }}
    for(int n=0; n<pm->nuser_history_output_; n++) { // user-defined history outputs
      if(pm->user_history_func_[n]!=NULL)
        data_sum[NHISTORY_VARS+n] += pm->user_history_func_[n](pmb, n);
    }
    pmb=pmb->next;
  }  // end loop over MeshBlocks

#ifdef MPI_PARALLEL
  // sum over all ranks
  if (Globals::my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, data_sum, nhistory_output, MPI_ATHENA_REAL, MPI_SUM, 0,
               MPI_COMM_WORLD);
  } else {
    MPI_Reduce(data_sum, data_sum, nhistory_output, MPI_ATHENA_REAL, MPI_SUM, 0,
               MPI_COMM_WORLD);
  }
#endif

  // only the master rank writes the file
  // create filename: "file_basename" + ".hst".  There is no file number.
  if (Globals::my_rank == 0) {
    std::string fname;
    fname.assign(output_params.file_basename);
    fname.append(".hst");

    // open file for output
    FILE *pfile;
    std::stringstream msg;
    if((pfile = fopen(fname.c_str(),"a")) == NULL){
      msg << "### FATAL ERROR in function [OutputType::HistoryFile]" << std::endl
          << "Output file '" << fname << "' could not be opened";
      throw std::runtime_error(msg.str().c_str());
    }

    // If this is the first output, write header
    int iout = 1;
    if (output_params.file_number == 0) {
      fprintf(pfile,"# Athena++ history data\n"); // descriptor is first line
      fprintf(pfile,"# [%d]=time     ", iout++);
      fprintf(pfile,"[%d]=dt       ", iout++);
      fprintf(pfile,"[%d]=mass     ", iout++);
      fprintf(pfile,"[%d]=1-mom    ", iout++);
      fprintf(pfile,"[%d]=2-mom    ", iout++);
      fprintf(pfile,"[%d]=3-mom    ", iout++);
      fprintf(pfile,"[%d]=1-KE     ", iout++);
      fprintf(pfile,"[%d]=2-KE     ", iout++);
      fprintf(pfile,"[%d]=3-KE     ", iout++);
      if (NON_BAROTROPIC_EOS) fprintf(pfile,"[%d]=tot-E   ", iout++);
      if (MAGNETIC_FIELDS_ENABLED) {
        fprintf(pfile,"[%d]=1-ME    ", iout++);
        fprintf(pfile,"[%d]=2-ME    ", iout++);
        fprintf(pfile,"[%d]=3-ME    ", iout++);
      }
      if (RADIATION_ENABLED){
        fprintf(pfile,"[%d]=Er    ", iout++);
        fprintf(pfile,"[%d]=1-Fr    ", iout++);
        fprintf(pfile,"[%d]=2-Fr    ", iout++);
        fprintf(pfile,"[%d]=3-Fr    ", iout++);
        fprintf(pfile,"[%d]=Er0    ", iout++);
        fprintf(pfile,"[%d]=1-Fr0    ", iout++);
        fprintf(pfile,"[%d]=2-Fr0    ", iout++);
        fprintf(pfile,"[%d]=3-Fr0    ", iout++);
        fprintf(pfile,"[%d]=Pr11    ", iout++);
        fprintf(pfile,"[%d]=Pr12    ", iout++);
        fprintf(pfile,"[%d]=Pr13    ", iout++);
        fprintf(pfile,"[%d]=Pr22    ", iout++);
        fprintf(pfile,"[%d]=Pr23    ", iout++);
        fprintf(pfile,"[%d]=Pr33    ", iout++);
      }
      if (CR_ENABLED){
        fprintf(pfile,"[%d]=Ec    ", iout++);
        fprintf(pfile,"[%d]=1-Fc    ", iout++);
        fprintf(pfile,"[%d]=2-Fc    ", iout++);
        fprintf(pfile,"[%d]=3-Fc    ", iout++);
      }
      for(int n=0; n<pm->nuser_history_output_; n++)
        fprintf(pfile,"[%d]=%-8s", iout++, pm->user_history_output_names_[n].c_str());
      fprintf(pfile,"\n");                              // terminate line
    }

    // write history variables
    fprintf(pfile, output_params.data_format.c_str(), pm->time);
    fprintf(pfile, output_params.data_format.c_str(), pm->dt);
    for (int n=0; n<nhistory_output; ++n)
      fprintf(pfile, output_params.data_format.c_str(), data_sum[n]);
    fprintf(pfile,"\n"); // terminate line
    fclose(pfile);
  }

  // increment counters, clean up
  output_params.file_number++;
  output_params.next_time += output_params.dt;
  pin->SetInteger(output_params.block_name, "file_number", output_params.file_number);
  pin->SetReal(output_params.block_name, "next_time", output_params.next_time);
  vol.DeleteAthenaArray();
  delete [] data_sum;
  return;
}
