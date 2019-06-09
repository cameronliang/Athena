//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file time_integrator.cpp
//  \brief derived class for time integrator task list.  Can create task lists for one
//  of many different time integrators (e.g. van Leer, RK2, RK3, etc.)

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "task_list.hpp"
#include "../athena.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../bvals/bvals.hpp"
#include "../eos/eos.hpp"
#include "../radiation/radiation.hpp"
#include "../radiation/integrators/rad_integrators.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../cr/cr.hpp"
#include "../cr/integrators/cr_integrators.hpp"
#include "../hydro/conduction/tc.hpp"
#include "../hydro/conduction/integrators/tc_integrators.hpp"

//----------------------------------------------------------------------------------------
//  TimeIntegratorTaskList constructor

TimeIntegratorTaskList::TimeIntegratorTaskList(ParameterInput *pin, Mesh *pm)
  : TaskList(pm)
{
  // First, set weights for each step of time-integration algorithm.  Each step is
  //    U^{2} = a*U^0 + b*U^1 + c*dt*Div(F), where U^0 and U^1 are previous steps 
  // a,b=(1-a),and c are weights that are different for each step and each integrator
  // These are stored as: time_int_wght1 = a, time_int_wght2 = b, time_int_wght3 = c

  integrator = pin->GetOrAddString("time","integrator","vl2");

  uint64_t step_flag;

  // second-order van Leer integrator (Gardiner & Stone, NewA 14, 139 2009)
  if (integrator == "vl2") {
    nsub_steps = 2;
    step_wghts[0].a = 1.0;
    step_wghts[0].b = 0.0;
    step_wghts[0].c = 0.5;

    step_wghts[1].a = 1.0;
    step_wghts[1].b = 0.0;
    step_wghts[1].c = 1.0;
  } else if (integrator == "rk2") {
    nsub_steps = 2;
    step_wghts[0].a = 1.0;
    step_wghts[0].b = 0.0;
    step_wghts[0].c = 1.0;

    step_wghts[1].a = 0.5;
    step_wghts[1].b = 0.5;
    step_wghts[1].c = 0.5;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in CreateTimeIntegrator" << std::endl
        << "integrator=" << integrator << " not valid time integrator" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // Now assemble list of tasks for each step of time integrator
  {using namespace HydroIntegratorTaskNames;
    AddTimeIntegratorTask(START_ALLRECV,NONE);

    // compute hydro fluxes, integrate hydro variables
    AddTimeIntegratorTask(CALC_HYDFLX,START_ALLRECV);
    if(RADIATION_ENABLED)
      AddTimeIntegratorTask(CALC_RADFLX,START_ALLRECV);
    if(CR_ENABLED)
      AddTimeIntegratorTask(CALC_CRFLX,START_ALLRECV);  
    if(NEW_TH_CON_ENABLED){
      AddTimeIntegratorTask(INI_TC,NONE);
      AddTimeIntegratorTask(CALC_TCFLX,START_ALLRECV|INI_TC);
    }   
    if(pm->multilevel==true) { // SMR or AMR
      step_flag = CALC_HYDFLX;
      if(RADIATION_ENABLED)
        step_flag = (step_flag|CALC_RADFLX);
      if(CR_ENABLED)
        step_flag = (step_flag|CALC_CRFLX);
      if(NEW_TH_CON_ENABLED)
        step_flag = (step_flag|CALC_TCFLX);

      AddTimeIntegratorTask(SEND_HYDFLX,step_flag);
      AddTimeIntegratorTask(RECV_HYDFLX,step_flag);

      AddTimeIntegratorTask(INT_HYD, RECV_HYDFLX);
      if(RADIATION_ENABLED)
        AddTimeIntegratorTask(INT_RAD, RECV_HYDFLX);
      if(CR_ENABLED)
        AddTimeIntegratorTask(INT_CR, RECV_HYDFLX);
      if(NEW_TH_CON_ENABLED)
        AddTimeIntegratorTask(INT_TC, RECV_HYDFLX);        
    } else {
      AddTimeIntegratorTask(INT_HYD, CALC_HYDFLX);
      if(RADIATION_ENABLED)
        AddTimeIntegratorTask(INT_RAD, CALC_RADFLX); 
      if(CR_ENABLED)
        AddTimeIntegratorTask(INT_CR, CALC_CRFLX); 
      if(NEW_TH_CON_ENABLED)
        AddTimeIntegratorTask(INT_TC, CALC_TCFLX);
    }
    AddTimeIntegratorTask(SRCTERM_HYD,INT_HYD);
    if(CR_ENABLED)
      AddTimeIntegratorTask(SRCTERM_CR,INT_CR);
    if(RADIATION_ENABLED)
      AddTimeIntegratorTask(SRCTERM_RAD,INT_RAD);
    if(NEW_TH_CON_ENABLED)
      AddTimeIntegratorTask(SRCTERM_TC,INT_TC);    

    step_flag =  SRCTERM_HYD;
    if(CR_ENABLED)
      step_flag = (step_flag|SRCTERM_CR);
    if(RADIATION_ENABLED)
      step_flag = (step_flag|SRCTERM_RAD);
    if(NEW_TH_CON_ENABLED)
      step_flag = (step_flag|SRCTERM_TC);

    AddTimeIntegratorTask(SEND_HYD,step_flag);
   

    AddTimeIntegratorTask(RECV_HYD,START_ALLRECV);

    // compute MHD fluxes, integrate field
    if (MAGNETIC_FIELDS_ENABLED) { // MHD
      AddTimeIntegratorTask(CALC_FLDFLX,CALC_HYDFLX);
      AddTimeIntegratorTask(SEND_FLDFLX,CALC_FLDFLX);
      AddTimeIntegratorTask(RECV_FLDFLX,SEND_FLDFLX);
      AddTimeIntegratorTask(INT_FLD, RECV_FLDFLX);
      AddTimeIntegratorTask(SEND_FLD,INT_FLD);
      AddTimeIntegratorTask(RECV_FLD,START_ALLRECV);
    }

    // prolongate, compute new primitives
    if (MAGNETIC_FIELDS_ENABLED) { // MHD
      if(pm->multilevel==true) { // SMR or AMR
        // Radiation MR boundary and prolongation is done with hydro together.
        AddTimeIntegratorTask(PROLONG, (SEND_HYD|RECV_HYD|SEND_FLD|RECV_FLD));
        AddTimeIntegratorTask(CON2PRIM,PROLONG);
      } else {
        AddTimeIntegratorTask(CON2PRIM,(INT_HYD|RECV_HYD|INT_FLD|RECV_FLD));
      }
    } else {  // HYDRO
      if(pm->multilevel==true) { // SMR or AMR
        AddTimeIntegratorTask(PROLONG,(SEND_HYD|RECV_HYD));
        AddTimeIntegratorTask(CON2PRIM,PROLONG);
      } else {
        AddTimeIntegratorTask(CON2PRIM,(INT_HYD|RECV_HYD));
      }
    }

    // Apply Physical Boundary Condition
    AddTimeIntegratorTask(PHY_BVAL,CON2PRIM);

    // Radiation physical boundary
    if(RADIATION_ENABLED){
      if(pm->multilevel==true) { // SMR or AMR
        AddTimeIntegratorTask(RADPHY_BVAL,PROLONG);
      }else{
        AddTimeIntegratorTask(RADPHY_BVAL,INT_RAD|RECV_HYD);
      }
      AddTimeIntegratorTask(RAD_MOMOPACITY,RADPHY_BVAL|PHY_BVAL);      
    }// End Radiation


    if(CR_ENABLED)
      AddTimeIntegratorTask(CR_VAOPACITY,PHY_BVAL); 

    if(NEW_TH_CON_ENABLED)
      AddTimeIntegratorTask(TC_OPACITY,PHY_BVAL);

    step_flag = PHY_BVAL;
    if(RADIATION_ENABLED)
      step_flag = (step_flag|RAD_MOMOPACITY);
    if(CR_ENABLED)
      step_flag = (step_flag|CR_VAOPACITY);
    if(NEW_TH_CON_ENABLED)
      step_flag = (step_flag|TC_OPACITY);

    AddTimeIntegratorTask(USERWORK,step_flag);

    AddTimeIntegratorTask(NEW_DT,USERWORK);
    if(pm->adaptive==true) {
      AddTimeIntegratorTask(AMR_FLAG,USERWORK);
      AddTimeIntegratorTask(CLEAR_ALLRECV,AMR_FLAG);
    } else {
      AddTimeIntegratorTask(CLEAR_ALLRECV,NEW_DT);
    }

  } // end of using namespace block
}

//----------------------------------------------------------------------------------------//! \fn
//  \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//  value of ntask.  

void TimeIntegratorTaskList::AddTimeIntegratorTask(uint64_t id, uint64_t dep)
{
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace HydroIntegratorTaskNames;
  switch((id)) {
    case (START_ALLRECV):
      task_list_[ntasks].TaskFunc= 
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::StartAllReceive);
      break;
    case (CLEAR_ALLRECV):
      task_list_[ntasks].TaskFunc= 
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ClearAllReceive);
      break;

    case (CALC_HYDFLX):
      task_list_[ntasks].TaskFunc= 
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CalculateFluxes);
      break;
    case (CALC_FLDFLX):
      task_list_[ntasks].TaskFunc= 
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CalculateEMF);
      break;

    case (SEND_HYDFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FluxCorrectSend);
      break;
    case (SEND_FLDFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::EMFCorrectSend);
      break;

    case (RECV_HYDFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FluxCorrectReceive);
      break;
    case (RECV_FLDFLX):
      task_list_[ntasks].TaskFunc= 
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::EMFCorrectReceive);
      break;

    case (INT_HYD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::HydroIntegrate);
      break;
    case (INT_FLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FieldIntegrate);
      break;

    case (SRCTERM_HYD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::HydroSourceTerms);
      break;

    case (SEND_HYD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::HydroSend);
      break;
    case (SEND_FLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FieldSend);
      break;

    case (RECV_HYD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::HydroReceive);
      break;
    case (RECV_FLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FieldReceive);
      break;

    case (PROLONG):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::Prolongation);
      break;
    case (CON2PRIM):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::Primitives);
      break;
    case (PHY_BVAL):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::PhysicalBoundary);
      break;
    case (USERWORK):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::UserWork);
      break;
    case (NEW_DT):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::NewBlockTimeStep);
      break;
    case (AMR_FLAG):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CheckRefinement);
      break;
    //case for radiation
    case (CALC_RADFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::RadFluxes);    
      break;  
    case (SRCTERM_RAD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::RadSourceTerms);
      break;
    case (INT_RAD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::RadIntegrate);
      break;
    case (RAD_MOMOPACITY):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::RadMomOpacity);
      break;
    case (RADPHY_BVAL):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::RadPhysicalBoundary);  
      break;
    // case for Cosmic Ray
    case (CALC_CRFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CRFluxes);
      break;
    case (INT_CR):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CRIntegrate);
      break;
    case (SRCTERM_CR):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CRSourceTerms);
      break;
    case (CR_VAOPACITY):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CRVAOpacity);
      break;
    //case for thermal conduction
    case (INI_TC):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::TCInitialize);
      break;
    case (CALC_TCFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::TCFluxes);
      break;
    case (INT_TC):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::TCIntegrate);
      break;
    case (SRCTERM_TC):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::TCSourceTerms);
      break;
    case (TC_OPACITY):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::TCOpacity);
      break;


    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in AddTask" << std::endl
          << "Invalid Task "<< id << " is specified" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  ntasks++;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief

//----------------------------------------------------------------------------------------
// Functions to start/end MPI communication

enum TaskStatus TimeIntegratorTaskList::StartAllReceive(MeshBlock *pmb, int step)
{
  pmb->pbval->StartReceivingAll();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::ClearAllReceive(MeshBlock *pmb, int step)
{
  pmb->pbval->ClearBoundaryAll();
  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions to calculates fluxes

enum TaskStatus TimeIntegratorTaskList::CalculateFluxes(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;

  if((step == 1) && (integrator == "vl2")) {
    phydro->CalculateFluxes(phydro->w,  pfield->b,  pfield->bcc, 1);
    return TASK_NEXT;
  }

  if((step == 1) && (integrator == "rk2")) {
    phydro->CalculateFluxes(phydro->w,  pfield->b,  pfield->bcc, 2);
    return TASK_NEXT;
  } 

  if(step == 2) {
    phydro->CalculateFluxes(phydro->w1, pfield->b1, pfield->bcc1, 2);
    return TASK_NEXT;
  }

  return TASK_FAIL;
}

enum TaskStatus TimeIntegratorTaskList::CalculateEMF(MeshBlock *pmb, int step)
{
  if(step == 1) {
    pmb->pfield->ComputeCornerE(pmb->phydro->w,  pmb->pfield->bcc);
    return TASK_NEXT;
  }

  if(step == 2) {
    pmb->pfield->ComputeCornerE(pmb->phydro->w1, pmb->pfield->bcc1);
    return TASK_NEXT;
  } 

  return TASK_FAIL;
}

//----------------------------------------------------------------------------------------
// Functions to communicate fluxes between MeshBlocks for flux correction step with AMR

enum TaskStatus TimeIntegratorTaskList::FluxCorrectSend(MeshBlock *pmb, int step)
{
  pmb->pbval->SendFluxCorrection(FLUX_HYDRO);
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::EMFCorrectSend(MeshBlock *pmb, int step)
{
  pmb->pbval->SendEMFCorrection();
  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions to receive fluxes between MeshBlocks

enum TaskStatus TimeIntegratorTaskList::FluxCorrectReceive(MeshBlock *pmb, int step)
{
  if(pmb->pbval->ReceiveFluxCorrection(FLUX_HYDRO) == true) {
    return TASK_NEXT;
  } else {
    return TASK_FAIL;
  }
}

enum TaskStatus TimeIntegratorTaskList::EMFCorrectReceive(MeshBlock *pmb, int step)
{
  if(pmb->pbval->ReceiveEMFCorrection() == true) {
    return TASK_NEXT;
  } else {
    return TASK_FAIL;
  }
}

//----------------------------------------------------------------------------------------
// Functions to integrate conserved variables

enum TaskStatus TimeIntegratorTaskList::HydroIntegrate(MeshBlock *pmb, int step)
{
  Hydro *ph=pmb->phydro;
  Field *pf=pmb->pfield;

  if(step == 1) {
    ph->AddFluxDivergenceToAverage(ph->u,ph->u,ph->w,pf->bcc,step_wghts[0],ph->u1);
    return TASK_NEXT;
  }

  if((step == 2) && (integrator == "vl2")) {
    ph->AddFluxDivergenceToAverage(ph->u,ph->u,ph->w1,pf->bcc1,step_wghts[1],ph->u);
    return TASK_NEXT;
  }

  if((step == 2) && (integrator == "rk2")) {
   ph->AddFluxDivergenceToAverage(ph->u,ph->u1,ph->w1,pf->bcc1,step_wghts[1],ph->u);
   return TASK_NEXT;
  }

  return TASK_FAIL;
}

enum TaskStatus TimeIntegratorTaskList::FieldIntegrate(MeshBlock *pmb, int step)
{
  if(step == 1) {
    pmb->pfield->CT(pmb->pfield->b, pmb->pfield->b, step_wghts[0], pmb->pfield->b1);
    return TASK_NEXT;
  }

  if((step == 2) && (integrator == "vl2")) {
    pmb->pfield->CT(pmb->pfield->b, pmb->pfield->b, step_wghts[1], pmb->pfield->b);
    return TASK_NEXT;
  }

  if((step == 2) && (integrator == "rk2")) {
    pmb->pfield->CT(pmb->pfield->b, pmb->pfield->b1, step_wghts[1], pmb->pfield->b);
    return TASK_NEXT;
  }

  return TASK_FAIL;
}

//----------------------------------------------------------------------------------------
// Functions to add source terms

enum TaskStatus TimeIntegratorTaskList::HydroSourceTerms(MeshBlock *pmb, int step)
{
  Hydro *ph=pmb->phydro;
  Field *pf=pmb->pfield;

  // return if there are no source terms to be added
  if (ph->psrc->hydro_sourceterms_defined == false) return TASK_NEXT;

  Real dt = (step_wghts[(step-1)].c)*(pmb->pmy_mesh->dt);
  Real time;
  // *** this must be changed for the RK3 integrator
  if(step == 1) {
    time=pmb->pmy_mesh->time;
    ph->psrc->AddHydroSourceTerms(time,dt,ph->flux,ph->w,pf->bcc,ph->u1);
  } else if(step == 2) {
    if      (integrator == "vl2") time=pmb->pmy_mesh->time + 0.5*pmb->pmy_mesh->dt;
    else if (integrator == "rk2") time=pmb->pmy_mesh->time +     pmb->pmy_mesh->dt;
    ph->psrc->AddHydroSourceTerms(time,dt,ph->flux,ph->w1,pf->bcc1,ph->u);
  } else {
    return TASK_FAIL;
  }

  return TASK_NEXT;
}

//----------------------------------------------------------------------------------------
// Functions to communicate conserved variables between MeshBlocks

enum TaskStatus TimeIntegratorTaskList::HydroSend(MeshBlock *pmb, int step)
{
  if(step == 1) {
    pmb->pbval->SendCellCenteredBoundaryBuffers(pmb->phydro->u1, 
        pmb->prad->ir1, pmb->pcr->u_cr1, pmb->phydro->ptc->u_tc1, HYDRO_CONS);
  } else if(step == 2) {
    pmb->pbval->SendCellCenteredBoundaryBuffers(pmb->phydro->u, 
           pmb->prad->ir, pmb->pcr->u_cr, pmb->phydro->ptc->u_tc, HYDRO_CONS);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::FieldSend(MeshBlock *pmb, int step)
{
  if(step == 1) {
    pmb->pbval->SendFieldBoundaryBuffers(pmb->pfield->b1);
  } else if(step == 2) {
    pmb->pbval->SendFieldBoundaryBuffers(pmb->pfield->b);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions to receive conserved variables between MeshBlocks

enum TaskStatus TimeIntegratorTaskList::HydroReceive(MeshBlock *pmb, int step)
{
  bool ret;
  if(step == 1) {
    ret=pmb->pbval->ReceiveCellCenteredBoundaryBuffers(pmb->phydro->u1, 
         pmb->prad->ir1, pmb->pcr->u_cr1, pmb->phydro->ptc->u_tc1, HYDRO_CONS);
  } else if(step == 2) {
    ret=pmb->pbval->ReceiveCellCenteredBoundaryBuffers(pmb->phydro->u, 
             pmb->prad->ir, pmb->pcr->u_cr, pmb->phydro->ptc->u_tc, HYDRO_CONS);
  } else {
    return TASK_FAIL;
  }
  if(ret==true) {
    return TASK_SUCCESS;
  } else {
    return TASK_FAIL;
  }
}

enum TaskStatus TimeIntegratorTaskList::FieldReceive(MeshBlock *pmb, int step)
{
  bool ret;
  if(step == 1) {
    ret=pmb->pbval->ReceiveFieldBoundaryBuffers(pmb->pfield->b1);
  } else if(step == 2) {
    ret=pmb->pbval->ReceiveFieldBoundaryBuffers(pmb->pfield->b);
  } else {
    return TASK_FAIL;
  }
  if(ret==true) {
    return TASK_SUCCESS;
  } else {
    return TASK_FAIL;
  }
}

//----------------------------------------------------------------------------------------
// Functions for everything else

enum TaskStatus TimeIntegratorTaskList::Prolongation(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  Radiation *prad=pmb->prad;
  CosmicRay *pcr=pmb->pcr;
  NewThermalConduction *ptc=pmb->phydro->ptc;
  BoundaryValues *pbval=pmb->pbval;
  Real dt;

  if(step == 1) {
    dt = (step_wghts[(step-1)].c)*(pmb->pmy_mesh->dt);
    pbval->ProlongateBoundaries(phydro->w1, phydro->u1, pfield->b1, pfield->bcc1,
                  prad->ir1, pcr->u_cr1, ptc->u_tc1,  pmb->pmy_mesh->time+dt, dt);
  } else if(step == 2) {
    dt=pmb->pmy_mesh->dt;
    pbval->ProlongateBoundaries(phydro->w,  phydro->u,  pfield->b,  pfield->bcc, 
                      prad->ir,  pcr->u_cr, ptc->u_tc, pmb->pmy_mesh->time+dt, dt);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::Primitives(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  if(pmb->nblevel[1][1][0]!=-1) is-=NGHOST;
  if(pmb->nblevel[1][1][2]!=-1) ie+=NGHOST;
  if(pmb->nblevel[1][0][1]!=-1) js-=NGHOST;
  if(pmb->nblevel[1][2][1]!=-1) je+=NGHOST;
  if(pmb->nblevel[0][1][1]!=-1) ks-=NGHOST;
  if(pmb->nblevel[2][1][1]!=-1) ke+=NGHOST;

  if(step == 1) {
    pmb->peos->ConservedToPrimitive(phydro->u1, phydro->w, pfield->b1,
                                    phydro->w1, pfield->bcc1, pmb->pcoord,
                                    is, ie, js, je, ks, ke);
  } else if(step == 2) {
    pmb->peos->ConservedToPrimitive(phydro->u, phydro->w1, pfield->b,
                                    phydro->w, pfield->bcc, pmb->pcoord,
                                    is, ie, js, je, ks, ke);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::PhysicalBoundary(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  CosmicRay *pcr=pmb->pcr;
  NewThermalConduction *ptc=pmb->phydro->ptc;
  BoundaryValues *pbval=pmb->pbval;
  Real dt;
  if(step == 1) {
    dt = (step_wghts[(step-1)].c)*(pmb->pmy_mesh->dt);
    pbval->ApplyPhysicalBoundaries(phydro->w1, phydro->u1, pfield->b1, pfield->bcc1,
                                pcr->u_cr1, ptc->u_tc1, pmb->pmy_mesh->time+dt, dt);
  } else if(step == 2) {
    dt=pmb->pmy_mesh->dt;
    pbval->ApplyPhysicalBoundaries(phydro->w,  phydro->u,  pfield->b,  pfield->bcc,
                                   pcr->u_cr, ptc->u_tc, pmb->pmy_mesh->time+dt, dt);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::UserWork(MeshBlock *pmb, int step)
{
  if (step != nsub_steps) return TASK_SUCCESS; // only do on last sub-step

  pmb->UserWorkInLoop();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::NewBlockTimeStep(MeshBlock *pmb, int step)
{
  if (step != nsub_steps) return TASK_SUCCESS; // only do on last sub-step

  pmb->phydro->NewBlockTimeStep();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::CheckRefinement(MeshBlock *pmb, int step)
{
  if (step != nsub_steps) return TASK_SUCCESS; // only do on last sub-step

  pmb->pmr->CheckRefinementCondition();
  return TASK_SUCCESS;
}


// Task Functions for Radiation
enum TaskStatus TimeIntegratorTaskList::RadPhysicalBoundary(MeshBlock *pmb, int step)
{
  Radiation *prad=pmb->prad;
  BoundaryValues *pbval=pmb->pbval;
  Real dt;
  if(step == 1) {
    dt = (step_wghts[(step-1)].c)*(pmb->pmy_mesh->dt);
    pbval->ApplyRadPhysicalBoundaries(prad->ir1, pmb->pmy_mesh->time+dt, dt);
  } else if(step == 2) {
    dt=pmb->pmy_mesh->dt;
    pbval->ApplyRadPhysicalBoundaries(prad->ir, pmb->pmy_mesh->time+dt, dt);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}


enum TaskStatus TimeIntegratorTaskList::RadFluxes(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Radiation *prad=pmb->prad;

  if((step == 1) && (integrator == "vl2")) {
    // copy ir to ir1
    prad->ir1 = prad->ir;
    prad->pradintegrator->CalculateFluxes(pmb, phydro->w, prad->ir, 1);
    return TASK_NEXT;
  }

  if((step == 1) && (integrator == "rk2")) {
    prad->ir1 = prad->ir;
    prad->pradintegrator->CalculateFluxes(pmb, phydro->w, prad->ir, 2);
    return TASK_NEXT;
  } 

  if(step == 2) {
    prad->pradintegrator->CalculateFluxes(pmb, phydro->w1, prad->ir1, 2);
    return TASK_NEXT;
  }

  return TASK_FAIL;
}

enum TaskStatus TimeIntegratorTaskList::RadIntegrate(MeshBlock *pmb, int step)
{
  Radiation *prad=pmb->prad;

  if(step == 1) {
    prad->pradintegrator->FluxDivergence(pmb, prad->ir, prad->ir, step_wghts[0],
                                                                      prad->ir1);
    return TASK_NEXT;
  }

  if((step == 2) && (integrator == "vl2")) {
    prad->pradintegrator->FluxDivergence(pmb, prad->ir, prad->ir, step_wghts[1], 
                                                                       prad->ir);
    return TASK_NEXT;
  }

  if((step == 2) && (integrator == "rk2")) {
   prad->pradintegrator->FluxDivergence(pmb, prad->ir, prad->ir1, step_wghts[1], 
                                                                       prad->ir);
   return TASK_NEXT;
  }

  return TASK_FAIL;
}


enum TaskStatus TimeIntegratorTaskList::RadSourceTerms(MeshBlock *pmb, int step)
{
  Hydro *ph=pmb->phydro;
  Radiation *prad=pmb->prad;

  Real dt = (step_wghts[(step-1)].c)*(pmb->pmy_mesh->dt);
  Real time;
  // *** this must be changed for the RK3 integrator
  if(step == 1) {
    prad->pradintegrator->AddSourceTerms(pmb,dt,ph->u1,ph->w,prad->ir1,1);
  } else if(step == 2) {
    if      (integrator == "vl2") time=pmb->pmy_mesh->time + 0.5*pmb->pmy_mesh->dt;
    else if (integrator == "rk2") time=pmb->pmy_mesh->time +     pmb->pmy_mesh->dt;
    prad->pradintegrator->AddSourceTerms(pmb,dt,ph->u,ph->w1,prad->ir,2);
  } else {
    return TASK_FAIL;
  }

  return TASK_NEXT;
}


enum TaskStatus TimeIntegratorTaskList::RadMomOpacity(MeshBlock *pmb, int step)
{
  Radiation *prad = pmb->prad;
  Hydro *phydro=pmb->phydro;
  
  if(step == 1) {
    prad->CalculateMoment(prad->ir1);
    prad->UpdateOpacity(pmb, phydro->w1);
  } else if(step == 2) {
    prad->CalculateMoment(prad->ir);
    prad->UpdateOpacity(pmb, phydro->w);
  } else {
    return TASK_FAIL;
  }

  return TASK_NEXT;
}


// Task functions for Cosmic Rays

enum TaskStatus TimeIntegratorTaskList::CRFluxes(MeshBlock *pmb, int step)
{
  CosmicRay *pcr=pmb->pcr;
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;

  if((step == 1) && (integrator == "vl2")) {
    // copy ir to ir1
    pcr->u_cr1 = pcr->u_cr;
    pcr->pcrintegrator->CalculateFluxes(pmb, phydro->w, pfield->bcc, 
    	                                pcr->u_cr, 1);
    return TASK_NEXT;
  }

  if((step == 1) && (integrator == "rk2")) {
    pcr->u_cr1 = pcr->u_cr;
    pcr->pcrintegrator->CalculateFluxes(pmb, phydro->w, pfield->bcc, 
    	                                pcr->u_cr, 2);
    return TASK_NEXT;
  } 

  if(step == 2) {
    pcr->pcrintegrator->CalculateFluxes(pmb, phydro->w1, pfield->bcc1, 
    	                                pcr->u_cr1, 2);
    return TASK_NEXT;
  }

  return TASK_FAIL;
}

enum TaskStatus TimeIntegratorTaskList::CRIntegrate(MeshBlock *pmb, int step)
{
  CosmicRay *pcr=pmb->pcr;
  Hydro *ph=pmb->phydro;
  Field *pfield=pmb->pfield;

  if(step == 1) {
    pcr->pcrintegrator->FluxDivergence(pmb, pcr->u_cr, pcr->u_cr, step_wghts[0],
                                       pcr->u_cr1, ph->u1,ph->w,pfield->bcc);
    return TASK_NEXT;
  }

  if((step == 2) && (integrator == "vl2")) {
    pcr->pcrintegrator->FluxDivergence(pmb, pcr->u_cr, pcr->u_cr, step_wghts[1], 
                                       pcr->u_cr, ph->u, ph->w1, pfield->bcc1);
    return TASK_NEXT;
  }

  if((step == 2) && (integrator == "rk2")) {
    pcr->pcrintegrator->FluxDivergence(pmb, pcr->u_cr, pcr->u_cr1, step_wghts[1], 
                                       pcr->u_cr, ph->u, ph->w1, pfield->bcc1);
    return TASK_NEXT;
  }

  return TASK_FAIL;
}


enum TaskStatus TimeIntegratorTaskList::CRSourceTerms(MeshBlock *pmb, int step)
{
  CosmicRay *pcr=pmb->pcr;
  Hydro *ph=pmb->phydro;
  Field *pfield=pmb->pfield;


  Real dt = (step_wghts[(step-1)].c)*(pmb->pmy_mesh->dt);
  Real time;
  // *** this must be changed for the RK3 integrator
  if(step == 1) {
    pcr->pcrintegrator->AddSourceTerms(pmb,dt,ph->u1,ph->w,pfield->bcc,
                                                           pcr->u_cr1,1);
  } else if(step == 2) {
    if      (integrator == "vl2") time=pmb->pmy_mesh->time + 0.5*pmb->pmy_mesh->dt;
    else if (integrator == "rk2") time=pmb->pmy_mesh->time +     pmb->pmy_mesh->dt;
    pcr->pcrintegrator->AddSourceTerms(pmb,dt,ph->u,ph->w1,pfield->bcc1,
                                                             pcr->u_cr,2);
  } else {
    return TASK_FAIL;
  }

  return TASK_NEXT;
}


enum TaskStatus TimeIntegratorTaskList::CRVAOpacity(MeshBlock *pmb, int step)
{
  CosmicRay *pcr = pmb->pcr;
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;

  Real dt = (step_wghts[(step-1)].c)*(pmb->pmy_mesh->dt);
  
  if(step == 1) {
    // Need to update the Eddington tensor first to get cosmic ray pressure
    pcr->UpdateCRTensor(pmb, phydro->w1);
    pcr->UpdateDiff(pmb, pcr->u_cr1,phydro->w1,pfield->bcc1,dt);
  } else if(step == 2) {
    pcr->UpdateCRTensor(pmb, phydro->w);
    pcr->UpdateDiff(pmb, pcr->u_cr,phydro->w,pfield->bcc,dt);
  } else {
    return TASK_FAIL;
  }

  return TASK_NEXT;

}


// tasks for thermal conduction


enum TaskStatus TimeIntegratorTaskList::TCInitialize(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  NewThermalConduction *ptc=phydro->ptc;

  if(step == 1) {
    ptc->Initialize(pmb, phydro->w, ptc->u_tc);
    ptc->u_tc1 = ptc->u_tc;
    return TASK_NEXT;
  }

  if(step == 2) {
    ptc->Initialize(pmb, phydro->w1, ptc->u_tc1);

    return TASK_NEXT;
  }

  return TASK_FAIL;
}



enum TaskStatus TimeIntegratorTaskList::TCFluxes(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  NewThermalConduction *ptc=phydro->ptc;

  if((step == 1) && (integrator == "vl2")) {
    ptc->pntcintegrator->CalculateFluxes(pmb, phydro->w, pfield->bcc, 
                                      ptc->u_tc, 1);
    return TASK_NEXT;
  }

  if((step == 1) && (integrator == "rk2")) {
    ptc->pntcintegrator->CalculateFluxes(pmb, phydro->w, pfield->bcc, 
                                      ptc->u_tc, 2);
    return TASK_NEXT;
  } 

  if(step == 2) {
    ptc->pntcintegrator->CalculateFluxes(pmb, phydro->w1, pfield->bcc1, 
                                      ptc->u_tc1, 2);
    return TASK_NEXT;
  }

  return TASK_FAIL;
}


enum TaskStatus TimeIntegratorTaskList::TCIntegrate(MeshBlock *pmb, int step)
{

  Hydro *ph=pmb->phydro;
  Field *pfield=pmb->pfield;
  NewThermalConduction *ptc=ph->ptc;

  if(step == 1) {
    ptc->pntcintegrator->FluxDivergence(pmb, ph->w, ptc->u_tc, ptc->u_tc, 
                                       step_wghts[0], ptc->u_tc1);
    return TASK_NEXT;
  }

  if((step == 2) && (integrator == "vl2")) {
    ptc->pntcintegrator->FluxDivergence(pmb, ph->w1, ptc->u_tc, ptc->u_tc, 
                                       step_wghts[1], ptc->u_tc);
    return TASK_NEXT;
  }

  if((step == 2) && (integrator == "rk2")) {
    ptc->pntcintegrator->FluxDivergence(pmb, ph->w1, ptc->u_tc, ptc->u_tc1, 
                                       step_wghts[1], ptc->u_tc);
    return TASK_NEXT;
  }

  return TASK_FAIL;
}

enum TaskStatus TimeIntegratorTaskList::TCSourceTerms(MeshBlock *pmb, int step)
{

  Hydro *ph=pmb->phydro;
  Field *pfield=pmb->pfield;
  NewThermalConduction *ptc=ph->ptc;

  Real dt = (step_wghts[(step-1)].c)*(pmb->pmy_mesh->dt);
  Real time;
  // *** this must be changed for the RK3 integrator
  if(step == 1) {
    ptc->pntcintegrator->AddSourceTerms(pmb,dt,ph->u1,ptc->u_tc1,1);
  } else if(step == 2) {
    if      (integrator == "vl2") time=pmb->pmy_mesh->time + 0.5*pmb->pmy_mesh->dt;
    else if (integrator == "rk2") time=pmb->pmy_mesh->time +     pmb->pmy_mesh->dt;
    ptc->pntcintegrator->AddSourceTerms(pmb,dt,ph->u,ptc->u_tc,2);
  } else {
    return TASK_FAIL;
  }

  return TASK_NEXT;
}


enum TaskStatus TimeIntegratorTaskList::TCOpacity(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  NewThermalConduction *ptc = phydro->ptc;

  Real dt = (step_wghts[(step-1)].c)*(pmb->pmy_mesh->dt);
  
  if(step == 1) {
    ptc->UpdateKappa(pmb, phydro->w1, pfield->bcc1);
  } else if(step == 2) {
    ptc->UpdateKappa(pmb, phydro->w,pfield->bcc);
  } else {
    return TASK_FAIL;
  }

  return TASK_NEXT;

}

