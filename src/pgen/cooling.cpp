// C/C++ headers
#include <cmath>
#include <stdio.h>
#include <string>
#include <cstdlib>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"
#include "../cr/cr.hpp"
#include "../cr/integrators/cr_integrators.hpp"

// Cooling Tables
static int nbins;
static Real const_factor;
static AthenaArray<Real> cool_t;
static AthenaArray<Real> cool_tef;
static AthenaArray<Real> cool_coef;
static AthenaArray<Real> cool_index;

// Define Code units in cgs
static const Real mh = 1.6605e-24;   // atomic mass unit (g)
static const Real kb = 1.380648e-16; // boltzmann constant (erg/K)

static const Real unit_len  = 3.086e18; // 1 pc
static const Real unit_temp = 1.0e5;
static const Real unit_n    = 1.0e-1;

static const Real unit_rho  = unit_n * mh;
static const Real unit_pres = kb * unit_n * unit_temp;
static const Real unit_vel  = sqrt(unit_pres / unit_rho);
static const Real unit_time = unit_len / unit_vel;

//------------------------------------------------------------------------------
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // Read/set problem parameters
  Real gm1 = peos->GetGamma() - 1.0;
  Real amp = pin->GetReal("problem","amp");
  Real rho = pin->GetReal("problem","rho");
  Real temp = pin->GetReal("problem","temp");

  // Initialize Variables
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    Real x = pcoord->x1v(i);
    Real y = pcoord->x2v(j);

    Real sigma = 0.2;
    //Real perturb = amp*exp(-(SQR(x-0.5) + SQR(y-0.5))/SQR(sigma)); // Gaussian
    Real perturb = amp*((double)rand()/(double)RAND_MAX-0.5); // Random

    phydro->u(IDN,k,j,i) = rho*(1.0 + perturb);

    phydro->u(IM1,k,j,i) = 0.0;
    phydro->u(IM2,k,j,i) = 0.0;
    phydro->u(IM3,k,j,i) = 0.0;

    phydro->u(IEN,k,j,i) = (temp*rho)/gm1;
    phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i)) +
                                 SQR(phydro->u(IM2,k,j,i)) +
                                 SQR(phydro->u(IM3,k,j,i)))
                               /phydro->u(IDN,k,j,i);
  }}}

  // initialize uniform interface B
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie+1; i++) {
      pfield->b.x1f(k,j,i) = 0.707106781;// pin->GetReal("problem","bx");
    }}}
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je+1; j++) {
    for (int i=is; i<=ie; i++) {
      pfield->b.x2f(k,j,i) = 0.707106781;//pin->GetReal("problem","by");
    }}}
    for (int k=ks; k<=ke+1; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      pfield->b.x3f(k,j,i) = 0.0;
    }}}
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      phydro->u(IEN,k,j,i) += 0.5*(SQR(pfield->b.x1f(k,j,i)) +
                                   SQR(pfield->b.x2f(k,j,i)) +
                                   SQR(pfield->b.x3f(k,j,i)));
    }}}
  }

  return;
}

// Exact Integration Scheme for Radiative Cooling from Townsend (2009)
// Returns: Temperature(K) at the next timestep after cooling
// Requires: - Input temperature, density, and timestep in cgs units
//           - T < cool_t(nbins-1) and T > cool_t(0)
//           - All cooling slopes are not equal to 1
Real townsend(Real temp, Real rho, const Real dt)
{
  // Set a temperature floor for Cooling
  Real t_floor = 1.0e4;

  // Check that temperature is above the cooling floor
  if (temp > t_floor) {
    // Get Reference values from the last bin
    Real t_n    = cool_t(nbins-1);
    Real coef_n = cool_coef(nbins-1);

    // Get the index of the right temperature bin
    int idx = 0;
    while ((idx < nbins-2) && (cool_t(idx+1) < temp)) { idx += 1; }

    // Look up the corresponding slope and coefficient
    Real t_i   = cool_t(idx);
    Real tef_i = cool_tef(idx);
    Real coef  = cool_coef(idx);
    Real slope = cool_index(idx);

    // Compute the Temporal Evolution Function Y(T) (Eq. A5)
    Real sm1 = slope - 1.0;
    Real tef = tef_i + (coef_n/coef)*(t_i/t_n)*(pow(t_i/temp,sm1)-1)/sm1;

    // Compute the adjusted TEF for new timestep (Eqn. 26)
    Real tef_adj = tef + rho*coef_n*const_factor*dt/t_n;

    // TEF is a strictly decreasing function and new_tef > tef
    // Check if the new TEF falls into a lower bin
    while ((idx > 0) && (tef_adj > cool_tef(idx))) {
      idx -= 1;
      // If so, update slopes and coefficients
      t_i   = cool_t(idx);
      tef_i = cool_tef(idx);
      coef  = cool_coef(idx);
      slope = cool_index(idx);
    }

    // Compute the Inverse Temporal Evolution Function Y^{-1}(Y) (Eq. A7)
    Real oms  = 1.0 - slope;
    Real itef = t_i*pow(1-oms*(coef/coef_n)*(t_n/t_i)*(tef_adj-tef_i),1/oms);

    // Return the new temperature if it is still above the temperature floor
    return std::max(itef,t_floor);
  }
  else { return t_floor; }
}

// User Defined Cooling Function
void cooling(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
             AthenaArray<Real> &cons)
{
  Real g = pmb->peos->GetGamma();

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        // Need to take density and temperature at time step n from cons, not
        // prim because we do not need intermediate step to calculate the
        // cooling function
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

        // T = P/rho
        Real temp = eint * (g-1.0)/rho;

        // Calculate new temperature using the Townsend Algorithm
        // The inputs should be given in cgs units
        Real temp_cgs = temp * unit_temp;
        Real rho_cgs  = rho  * unit_rho;
        Real dt_cgs   = dt   * unit_time;
        Real temp_new = townsend(temp_cgs,rho_cgs,dt_cgs)/unit_temp;

        // Update energy based on change in temperature
        cons(IEN,k,j,i) += (temp_new - temp) * (rho/(g-1.0));

        // Store the change in energy in a user defined output variable
        pmb->user_out_var(0,k,j,i) = (temp_new - temp) * rho/(g-1.0);

        // Update energy based on a constant heating term
        Real heating_rate = 0.1;
        if (temp_new * unit_temp < 1.0e6) {
          cons(IEN,k,j,i) += heating_rate * (dt/(g-1.0));
        }
      }
    }
  }

  return;
}

// Register our user defined cooling function with Athena
void Mesh::InitUserMeshData(ParameterInput *pin)
{
  EnrollUserExplicitSourceFunction(cooling);

  // Initialize all values in the cooling table
  nbins = 7;
  cool_t.NewAthenaArray(nbins);
  cool_tef.NewAthenaArray(nbins);
  cool_coef.NewAthenaArray(nbins);
  cool_index.NewAthenaArray(nbins);

  // Last Temperature should be the ceiling temperature
  // [K]
  cool_t(0) = 10.000000;
  cool_t(1) = 138.949549;
  cool_t(2) = 1930.697729;
  cool_t(3) = 26826.957953;
  cool_t(4) = 372759.372031;
  cool_t(5) = 5179474.679231;
  cool_t(6) = 71968567.300115;

  // [1e-23 ergs*cm3/s]
  cool_coef(0) = 0.000062;
  cool_coef(1) = 0.017520;
  cool_coef(2) = 0.028247;
  cool_coef(3) = 5.059731;
  cool_coef(4) = 18.207846;
  cool_coef(5) = 1.302815;
  cool_coef(6) = 1.592144;

  cool_index(0) = 2.145482;
  cool_index(1) = 0.181502;
  cool_index(2) = 1.971516;
  cool_index(3) = 0.486615;
  cool_index(4) = -1.002204;
  cool_index(5) = 0.076212;
  cool_index(6) = 0.167512;

  // Calculate All TEFs Y_k recursively (Eq. A6)
  cool_tef(nbins-1) = 0.0; // Last Y_N = 0
  for (int i=nbins-2; i>=0; i--) {
    Real t_n    = cool_t(nbins-1);
    Real coef_n = cool_coef(nbins-1);

    Real t_i   = cool_t(i);
    Real tef_i = cool_tef(i);
    Real coef  = cool_coef(i);
    Real slope = cool_index(i);

    Real sm1  = slope - 1.0;
    Real step = (coef_n/coef)*(t_i/t_n)*(pow(t_i/cool_t(i+1),sm1)-1)/sm1;
    cool_tef(i) = cool_tef(i+1) - step;
  }

  // Compute the constant factor needed to compute new TEFs
  Real g  = 5.0/3.0; // adiabatic index
  Real X = 0.7; Real Z = 0.02; // Fully ionized, solar abundances
  Real mu   = 1.0/(2.0*X + 0.75*(1.0-X-Z) + Z/2.0);
  Real mu_e = 1.0/(1.0+X);
  Real mu_h = 1.0/X;

  const_factor = (1.0e-23)*(g-1.0)*mu/(kb*mu_e*mu_h*mh);
}

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  cool_t.DeleteAthenaArray();
  cool_tef.DeleteAthenaArray();
  cool_coef.DeleteAthenaArray();
  cool_index.DeleteAthenaArray();
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(1);
}
