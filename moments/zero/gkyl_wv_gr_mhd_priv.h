#pragma once

// Private header, not for direct use in user-facing code.

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct wv_gr_mhd {
  struct gkyl_wv_eqn eqn; // Base equation object.
  struct gkyl_gr_spacetime *spacetime; // Pointer to base spacetime object.
  double gas_gamma; // Adiabatic index.

  enum gkyl_spacetime_gauge spacetime_gauge; // Spacetime gauge choice.
  int reinit_freq; // Spacetime reinitialization frequency.
};

/**
* Compute flux vector. Assumes rotation to local coordinate system.
*
* @param gas_gamma Adiabatic index.
* @param q Conserved variable vector.
* @param flux Flux vector in direction 'dir' (output).
*/
GKYL_CU_D
void
gkyl_gr_mhd_flux(double gas_gamma, const double q[74], double flux[74]);

/**
* Compute primitive variables given the conserved variables.
*
* @param gas_gamma Adiabatic index.
* @param q Conserved variable vector.
* @param v Primitive variable vector (output).
*/
GKYL_CU_D
void
gkyl_gr_mhd_prim_vars(double gas_gamma, const double q[74], double v[74]);

/**
* Compute inverse spatial metric tensor (in covariant component form) given the conserved variables.
*
* @param q Conserved variable vector.
* @param inv_spatial_metric Inverse spatial metric tensor (output).
*/
GKYL_CU_D
void
gkyl_gr_mhd_inv_spatial_metric(const double q[74], double ***inv_spatial_metric);

/**
* Compute ideal magnetohydrodynamic stress-energy tensor (in contravariant component form) given the conserved variables.
*
* @param gas_gamma Adiabatic index.
* @param q Conserved variable vector.
* @param stress_energy Stress-energy tensor (output).
*/
GKYL_CU_D
void
gkyl_gr_mhd_stress_energy_tensor(double gas_gamma, const double q[74], double ***stress_energy);

/**
* Compute Riemann variables given the conserved variables.
*
* @param eqn Base equation object.
* @param qstate Current state vector.
* @param qin Conserved variable vector (input).
* @param wout Riemann variable vector (output).
*/
GKYL_CU_D
static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout);

/**
* Compute conserved variables given the Riemann variables.
*
* @param eqn Base equation object.
* @param qstate Current state vector.
* @param win Riemann variable vector (input).
* @param qout Conserved variable vector (output).
*/
GKYL_CU_D
static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double *qout);

/**
* Rotate state vector from global to local coordinate frame.
*
* @param eqn Base equation object.
* @param tau1 First tangent vector of the coordinate frame.
* @param tau2 Second tangent vector of the coordinate frame.
* @param norm Normal vector of the coordinate frame.
* @param qglobal State vector in global coordinate frame (input).
* @param qlocal State vector in local coordinate frame (output).
*/
GKYL_CU_D
static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal,
  double* GKYL_RESTRICT qlocal);

/**
* Rotate state vector from local to global coordinate frame.
*
* @param eqn Base equation object.
* @param tau1 First tangent vector of the coordinate frame.
* @param tau2 Second tangent vector of the coordinate frame.
* @param norm Normal vector of the coordinate frame.
* @param qlocal State vector in local coordinate frame (input).
* @param qglobal State vector in global coordinate frame (output).
*/
GKYL_CU_D
static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal,
  double* GKYL_RESTRICT qglobal);

/**
* Free general relativistic magnetohydrodynamics equations object with ideal gas equation of state.
*
* @param ref Reference counter for general relativistic magnetohydrodynamics equations with ideal gas equation of state.
*/
void gkyl_gr_mhd_free(const struct gkyl_ref_count* ref);