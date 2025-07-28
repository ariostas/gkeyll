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
* Free general relativistic magnetohydrodynamics equations object with ideal gas equation of state.
*
* @param ref Reference counter for general relativistic magnetohydrodynamics equations with ideal gas equation of state.
*/
void gkyl_gr_mhd_free(const struct gkyl_ref_count* ref);