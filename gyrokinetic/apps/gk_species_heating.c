#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static double
gk_heating_volume_integrate(gkyl_gyrokinetic_app *app, struct gk_heating *src, const struct gkyl_array *arrin)
{
  // Compute the volume integral of arrin.
  gkyl_array_integrate_advance(src->vol_integ_op, arrin, 1.0,
    0, &app->local, 0, src->volint_local);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, src->volint_local, src->volint_global);
  double volint_global = 0.0;
  if (app->use_gpu)
    gkyl_cu_memcpy(&volint_global, src->volint_global, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(&volint_global, src->volint_global, sizeof(double));
  return volint_global;
}

void
gk_species_heating_rhs_disabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_heating *src, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
}

void
gk_species_heating_rhs_enabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_heating *src, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  // Compute Maxwellian moments (n, u_par, T/m).
  gk_species_moment_calc(&species->lte.moms, species->local, app->local, fin);
  gkyl_dg_div_op_range(species->lte.moms.mem_geo, app->basis, 0, species->lte.moms.marr, 
    0, species->lte.moms.marr, 0, app->gk_geom->jacobgeo, &app->local);  

  // Volume integrate Jrate times the thermal M2.
  gkyl_dg_mul_op_range(app->basis, 0, src->Jrate_mom, 0, src->Jrate, 0, species->lte.moms.marr, &app->local);
  gkyl_dg_mul_op_range(app->basis, 0, src->Jrate_mom, 0, src->Jrate_mom, 2, species->lte.moms.marr, &app->local);
  double Jrate_M2thermal_int = GKYL_MAX2(0.0, gk_heating_volume_integrate(app, src, src->Jrate_mom));

  // Volume integrate Jrate times the vtsq_shape time M0.
  gkyl_dg_mul_op_range(app->basis, 0, src->Jrate_mom, 0, src->Jrate_vtsq_shape, 0, species->lte.moms.marr, &app->local);
  double Jrate_vtsq_shape_M0_int = GKYL_MAX2(0.0, gk_heating_volume_integrate(app, src, src->Jrate_mom));

  // Thermal speed squared of the Maxwellian.
  double vtsq_amplitude = (src->norm_power + Jrate_M2thermal_int)/Jrate_vtsq_shape_M0_int;
  gkyl_array_set_offset_range(species->lte.moms.marr, vtsq_amplitude, src->vtsq_shape, 2*app->basis.num_basis, &app->local);

  // Compute the Maxwellian.
  gk_species_lte_from_moms(app, species, &species->lte, species->lte.moms.marr);

  // Multiply the Maxwellian by Jrate.
  gkyl_dg_mul_conf_phase_op_range(&app->basis, &species->basis, src->Jrate_fmax, 
    src->Jrate, species->lte.f_lte, &app->local, &species->local);

  // Assemble the BGK-like term and add it to rhs.
  gkyl_bgk_collisions_advance(src->bgk_op, &app->local, &species->local, 
    src->rate, src->Jrate_fmax, fin, src->implicit_step, src->dt_implicit, rhs, species->cflrate);

  app->stat.species_heat_tm += gkyl_time_diff_now_sec(wst);
}

void 
gk_species_heating_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_heating *src)
{
  src->heating_id = s->info.heating.heating_id;
  src->write_diagnostics = s->info.heating.write_diagnostics;

  if (src->heating_id) {
    src->norm_power = 2.0*s->info.heating.power/((app->vdim == 1? 1.0 : 3.0)*s->info.mass);

    // Heating rate.
    src->rate = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    struct gkyl_array *rate_host = app->use_gpu? mkarr(false, src->rate->ncomp, src->rate->size)
                                               : gkyl_array_acquire(src->rate);
    gkyl_proj_on_basis *proj_rate = gkyl_proj_on_basis_new(&app->grid, &app->basis,
      app->poly_order+1, 1, s->info.heating.rate_profile, s->info.heating.rate_profile_ctx);
    gkyl_proj_on_basis_advance(proj_rate, 0.0, &app->local, rate_host);
    gkyl_array_copy(src->rate, rate_host);
    gkyl_proj_on_basis_release(proj_rate);
    gkyl_array_release(rate_host);
    // Multiply the rate by the conf-space Jacobian.
    src->Jrate = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkyl_dg_mul_op_range(app->basis, 0, src->Jrate, 0, app->gk_geom->jacobgeo, 0, src->rate, &app->local);

    // Heating rate.
    src->vtsq_shape = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    struct gkyl_array *vtsq_shape_host = app->use_gpu? mkarr(false, src->vtsq_shape->ncomp, src->vtsq_shape->size)
                                                     : gkyl_array_acquire(src->vtsq_shape);
    gkyl_proj_on_basis *proj_vtsq_shape = gkyl_proj_on_basis_new(&app->grid, &app->basis,
      app->poly_order+1, 1, s->info.heating.temp_shape, s->info.heating.temp_shape_ctx);
    gkyl_proj_on_basis_advance(proj_vtsq_shape, 0.0, &app->local, vtsq_shape_host);
    gkyl_array_copy(src->vtsq_shape, vtsq_shape_host);
    gkyl_proj_on_basis_release(proj_vtsq_shape);
    gkyl_array_release(vtsq_shape_host);

    // Multiply Jrate by the shape of v_t^2.
    src->Jrate_vtsq_shape = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkyl_dg_mul_op_range(app->basis, 0, src->Jrate_vtsq_shape, 0, src->Jrate, 0, src->vtsq_shape, &app->local);

    // Rate times the Maxwellian.
    src->Jrate_fmax = mkarr(app->use_gpu, s->basis.num_basis, s->local_ext.volume);
    // Rate times a velocity moment.
    src->Jrate_mom = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

    // Volume integrator.
    src->vol_integ_op = gkyl_array_integrate_new(&app->grid, &app->basis, 1, GKYL_ARRAY_INTEGRATE_OP_NONE, app->use_gpu);
    if (app->use_gpu) {
      src->volint_local = gkyl_cu_malloc(sizeof(double));
      src->volint_global = gkyl_cu_malloc(sizeof(double));
    } 
    else {
      src->volint_local = gkyl_malloc(sizeof(double));
      src->volint_global = gkyl_malloc(sizeof(double));
    }

    // BGK operator.
    src->bgk_op = gkyl_bgk_collisions_new(&app->basis, &s->basis, app->use_gpu);
    src->implicit_step = false;
    src->dt_implicit = 1e9;

    // Methods chosen at runtime.
    src->rhs_func = gk_species_heating_rhs_enabled;
  }
  else {
    src->rhs_func = gk_species_heating_rhs_disabled;
  }
}

void
gk_species_heating_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_heating *src, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  src->rhs_func(app, species, src, fin, rhs);
}

void
gk_species_heating_release(const struct gkyl_gyrokinetic_app *app, const struct gk_heating *src)
{
  if (src->heating_id) {
    gkyl_array_release(src->rate);
    gkyl_array_release(src->Jrate);
    gkyl_array_release(src->vtsq_shape);
    gkyl_array_release(src->Jrate_vtsq_shape);
    gkyl_array_release(src->Jrate_fmax);
    gkyl_array_release(src->Jrate_mom);

    gkyl_array_integrate_release(src->vol_integ_op);
    if (app->use_gpu) {
      gkyl_cu_free(src->volint_local);
      gkyl_cu_free(src->volint_global);
    }
    else {
      gkyl_free(src->volint_local);
      gkyl_free(src->volint_global);
    }

    gkyl_bgk_collisions_release(src->bgk_op);
  }
}
