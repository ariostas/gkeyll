#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_gr_mhd.h>
#include <gkyl_wv_gr_mhd_priv.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_gr_blackhole.h>

void
test_gr_mhd_basic_minkowski()
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  struct gkyl_wv_eqn *gr_mhd = gkyl_wv_gr_mhd_new(gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  TEST_CHECK( gr_mhd->num_equations == 74 );
  // TEST_CHECK( gr_mhd->num_waves == 2 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, p = 1.5;
      double mag_x = 0.3, mag_y = 0.2, mag_z = 0.1;

      double spatial_det, lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der = gkyl_malloc(sizeof(double[3]));
      double **shift_der = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_det);
      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);
      
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

      spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
      spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

      double *vel = gkyl_malloc(sizeof(double[3]));
      vel[0] = u; vel[1] = v; vel[2] = w;

      double v_sq = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq += spatial_metric[i][j] * vel[i] * vel[j];
        }
      }

      double W = 1.0 / sqrt(1.0 - v_sq);

      double *mag = gkyl_malloc(sizeof(double[3]));
      mag[0] = mag_x; mag[1] = mag_y; mag[2] = mag_z;

      double *cov_mag = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        cov_mag[i] = 0.0;

        for (int j = 0; j < 3; j++) {
          cov_mag[i] += spatial_metric[i][j] * mag[j];
        }
      }

      double *cov_vel = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        cov_vel[i] = 0.0;

        for (int j = 0; j < 3; j++) {
          cov_vel[i] += spatial_metric[i][j] * vel[j];
        }
      }
      
      double b0 = 0.0;
      for (int i = 0; i < 3; i++) {
        b0 += W * mag[i] * (cov_vel[i] / lapse);
      }

      double *spacetime_vel = gkyl_malloc(sizeof(double[4]));
      spacetime_vel[0] = W / lapse;
      for (int i = 0; i < 3; i++) {
        spacetime_vel[i + 1] = (W * vel[i]) - (shift[i] * (W / lapse));
      }

      double *b = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        b[i] = (mag[i] + (lapse * b0 * spacetime_vel[i + 1])) / W;
      }

      double b_sq = 0.0;
      for (int i = 0; i < 3; i++) {
        b_sq += (mag[i] * cov_mag[i]) / (W * W);
      }
      b_sq += ((lapse * lapse) * (b0 * b0)) / (W * W);

      double *cov_b = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        cov_b[i] = 0.0;

        for (int j = 0; j < 3; j++) {
          cov_b[i] += spatial_metric[i][j] * b[j];
        }
      }

      double h_star = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0))) + (b_sq / rho);
      double p_star = p + (0.5 * b_sq);

      double q[74];
      q[0] = sqrt(spatial_det) * rho * W;
      q[1] = sqrt(spatial_det) * ((rho * h_star * (W * W) * cov_vel[0]) - (lapse * b0 * cov_b[0]));
      q[2] = sqrt(spatial_det) * ((rho * h_star * (W * W) * cov_vel[1]) - (lapse * b0 * cov_b[1]));
      q[3] = sqrt(spatial_det) * ((rho * h_star * (W * W) * cov_vel[2]) - (lapse * b0 * cov_b[2]));
      q[4] = sqrt(spatial_det) * ((rho * h_star * (W * W)) - p_star - ((lapse * lapse) * (b0 * b0)) - (rho * W));

      q[5] = sqrt(spatial_det) * mag_x;
      q[6] = sqrt(spatial_det) * mag_y;
      q[7] = sqrt(spatial_det) * mag_z;

      q[8] = lapse;
      q[9] = shift[0]; q[10] = shift[1]; q[11] = shift[2];

      q[12] = spatial_metric[0][0]; q[13] = spatial_metric[0][1]; q[14] = spatial_metric[0][2];
      q[15] = spatial_metric[1][0]; q[16] = spatial_metric[1][1]; q[17] = spatial_metric[1][2];
      q[18] = spatial_metric[2][0]; q[19] = spatial_metric[2][1]; q[20] = spatial_metric[2][2];

      q[21] = extrinsic_curvature[0][0]; q[22] = extrinsic_curvature[0][1]; q[23] = extrinsic_curvature[0][2];
      q[24] = extrinsic_curvature[1][0]; q[25] = extrinsic_curvature[1][1]; q[26] = extrinsic_curvature[1][2];
      q[27] = extrinsic_curvature[2][0]; q[28] = extrinsic_curvature[2][1]; q[29] = extrinsic_curvature[2][2];

      q[30] = 1.0;

      q[31] = lapse_der[0]; q[32] = lapse_der[1]; q[33] = lapse_der[2];
      q[34] = shift_der[0][0]; q[35] = shift_der[0][1]; q[36] = shift_der[0][2];
      q[37] = shift_der[1][0]; q[38] = shift_der[1][1]; q[39] = shift_der[1][2];
      q[40] = shift_der[2][0]; q[41] = shift_der[2][1]; q[42] = shift_der[2][2];

      q[43] = spatial_metric_der[0][0][0]; q[44] = spatial_metric_der[0][0][1]; q[45] = spatial_metric_der[0][0][2];
      q[46] = spatial_metric_der[0][1][0]; q[47] = spatial_metric_der[0][1][1]; q[48] = spatial_metric_der[0][1][2];
      q[49] = spatial_metric_der[0][2][0]; q[50] = spatial_metric_der[0][2][1]; q[51] = spatial_metric_der[0][2][2];

      q[52] = spatial_metric_der[1][0][0]; q[53] = spatial_metric_der[1][0][1]; q[54] = spatial_metric_der[1][0][2];
      q[55] = spatial_metric_der[1][1][0]; q[56] = spatial_metric_der[1][1][1]; q[57] = spatial_metric_der[1][1][2];
      q[58] = spatial_metric_der[1][2][0]; q[59] = spatial_metric_der[1][2][1]; q[60] = spatial_metric_der[1][2][2];

      q[61] = spatial_metric_der[2][0][0]; q[62] = spatial_metric_der[2][0][1]; q[63] = spatial_metric_der[2][0][2];
      q[64] = spatial_metric_der[2][1][0]; q[65] = spatial_metric_der[2][1][1]; q[66] = spatial_metric_der[2][1][2];
      q[67] = spatial_metric_der[2][2][0]; q[68] = spatial_metric_der[2][2][1]; q[69] = spatial_metric_der[2][2][2];

      q[70] = 0.0;
      q[71] = x; q[72] = y; q[73] = 0.0;

      double prims[74];
      gkyl_gr_mhd_prim_vars(gas_gamma, q, prims);
      
      TEST_CHECK( gkyl_compare(prims[0], rho, 1e-1) );
      TEST_CHECK( gkyl_compare(prims[1], u, 1e-1) );
      TEST_CHECK( gkyl_compare(prims[2], v, 1e-1) );
      TEST_CHECK( gkyl_compare(prims[3], w, 1e-1) );
      TEST_CHECK( gkyl_compare(prims[4], p, 1e-1) );

      TEST_CHECK( gkyl_compare(prims[5], mag_x, 1e-1) );
      TEST_CHECK( gkyl_compare(prims[6], mag_y, 1e-1) );
      TEST_CHECK( gkyl_compare(prims[7], mag_z, 1e-1) );

      double D = rho * W;
      double Sx = (rho * h_star * (W * W) * cov_vel[0]) - (lapse * b0 * cov_b[0]);
      double Sy = (rho * h_star * (W * W) * cov_vel[1]) - (lapse * b0 * cov_b[1]);
      double Sz = (rho * h_star * (W * W) * cov_vel[2]) - (lapse * b0 * cov_b[2]);
      double Etot = (rho * h_star * (W * W)) - p_star - ((lapse * lapse) * (b0 * b0)) - (rho * W);

      double fluxes[3][8] = {
        { (lapse * sqrt(spatial_det)) * (D * (vel[0] - (shift[0] / lapse))),
          (lapse * sqrt(spatial_det)) * ((Sx * (vel[0] - (shift[0] / lapse))) + p_star - ((cov_b[0] * mag[0]) / W)),
          (lapse * sqrt(spatial_det)) * ((Sy * (vel[0] - (shift[0] / lapse))) - ((cov_b[1] * mag[0]) / W)),
          (lapse * sqrt(spatial_det)) * ((Sz * (vel[0] - (shift[0] / lapse))) - ((cov_b[2] * mag[0]) / W)),
          (lapse * sqrt(spatial_det)) * ((Etot * (vel[0] - (shift[0] / lapse))) + (p_star * vel[0]) - ((lapse * b0 * mag[0]) / W)),
          (lapse * sqrt(spatial_det)) * (((vel[0] - (shift[0] / lapse)) * mag[0]) - ((vel[0] - (shift[0] / lapse)) * mag[0])),
          (lapse * sqrt(spatial_det)) * (((vel[0] - (shift[0] / lapse)) * mag[1]) - ((vel[1] - (shift[1] / lapse)) * mag[0])),
          (lapse * sqrt(spatial_det)) * (((vel[0] - (shift[0] / lapse)) * mag[2]) - ((vel[2] - (shift[2] / lapse)) * mag[0])) },
        { (lapse * sqrt(spatial_det)) * (D * (vel[1] - (shift[1] / lapse))),
          (lapse * sqrt(spatial_det)) * ((Sx * (vel[1] - (shift[1] / lapse))) - ((cov_b[0] * mag[1]) / W)),
          (lapse * sqrt(spatial_det)) * ((Sy * (vel[1] - (shift[1] / lapse))) + p_star - ((cov_b[1] * mag[1]) / W)),
          (lapse * sqrt(spatial_det)) * ((Sz * (vel[1] - (shift[1] / lapse))) - ((cov_b[2] * mag[1]) / W)),
          (lapse * sqrt(spatial_det)) * ((Etot * (vel[1] - (shift[1] / lapse))) + (p_star * vel[1]) - ((lapse * b0 * mag[1]) / W)),
          (lapse * sqrt(spatial_det)) * (((vel[1] - (shift[1] / lapse)) * mag[0]) - ((vel[0] - (shift[0] / lapse)) * mag[1])),
          (lapse * sqrt(spatial_det)) * (((vel[1] - (shift[1] / lapse)) * mag[1]) - ((vel[1] - (shift[1] / lapse)) * mag[1])),
          (lapse * sqrt(spatial_det)) * (((vel[1] - (shift[1] / lapse)) * mag[2]) - ((vel[2] - (shift[2] / lapse)) * mag[1])) },
        { (lapse * sqrt(spatial_det)) * (D * (vel[2] - (shift[2] / lapse))),
          (lapse * sqrt(spatial_det)) * ((Sx * (vel[2] - (shift[2] / lapse))) - ((cov_b[0] * mag[2]) / W)),
          (lapse * sqrt(spatial_det)) * ((Sy * (vel[2] - (shift[2] / lapse))) - ((cov_b[1] * mag[2]) / W)),
          (lapse * sqrt(spatial_det)) * ((Sz * (vel[2] - (shift[2] / lapse))) + p_star - ((cov_b[2] * mag[2]) / W)),
          (lapse * sqrt(spatial_det)) * ((Etot * (vel[2] - (shift[2] / lapse))) + (p_star * vel[2]) - ((lapse * b0 * mag[2]) / W)),
          (lapse * sqrt(spatial_det)) * (((vel[2] - (shift[2] / lapse)) * mag[0]) - ((vel[0] - (shift[0] / lapse)) * mag[2])),
          (lapse * sqrt(spatial_det)) * (((vel[2] - (shift[2] / lapse)) * mag[1]) - ((vel[1] - (shift[1] / lapse)) * mag[2])),
          (lapse * sqrt(spatial_det)) * (((vel[2] - (shift[2] / lapse)) * mag[2]) - ((vel[2] - (shift[2] / lapse)) * mag[2])) },
      };

      double norm[3][3] = {
        { 1.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 },
      };

      double tau1[3][3] = {
        { 0.0, 1.0, 0.0 },
        { 1.0, 0.0, 0.0 },
        { 1.0, 0.0, 0.0 },
      };

      double tau2[3][3] = {
        { 0.0, 0.0, 1.0 },
        { 0.0, 0.0, -1.0 },
        { 0.0, 1.0, 0.0 },
      };

      double q_local[74], flux_local[74], flux[74];
      for (int d = 0; d < 3; d++) {
        gr_mhd->rotate_to_local_func(gr_mhd, tau1[d], tau2[d], norm[d], q, q_local);
        gkyl_gr_mhd_flux(gas_gamma, q_local, flux_local);
        gr_mhd->rotate_to_global_func(gr_mhd, tau1[d], tau2[d], norm[d], flux_local, flux);

        for (int i = 0; i < 8; i++) {
          TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-1) );
        }
      }

      double q_l[74], q_g[74];
      for (int d = 0; d < 3; d++) {
        gkyl_wv_eqn_rotate_to_local(gr_mhd, tau1[d], tau2[d], norm[d], q, q_l);
        gkyl_wv_eqn_rotate_to_global(gr_mhd, tau1[d], tau2[d], norm[d], q_l, q_g);

        for (int i = 0; i < 8; i++) {
          TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
        }

        double w1[74], q1[74];
        gr_mhd->cons_to_riem(gr_mhd, q_local, q_local, w1);
        gr_mhd->riem_to_cons(gr_mhd, q_local, w1, q1);

        for (int i = 0; i < 8; i++) {
          TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(extrinsic_curvature[i]);
        gkyl_free(shift_der[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der[i][j]);
        }
        gkyl_free(spatial_metric_der[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(extrinsic_curvature);
      gkyl_free(shift);
      gkyl_free(vel);
      gkyl_free(mag);
      gkyl_free(cov_mag);
      gkyl_free(cov_vel);
      gkyl_free(spacetime_vel);
      gkyl_free(b);
      gkyl_free(cov_b);
      gkyl_free(lapse_der);
      gkyl_free(shift_der);
      gkyl_free(spatial_metric_der);
    }
  }

  gkyl_wv_eqn_release(gr_mhd);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_mhd_basic_schwarzschild()
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_mhd = gkyl_wv_gr_mhd_new(gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  TEST_CHECK( gr_mhd->num_equations == 74 );
  // TEST_CHECK( gr_mhd->num_waves == 2 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, p = 1.5;
      double mag_x = 0.3, mag_y = 0.2, mag_z = 0.1;

      double spatial_det, lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der = gkyl_malloc(sizeof(double[3]));
      double **shift_der = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_det);
      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);
      
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

      spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
      spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

      double *vel = gkyl_malloc(sizeof(double[3]));
      vel[0] = u; vel[1] = v; vel[2] = w;

      double v_sq = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq += spatial_metric[i][j] * vel[i] * vel[j];
        }
      }

      double W = 1.0 / sqrt(1.0 - v_sq);

      double *mag = gkyl_malloc(sizeof(double[3]));
      mag[0] = mag_x; mag[1] = mag_y; mag[2] = mag_z;

      double *cov_mag = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        cov_mag[i] = 0.0;

        for (int j = 0; j < 3; j++) {
          cov_mag[i] += spatial_metric[i][j] * mag[j];
        }
      }

      double *cov_vel = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        cov_vel[i] = 0.0;

        for (int j = 0; j < 3; j++) {
          cov_vel[i] += spatial_metric[i][j] * vel[j];
        }
      }
      
      double b0 = 0.0;
      for (int i = 0; i < 3; i++) {
        b0 += W * mag[i] * (cov_vel[i] / lapse);
      }

      double *spacetime_vel = gkyl_malloc(sizeof(double[4]));
      spacetime_vel[0] = W / lapse;
      for (int i = 0; i < 3; i++) {
        spacetime_vel[i + 1] = (W * vel[i]) - (shift[i] * (W / lapse));
      }

      double *b = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        b[i] = (mag[i] + (lapse * b0 * spacetime_vel[i + 1])) / W;
      }

      double b_sq = 0.0;
      for (int i = 0; i < 3; i++) {
        b_sq += (mag[i] * cov_mag[i]) / (W * W);
      }
      b_sq += ((lapse * lapse) * (b0 * b0)) / (W * W);

      double *cov_b = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        cov_b[i] = 0.0;

        for (int j = 0; j < 3; j++) {
          cov_b[i] += spatial_metric[i][j] * b[j];
        }
      }

      double h_star = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0))) + (b_sq / rho);
      double p_star = p + (0.5 * b_sq);

      if (!in_excision_region) {
        double q[74];
        q[0] = sqrt(spatial_det) * rho * W;
        q[1] = sqrt(spatial_det) * ((rho * h_star * (W * W) * cov_vel[0]) - (lapse * b0 * cov_b[0]));
        q[2] = sqrt(spatial_det) * ((rho * h_star * (W * W) * cov_vel[1]) - (lapse * b0 * cov_b[1]));
        q[3] = sqrt(spatial_det) * ((rho * h_star * (W * W) * cov_vel[2]) - (lapse * b0 * cov_b[2]));
        q[4] = sqrt(spatial_det) * ((rho * h_star * (W * W)) - p_star - ((lapse * lapse) * (b0 * b0)) - (rho * W));

        q[5] = sqrt(spatial_det) * mag_x;
        q[6] = sqrt(spatial_det) * mag_y;
        q[7] = sqrt(spatial_det) * mag_z;

        q[8] = lapse;
        q[9] = shift[0]; q[10] = shift[1]; q[11] = shift[2];

        q[12] = spatial_metric[0][0]; q[13] = spatial_metric[0][1]; q[14] = spatial_metric[0][2];
        q[15] = spatial_metric[1][0]; q[16] = spatial_metric[1][1]; q[17] = spatial_metric[1][2];
        q[18] = spatial_metric[2][0]; q[19] = spatial_metric[2][1]; q[20] = spatial_metric[2][2];

        q[21] = extrinsic_curvature[0][0]; q[22] = extrinsic_curvature[0][1]; q[23] = extrinsic_curvature[0][2];
        q[24] = extrinsic_curvature[1][0]; q[25] = extrinsic_curvature[1][1]; q[26] = extrinsic_curvature[1][2];
        q[27] = extrinsic_curvature[2][0]; q[28] = extrinsic_curvature[2][1]; q[29] = extrinsic_curvature[2][2];

        q[30] = 1.0;

        q[31] = lapse_der[0]; q[32] = lapse_der[1]; q[33] = lapse_der[2];
        q[34] = shift_der[0][0]; q[35] = shift_der[0][1]; q[36] = shift_der[0][2];
        q[37] = shift_der[1][0]; q[38] = shift_der[1][1]; q[39] = shift_der[1][2];
        q[40] = shift_der[2][0]; q[41] = shift_der[2][1]; q[42] = shift_der[2][2];

        q[43] = spatial_metric_der[0][0][0]; q[44] = spatial_metric_der[0][0][1]; q[45] = spatial_metric_der[0][0][2];
        q[46] = spatial_metric_der[0][1][0]; q[47] = spatial_metric_der[0][1][1]; q[48] = spatial_metric_der[0][1][2];
        q[49] = spatial_metric_der[0][2][0]; q[50] = spatial_metric_der[0][2][1]; q[51] = spatial_metric_der[0][2][2];

        q[52] = spatial_metric_der[1][0][0]; q[53] = spatial_metric_der[1][0][1]; q[54] = spatial_metric_der[1][0][2];
        q[55] = spatial_metric_der[1][1][0]; q[56] = spatial_metric_der[1][1][1]; q[57] = spatial_metric_der[1][1][2];
        q[58] = spatial_metric_der[1][2][0]; q[59] = spatial_metric_der[1][2][1]; q[60] = spatial_metric_der[1][2][2];

        q[61] = spatial_metric_der[2][0][0]; q[62] = spatial_metric_der[2][0][1]; q[63] = spatial_metric_der[2][0][2];
        q[64] = spatial_metric_der[2][1][0]; q[65] = spatial_metric_der[2][1][1]; q[66] = spatial_metric_der[2][1][2];
        q[67] = spatial_metric_der[2][2][0]; q[68] = spatial_metric_der[2][2][1]; q[69] = spatial_metric_der[2][2][2];

        q[70] = 0.0;
        q[71] = x; q[72] = y; q[73] = 0.0;

        double prims[74];
        gkyl_gr_mhd_prim_vars(gas_gamma, q, prims);
        
        TEST_CHECK( gkyl_compare(prims[0], rho, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[1], u, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[2], v, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[3], w, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[4], p, 1e-1) );

        TEST_CHECK( gkyl_compare(prims[5], mag_x, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[6], mag_y, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[7], mag_z, 1e-1) );

        double D = rho * W;
        double Sx = (rho * h_star * (W * W) * cov_vel[0]) - (lapse * b0 * cov_b[0]);
        double Sy = (rho * h_star * (W * W) * cov_vel[1]) - (lapse * b0 * cov_b[1]);
        double Sz = (rho * h_star * (W * W) * cov_vel[2]) - (lapse * b0 * cov_b[2]);
        double Etot = (rho * h_star * (W * W)) - p_star - ((lapse * lapse) * (b0 * b0)) - (rho * W);

        double fluxes[3][8] = {
          { (lapse * sqrt(spatial_det)) * (D * (vel[0] - (shift[0] / lapse))),
            (lapse * sqrt(spatial_det)) * ((Sx * (vel[0] - (shift[0] / lapse))) + p_star - ((cov_b[0] * mag[0]) / W)),
            (lapse * sqrt(spatial_det)) * ((Sy * (vel[0] - (shift[0] / lapse))) - ((cov_b[1] * mag[0]) / W)),
            (lapse * sqrt(spatial_det)) * ((Sz * (vel[0] - (shift[0] / lapse))) - ((cov_b[2] * mag[0]) / W)),
            (lapse * sqrt(spatial_det)) * ((Etot * (vel[0] - (shift[0] / lapse))) + (p_star * vel[0]) - ((lapse * b0 * mag[0]) / W)),
            (lapse * sqrt(spatial_det)) * (((vel[0] - (shift[0] / lapse)) * mag[0]) - ((vel[0] - (shift[0] / lapse)) * mag[0])),
            (lapse * sqrt(spatial_det)) * (((vel[0] - (shift[0] / lapse)) * mag[1]) - ((vel[1] - (shift[1] / lapse)) * mag[0])),
            (lapse * sqrt(spatial_det)) * (((vel[0] - (shift[0] / lapse)) * mag[2]) - ((vel[2] - (shift[2] / lapse)) * mag[0])) },
          { (lapse * sqrt(spatial_det)) * (D * (vel[1] - (shift[1] / lapse))),
            (lapse * sqrt(spatial_det)) * ((Sx * (vel[1] - (shift[1] / lapse))) - ((cov_b[0] * mag[1]) / W)),
            (lapse * sqrt(spatial_det)) * ((Sy * (vel[1] - (shift[1] / lapse))) + p_star - ((cov_b[1] * mag[1]) / W)),
            (lapse * sqrt(spatial_det)) * ((Sz * (vel[1] - (shift[1] / lapse))) - ((cov_b[2] * mag[1]) / W)),
            (lapse * sqrt(spatial_det)) * ((Etot * (vel[1] - (shift[1] / lapse))) + (p_star * vel[1]) - ((lapse * b0 * mag[1]) / W)),
            (lapse * sqrt(spatial_det)) * (((vel[1] - (shift[1] / lapse)) * mag[0]) - ((vel[0] - (shift[0] / lapse)) * mag[1])),
            (lapse * sqrt(spatial_det)) * (((vel[1] - (shift[1] / lapse)) * mag[1]) - ((vel[1] - (shift[1] / lapse)) * mag[1])),
            (lapse * sqrt(spatial_det)) * (((vel[1] - (shift[1] / lapse)) * mag[2]) - ((vel[2] - (shift[2] / lapse)) * mag[1])) },
          { (lapse * sqrt(spatial_det)) * (D * (vel[2] - (shift[2] / lapse))),
            (lapse * sqrt(spatial_det)) * ((Sx * (vel[2] - (shift[2] / lapse))) - ((cov_b[0] * mag[2]) / W)),
            (lapse * sqrt(spatial_det)) * ((Sy * (vel[2] - (shift[2] / lapse))) - ((cov_b[1] * mag[2]) / W)),
            (lapse * sqrt(spatial_det)) * ((Sz * (vel[2] - (shift[2] / lapse))) + p_star - ((cov_b[2] * mag[2]) / W)),
            (lapse * sqrt(spatial_det)) * ((Etot * (vel[2] - (shift[2] / lapse))) + (p_star * vel[2]) - ((lapse * b0 * mag[2]) / W)),
            (lapse * sqrt(spatial_det)) * (((vel[2] - (shift[2] / lapse)) * mag[0]) - ((vel[0] - (shift[0] / lapse)) * mag[2])),
            (lapse * sqrt(spatial_det)) * (((vel[2] - (shift[2] / lapse)) * mag[1]) - ((vel[1] - (shift[1] / lapse)) * mag[2])),
            (lapse * sqrt(spatial_det)) * (((vel[2] - (shift[2] / lapse)) * mag[2]) - ((vel[2] - (shift[2] / lapse)) * mag[2])) },
        };

        double norm[3][3] = {
          { 1.0, 0.0, 0.0 },
          { 0.0, 1.0, 0.0 },
          { 0.0, 0.0, 1.0 },
        };

        double tau1[3][3] = {
          { 0.0, 1.0, 0.0 },
          { 1.0, 0.0, 0.0 },
          { 1.0, 0.0, 0.0 },
        };

        double tau2[3][3] = {
          { 0.0, 0.0, 1.0 },
          { 0.0, 0.0, -1.0 },
          { 0.0, 1.0, 0.0 },
        };

        double q_local[74], flux_local[74], flux[74];
        for (int d = 0; d < 3; d++) {
          gr_mhd->rotate_to_local_func(gr_mhd, tau1[d], tau2[d], norm[d], q, q_local);
          gkyl_gr_mhd_flux(gas_gamma, q_local, flux_local);
          gr_mhd->rotate_to_global_func(gr_mhd, tau1[d], tau2[d], norm[d], flux_local, flux);

          for (int i = 0; i < 8; i++) {
            TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-1) );
          }
        }

        double q_l[74], q_g[74];
        for (int d = 0; d < 3; d++) {
          gkyl_wv_eqn_rotate_to_local(gr_mhd, tau1[d], tau2[d], norm[d], q, q_l);
          gkyl_wv_eqn_rotate_to_global(gr_mhd, tau1[d], tau2[d], norm[d], q_l, q_g);

          for (int i = 0; i < 8; i++) {
            TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
          }

          double w1[74], q1[74];
          gr_mhd->cons_to_riem(gr_mhd, q_local, q_local, w1);
          gr_mhd->riem_to_cons(gr_mhd, q_local, w1, q1);

          for (int i = 0; i < 8; i++) {
            TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
          }
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(extrinsic_curvature[i]);
        gkyl_free(shift_der[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der[i][j]);
        }
        gkyl_free(spatial_metric_der[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(extrinsic_curvature);
      gkyl_free(shift);
      gkyl_free(vel);
      gkyl_free(mag);
      gkyl_free(cov_mag);
      gkyl_free(cov_vel);
      gkyl_free(spacetime_vel);
      gkyl_free(b);
      gkyl_free(cov_b);
      gkyl_free(lapse_der);
      gkyl_free(shift_der);
      gkyl_free(spatial_metric_der);
    }
  }

  gkyl_wv_eqn_release(gr_mhd);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_mhd_basic_kerr()
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.8, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_mhd = gkyl_wv_gr_mhd_new(gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  TEST_CHECK( gr_mhd->num_equations == 74 );
  // TEST_CHECK( gr_mhd->num_waves == 2 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, p = 1.5;
      double mag_x = 0.3, mag_y = 0.2, mag_z = 0.1;

      double spatial_det, lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der = gkyl_malloc(sizeof(double[3]));
      double **shift_der = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_det);
      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);
      
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

      spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
      spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

      double *vel = gkyl_malloc(sizeof(double[3]));
      vel[0] = u; vel[1] = v; vel[2] = w;

      double v_sq = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq += spatial_metric[i][j] * vel[i] * vel[j];
        }
      }

      double W = 1.0 / sqrt(1.0 - v_sq);

      double *mag = gkyl_malloc(sizeof(double[3]));
      mag[0] = mag_x; mag[1] = mag_y; mag[2] = mag_z;

      double *cov_mag = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        cov_mag[i] = 0.0;

        for (int j = 0; j < 3; j++) {
          cov_mag[i] += spatial_metric[i][j] * mag[j];
        }
      }

      double *cov_vel = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        cov_vel[i] = 0.0;

        for (int j = 0; j < 3; j++) {
          cov_vel[i] += spatial_metric[i][j] * vel[j];
        }
      }
      
      double b0 = 0.0;
      for (int i = 0; i < 3; i++) {
        b0 += W * mag[i] * (cov_vel[i] / lapse);
      }

      double *spacetime_vel = gkyl_malloc(sizeof(double[4]));
      spacetime_vel[0] = W / lapse;
      for (int i = 0; i < 3; i++) {
        spacetime_vel[i + 1] = (W * vel[i]) - (shift[i] * (W / lapse));
      }

      double *b = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        b[i] = (mag[i] + (lapse * b0 * spacetime_vel[i + 1])) / W;
      }

      double b_sq = 0.0;
      for (int i = 0; i < 3; i++) {
        b_sq += (mag[i] * cov_mag[i]) / (W * W);
      }
      b_sq += ((lapse * lapse) * (b0 * b0)) / (W * W);

      double *cov_b = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        cov_b[i] = 0.0;

        for (int j = 0; j < 3; j++) {
          cov_b[i] += spatial_metric[i][j] * b[j];
        }
      }

      double h_star = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0))) + (b_sq / rho);
      double p_star = p + (0.5 * b_sq);

      if (!in_excision_region) {
        double q[74];
        q[0] = sqrt(spatial_det) * rho * W;
        q[1] = sqrt(spatial_det) * ((rho * h_star * (W * W) * cov_vel[0]) - (lapse * b0 * cov_b[0]));
        q[2] = sqrt(spatial_det) * ((rho * h_star * (W * W) * cov_vel[1]) - (lapse * b0 * cov_b[1]));
        q[3] = sqrt(spatial_det) * ((rho * h_star * (W * W) * cov_vel[2]) - (lapse * b0 * cov_b[2]));
        q[4] = sqrt(spatial_det) * ((rho * h_star * (W * W)) - p_star - ((lapse * lapse) * (b0 * b0)) - (rho * W));

        q[5] = sqrt(spatial_det) * mag_x;
        q[6] = sqrt(spatial_det) * mag_y;
        q[7] = sqrt(spatial_det) * mag_z;

        q[8] = lapse;
        q[9] = shift[0]; q[10] = shift[1]; q[11] = shift[2];

        q[12] = spatial_metric[0][0]; q[13] = spatial_metric[0][1]; q[14] = spatial_metric[0][2];
        q[15] = spatial_metric[1][0]; q[16] = spatial_metric[1][1]; q[17] = spatial_metric[1][2];
        q[18] = spatial_metric[2][0]; q[19] = spatial_metric[2][1]; q[20] = spatial_metric[2][2];

        q[21] = extrinsic_curvature[0][0]; q[22] = extrinsic_curvature[0][1]; q[23] = extrinsic_curvature[0][2];
        q[24] = extrinsic_curvature[1][0]; q[25] = extrinsic_curvature[1][1]; q[26] = extrinsic_curvature[1][2];
        q[27] = extrinsic_curvature[2][0]; q[28] = extrinsic_curvature[2][1]; q[29] = extrinsic_curvature[2][2];

        q[30] = 1.0;

        q[31] = lapse_der[0]; q[32] = lapse_der[1]; q[33] = lapse_der[2];
        q[34] = shift_der[0][0]; q[35] = shift_der[0][1]; q[36] = shift_der[0][2];
        q[37] = shift_der[1][0]; q[38] = shift_der[1][1]; q[39] = shift_der[1][2];
        q[40] = shift_der[2][0]; q[41] = shift_der[2][1]; q[42] = shift_der[2][2];

        q[43] = spatial_metric_der[0][0][0]; q[44] = spatial_metric_der[0][0][1]; q[45] = spatial_metric_der[0][0][2];
        q[46] = spatial_metric_der[0][1][0]; q[47] = spatial_metric_der[0][1][1]; q[48] = spatial_metric_der[0][1][2];
        q[49] = spatial_metric_der[0][2][0]; q[50] = spatial_metric_der[0][2][1]; q[51] = spatial_metric_der[0][2][2];

        q[52] = spatial_metric_der[1][0][0]; q[53] = spatial_metric_der[1][0][1]; q[54] = spatial_metric_der[1][0][2];
        q[55] = spatial_metric_der[1][1][0]; q[56] = spatial_metric_der[1][1][1]; q[57] = spatial_metric_der[1][1][2];
        q[58] = spatial_metric_der[1][2][0]; q[59] = spatial_metric_der[1][2][1]; q[60] = spatial_metric_der[1][2][2];

        q[61] = spatial_metric_der[2][0][0]; q[62] = spatial_metric_der[2][0][1]; q[63] = spatial_metric_der[2][0][2];
        q[64] = spatial_metric_der[2][1][0]; q[65] = spatial_metric_der[2][1][1]; q[66] = spatial_metric_der[2][1][2];
        q[67] = spatial_metric_der[2][2][0]; q[68] = spatial_metric_der[2][2][1]; q[69] = spatial_metric_der[2][2][2];

        q[70] = 0.0;
        q[71] = x; q[72] = y; q[73] = 0.0;

        double prims[74];
        gkyl_gr_mhd_prim_vars(gas_gamma, q, prims);
        
        TEST_CHECK( gkyl_compare(prims[0], rho, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[1], u, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[2], v, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[3], w, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[4], p, 1e-1) );

        TEST_CHECK( gkyl_compare(prims[5], mag_x, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[6], mag_y, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[7], mag_z, 1e-1) );

        double D = rho * W;
        double Sx = (rho * h_star * (W * W) * cov_vel[0]) - (lapse * b0 * cov_b[0]);
        double Sy = (rho * h_star * (W * W) * cov_vel[1]) - (lapse * b0 * cov_b[1]);
        double Sz = (rho * h_star * (W * W) * cov_vel[2]) - (lapse * b0 * cov_b[2]);
        double Etot = (rho * h_star * (W * W)) - p_star - ((lapse * lapse) * (b0 * b0)) - (rho * W);

        double fluxes[3][8] = {
          { (lapse * sqrt(spatial_det)) * (D * (vel[0] - (shift[0] / lapse))),
            (lapse * sqrt(spatial_det)) * ((Sx * (vel[0] - (shift[0] / lapse))) + p_star - ((cov_b[0] * mag[0]) / W)),
            (lapse * sqrt(spatial_det)) * ((Sy * (vel[0] - (shift[0] / lapse))) - ((cov_b[1] * mag[0]) / W)),
            (lapse * sqrt(spatial_det)) * ((Sz * (vel[0] - (shift[0] / lapse))) - ((cov_b[2] * mag[0]) / W)),
            (lapse * sqrt(spatial_det)) * ((Etot * (vel[0] - (shift[0] / lapse))) + (p_star * vel[0]) - ((lapse * b0 * mag[0]) / W)),
            (lapse * sqrt(spatial_det)) * (((vel[0] - (shift[0] / lapse)) * mag[0]) - ((vel[0] - (shift[0] / lapse)) * mag[0])),
            (lapse * sqrt(spatial_det)) * (((vel[0] - (shift[0] / lapse)) * mag[1]) - ((vel[1] - (shift[1] / lapse)) * mag[0])),
            (lapse * sqrt(spatial_det)) * (((vel[0] - (shift[0] / lapse)) * mag[2]) - ((vel[2] - (shift[2] / lapse)) * mag[0])) },
          { (lapse * sqrt(spatial_det)) * (D * (vel[1] - (shift[1] / lapse))),
            (lapse * sqrt(spatial_det)) * ((Sx * (vel[1] - (shift[1] / lapse))) - ((cov_b[0] * mag[1]) / W)),
            (lapse * sqrt(spatial_det)) * ((Sy * (vel[1] - (shift[1] / lapse))) + p_star - ((cov_b[1] * mag[1]) / W)),
            (lapse * sqrt(spatial_det)) * ((Sz * (vel[1] - (shift[1] / lapse))) - ((cov_b[2] * mag[1]) / W)),
            (lapse * sqrt(spatial_det)) * ((Etot * (vel[1] - (shift[1] / lapse))) + (p_star * vel[1]) - ((lapse * b0 * mag[1]) / W)),
            (lapse * sqrt(spatial_det)) * (((vel[1] - (shift[1] / lapse)) * mag[0]) - ((vel[0] - (shift[0] / lapse)) * mag[1])),
            (lapse * sqrt(spatial_det)) * (((vel[1] - (shift[1] / lapse)) * mag[1]) - ((vel[1] - (shift[1] / lapse)) * mag[1])),
            (lapse * sqrt(spatial_det)) * (((vel[1] - (shift[1] / lapse)) * mag[2]) - ((vel[2] - (shift[2] / lapse)) * mag[1])) },
          { (lapse * sqrt(spatial_det)) * (D * (vel[2] - (shift[2] / lapse))),
            (lapse * sqrt(spatial_det)) * ((Sx * (vel[2] - (shift[2] / lapse))) - ((cov_b[0] * mag[2]) / W)),
            (lapse * sqrt(spatial_det)) * ((Sy * (vel[2] - (shift[2] / lapse))) - ((cov_b[1] * mag[2]) / W)),
            (lapse * sqrt(spatial_det)) * ((Sz * (vel[2] - (shift[2] / lapse))) + p_star - ((cov_b[2] * mag[2]) / W)),
            (lapse * sqrt(spatial_det)) * ((Etot * (vel[2] - (shift[2] / lapse))) + (p_star * vel[2]) - ((lapse * b0 * mag[2]) / W)),
            (lapse * sqrt(spatial_det)) * (((vel[2] - (shift[2] / lapse)) * mag[0]) - ((vel[0] - (shift[0] / lapse)) * mag[2])),
            (lapse * sqrt(spatial_det)) * (((vel[2] - (shift[2] / lapse)) * mag[1]) - ((vel[1] - (shift[1] / lapse)) * mag[2])),
            (lapse * sqrt(spatial_det)) * (((vel[2] - (shift[2] / lapse)) * mag[2]) - ((vel[2] - (shift[2] / lapse)) * mag[2])) },
        };

        double norm[3][3] = {
          { 1.0, 0.0, 0.0 },
          { 0.0, 1.0, 0.0 },
          { 0.0, 0.0, 1.0 },
        };

        double tau1[3][3] = {
          { 0.0, 1.0, 0.0 },
          { 1.0, 0.0, 0.0 },
          { 1.0, 0.0, 0.0 },
        };

        double tau2[3][3] = {
          { 0.0, 0.0, 1.0 },
          { 0.0, 0.0, -1.0 },
          { 0.0, 1.0, 0.0 },
        };

        double q_local[74], flux_local[74], flux[74];
        for (int d = 0; d < 3; d++) {
          gr_mhd->rotate_to_local_func(gr_mhd, tau1[d], tau2[d], norm[d], q, q_local);
          gkyl_gr_mhd_flux(gas_gamma, q_local, flux_local);
          gr_mhd->rotate_to_global_func(gr_mhd, tau1[d], tau2[d], norm[d], flux_local, flux);

          for (int i = 0; i < 8; i++) {
            TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-1) );
          }
        }

        double q_l[74], q_g[74];
        for (int d = 0; d < 3; d++) {
          gkyl_wv_eqn_rotate_to_local(gr_mhd, tau1[d], tau2[d], norm[d], q, q_l);
          gkyl_wv_eqn_rotate_to_global(gr_mhd, tau1[d], tau2[d], norm[d], q_l, q_g);

          for (int i = 0; i < 8; i++) {
            TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
          }

          double w1[74], q1[74];
          gr_mhd->cons_to_riem(gr_mhd, q_local, q_local, w1);
          gr_mhd->riem_to_cons(gr_mhd, q_local, w1, q1);

          for (int i = 0; i < 8; i++) {
            TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
          }
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(extrinsic_curvature[i]);
        gkyl_free(shift_der[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der[i][j]);
        }
        gkyl_free(spatial_metric_der[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(extrinsic_curvature);
      gkyl_free(shift);
      gkyl_free(vel);
      gkyl_free(mag);
      gkyl_free(cov_mag);
      gkyl_free(cov_vel);
      gkyl_free(spacetime_vel);
      gkyl_free(b);
      gkyl_free(cov_b);
      gkyl_free(lapse_der);
      gkyl_free(shift_der);
      gkyl_free(spatial_metric_der);
    }
  }

  gkyl_wv_eqn_release(gr_mhd);
  gkyl_gr_spacetime_release(spacetime);
}

TEST_LIST = {
  { "gr_mhd_basic_minkowski" , test_gr_mhd_basic_minkowski},
  { "gr_mhd_basic_schwarzschild", test_gr_mhd_basic_schwarzschild },
  { "gr_mhd_basic_kerr", test_gr_mhd_basic_kerr },
  { NULL, NULL },
};