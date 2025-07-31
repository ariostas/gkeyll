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
  TEST_CHECK( gr_mhd->num_waves == 2 );

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
      
      TEST_CHECK( gkyl_compare(prims[0], rho, 1e-8) );
      TEST_CHECK( gkyl_compare(prims[1], u, 1e-8) );
      TEST_CHECK( gkyl_compare(prims[2], v, 1e-8) );
      TEST_CHECK( gkyl_compare(prims[3], w, 1e-8) );
      TEST_CHECK( gkyl_compare(prims[4], p, 1e-8) );

      TEST_CHECK( gkyl_compare(prims[5], mag_x, 1e-8) );
      TEST_CHECK( gkyl_compare(prims[6], mag_y, 1e-8) );
      TEST_CHECK( gkyl_compare(prims[7], mag_z, 1e-8) );

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
          TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-8) );
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
  TEST_CHECK( gr_mhd->num_waves == 2 );

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

        TEST_CHECK( gkyl_compare(prims[5], mag_x, 1e-8) );
        TEST_CHECK( gkyl_compare(prims[6], mag_y, 1e-8) );
        TEST_CHECK( gkyl_compare(prims[7], mag_z, 1e-8) );

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
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_mhd = gkyl_wv_gr_mhd_new(gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  TEST_CHECK( gr_mhd->num_equations == 74 );
  TEST_CHECK( gr_mhd->num_waves == 2 );

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

        TEST_CHECK( gkyl_compare(prims[5], mag_x, 1e-8) );
        TEST_CHECK( gkyl_compare(prims[6], mag_y, 1e-8) );
        TEST_CHECK( gkyl_compare(prims[7], mag_z, 1e-8) );

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
test_gr_mhd_waves_minkowski()
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  struct gkyl_wv_eqn *gr_mhd = gkyl_wv_gr_mhd_new(gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho_l = 1.0, u_l = 0.1, v_l = 0.2, w_l = 0.3, p_l = 1.5;
      double mag_x_l = 0.3, mag_y_l = 0.2, mag_z_l = 0.1;
      double rho_r = 0.1, u_r = 0.2, v_r = 0.3, w_r = 0.4, p_r = 0.15;
      double mag_x_r = 0.4, mag_y_r = 0.3, mag_z_r = 0.2;

      double spatial_det_l, spatial_det_r;
      double lapse_l, lapse_r;
      double *shift_l = gkyl_malloc(sizeof(double[3]));
      double *shift_r = gkyl_malloc(sizeof(double[3]));

      double **spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
      double **spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature_l = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
        spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature_l[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der_l = gkyl_malloc(sizeof(double[3]));
      double *lapse_der_r = gkyl_malloc(sizeof(double[3]));
      double **shift_der_l = gkyl_malloc(sizeof(double*[3]));
      double **shift_der_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der_l[i] = gkyl_malloc(sizeof(double[3]));
        shift_der_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der_l = gkyl_malloc(sizeof(double**[3]));
      double ***spatial_metric_der_r = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der_l[i] = gkyl_malloc(sizeof(double*[3]));
        spatial_metric_der_r[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der_l[i][j] = gkyl_malloc(sizeof(double[3]));
          spatial_metric_der_r[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_det_l);
      spacetime->spatial_metric_det_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_det_r);
      spacetime->lapse_function_func(spacetime, 0.0, x - 0.1, y, 0.0, &lapse_l);
      spacetime->lapse_function_func(spacetime, 0.0, x + 0.1, y, 0.0, &lapse_r);
      spacetime->shift_vector_func(spacetime, 0.0, x - 0.1, y, 0.0, &shift_l);
      spacetime->shift_vector_func(spacetime, 0.0, x + 0.1, y, 0.0, &shift_r);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_metric_l);
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_metric_r);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_l);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_r);

      spacetime->lapse_function_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der_l);
      spacetime->lapse_function_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der_r);
      spacetime->shift_vector_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der_l);
      spacetime->shift_vector_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der_r);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der_l);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der_r);

      double *vel_l = gkyl_malloc(sizeof(double[3]));
      double *vel_r = gkyl_malloc(sizeof(double[3]));
      vel_l[0] = u_l; vel_l[1] = v_l; vel_l[2] = w_l;
      vel_r[0] = u_r; vel_r[1] = v_r; vel_r[2] = w_r;

      double v_sq_l = 0.0, v_sq_r = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq_l += spatial_metric_l[i][j] * vel_l[i] * vel_l[j];
          v_sq_r += spatial_metric_r[i][j] * vel_r[i] * vel_r[j];
        }
      }

      double W_l = 1.0 / sqrt(1.0 - v_sq_l);
      double W_r = 1.0 / sqrt(1.0 - v_sq_r);

      double *mag_l = gkyl_malloc(sizeof(double[3]));
      double *mag_r = gkyl_malloc(sizeof(double[3]));
      mag_l[0] = mag_x_l; mag_l[1] = mag_y_l; mag_l[2] = mag_z_l;
      mag_r[0] = mag_x_r; mag_r[1] = mag_y_r; mag_r[2] = mag_z_r;

      double *cov_mag_l = gkyl_malloc(sizeof(double[3]));
      double *cov_mag_r = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        cov_mag_l[i] = 0.0;
        cov_mag_r[i] = 0.0;

        for (int j = 0; j < 3; j++) {
          cov_mag_l[i] += spatial_metric_l[i][j] * mag_l[j];
          cov_mag_r[i] += spatial_metric_r[i][j] * mag_r[j];
        }
      }

      double *cov_vel_l = gkyl_malloc(sizeof(double[3]));
      double *cov_vel_r = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        cov_vel_l[i] = 0.0;
        cov_vel_r[i] = 0.0;

        for (int j = 0; j < 3; j++) {
          cov_vel_l[i] += spatial_metric_l[i][j] * vel_l[j];
          cov_vel_r[i] += spatial_metric_r[i][j] * vel_r[j];
        }
      }
      
      double b0_l = 0.0, b0_r = 0.0;
      for (int i = 0; i < 3; i++) {
        b0_l += W_l * mag_l[i] * (cov_vel_l[i] / lapse_l);
        b0_r += W_r * mag_r[i] * (cov_vel_r[i] / lapse_r);
      }

      double *spacetime_vel_l = gkyl_malloc(sizeof(double[4]));
      double *spacetime_vel_r = gkyl_malloc(sizeof(double[4]));
      spacetime_vel_l[0] = W_l / lapse_l;
      spacetime_vel_r[0] = W_r / lapse_r;
      for (int i = 0; i < 3; i++) {
        spacetime_vel_l[i + 1] = (W_l * vel_l[i]) - (shift_l[i] * (W_l / lapse_l));
        spacetime_vel_r[i + 1] = (W_r * vel_r[i]) - (shift_r[i] * (W_r / lapse_r));
      }

      double *b_l = gkyl_malloc(sizeof(double[3]));
      double *b_r = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        b_l[i] = (mag_l[i] + (lapse_l * b0_l * spacetime_vel_l[i + 1])) / W_l;
        b_r[i] = (mag_r[i] + (lapse_r * b0_r * spacetime_vel_r[i + 1])) / W_r;
      }

      double b_sq_l = 0.0, b_sq_r = 0.0;
      for (int i = 0; i < 3; i++) {
        b_sq_l += (mag_l[i] * cov_mag_l[i]) / (W_l * W_l);
        b_sq_r += (mag_r[i] * cov_mag_r[i]) / (W_r * W_r);
      }
      b_sq_l += ((lapse_l * lapse_l) * (b0_l * b0_l)) / (W_l * W_l);
      b_sq_r += ((lapse_r * lapse_r) * (b0_r * b0_r)) / (W_r * W_r);

      double *cov_b_l = gkyl_malloc(sizeof(double[3]));
      double *cov_b_r = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        cov_b_l[i] = 0.0;
        cov_b_r[i] = 0.0;

        for (int j = 0; j < 3; j++) {
          cov_b_l[i] += spatial_metric_l[i][j] * b_l[j];
          cov_b_r[i] += spatial_metric_r[i][j] * b_r[j];
        }
      }

      double h_star_l = 1.0 + ((p_l / rho_l) * (gas_gamma / (gas_gamma - 1.0))) + (b_sq_l / rho_l);
      double h_star_r = 1.0 + ((p_r / rho_r) * (gas_gamma / (gas_gamma - 1.0))) + (b_sq_r / rho_r);
      double p_star_l = p_l + (0.5 * b_sq_l);
      double p_star_r = p_r + (0.5 * b_sq_r);

      double ql[74], qr[74];
      ql[0] = sqrt(spatial_det_l) * rho_l * W_l;
      ql[1] = sqrt(spatial_det_l) * ((rho_l * h_star_l * (W_l * W_l) * cov_vel_l[0]) - (lapse_l * b0_l * cov_b_l[0]));
      ql[2] = sqrt(spatial_det_l) * ((rho_l * h_star_l * (W_l * W_l) * cov_vel_l[1]) - (lapse_l * b0_l * cov_b_l[1]));
      ql[3] = sqrt(spatial_det_l) * ((rho_l * h_star_l * (W_l * W_l) * cov_vel_l[2]) - (lapse_l * b0_l * cov_b_l[2]));
      ql[4] = sqrt(spatial_det_l) * ((rho_l * h_star_l * (W_l * W_l)) - p_star_l - ((lapse_l * lapse_l) * (b0_l * b0_l)) - (rho_l * W_l));

      ql[5] = sqrt(spatial_det_l) * mag_x_l;
      ql[6] = sqrt(spatial_det_l) * mag_y_l;
      ql[7] = sqrt(spatial_det_l) * mag_z_l;

      ql[8] = lapse_l;
      ql[9] = shift_l[0]; ql[10] = shift_l[1]; ql[11] = shift_l[2];

      ql[12] = spatial_metric_l[0][0]; ql[13] = spatial_metric_l[0][1]; ql[14] = spatial_metric_l[0][2];
      ql[15] = spatial_metric_l[1][0]; ql[16] = spatial_metric_l[1][1]; ql[17] = spatial_metric_l[1][2];
      ql[18] = spatial_metric_l[2][0]; ql[19] = spatial_metric_l[2][1]; ql[20] = spatial_metric_l[2][2];

      ql[21] = extrinsic_curvature_l[0][0]; ql[22] = extrinsic_curvature_l[0][1]; ql[23] = extrinsic_curvature_l[0][2];
      ql[24] = extrinsic_curvature_l[1][0]; ql[25] = extrinsic_curvature_l[1][1]; ql[26] = extrinsic_curvature_l[1][2];
      ql[27] = extrinsic_curvature_l[2][0]; ql[28] = extrinsic_curvature_l[2][1]; ql[29] = extrinsic_curvature_l[2][2];

      ql[30] = 1.0;

      ql[31] = lapse_der_l[0]; ql[32] = lapse_der_l[1]; ql[33] = lapse_der_l[2];
      ql[34] = shift_der_l[0][0]; ql[35] = shift_der_l[0][1]; ql[36] = shift_der_l[0][2];
      ql[37] = shift_der_l[1][0]; ql[38] = shift_der_l[1][1]; ql[39] = shift_der_l[1][2];
      ql[40] = shift_der_l[2][0]; ql[41] = shift_der_l[2][1]; ql[42] = shift_der_l[2][2];

      ql[43] = spatial_metric_der_l[0][0][0]; ql[44] = spatial_metric_der_l[0][0][1]; ql[45] = spatial_metric_der_l[0][0][2];
      ql[46] = spatial_metric_der_l[0][1][0]; ql[47] = spatial_metric_der_l[0][1][1]; ql[48] = spatial_metric_der_l[0][1][2];
      ql[49] = spatial_metric_der_l[0][2][0]; ql[50] = spatial_metric_der_l[0][2][1]; ql[51] = spatial_metric_der_l[0][2][2];

      ql[52] = spatial_metric_der_l[1][0][0]; ql[53] = spatial_metric_der_l[1][0][1]; ql[54] = spatial_metric_der_l[1][0][2];
      ql[55] = spatial_metric_der_l[1][1][0]; ql[56] = spatial_metric_der_l[1][1][1]; ql[57] = spatial_metric_der_l[1][1][2];
      ql[58] = spatial_metric_der_l[1][2][0]; ql[59] = spatial_metric_der_l[1][2][1]; ql[60] = spatial_metric_der_l[1][2][2];

      ql[61] = spatial_metric_der_l[2][0][0]; ql[62] = spatial_metric_der_l[2][0][1]; ql[63] = spatial_metric_der_l[2][0][2];
      ql[64] = spatial_metric_der_l[2][1][0]; ql[65] = spatial_metric_der_l[2][1][1]; ql[66] = spatial_metric_der_l[2][1][2];
      ql[67] = spatial_metric_der_l[2][2][0]; ql[68] = spatial_metric_der_l[2][2][1]; ql[69] = spatial_metric_der_l[2][2][2];

      ql[70] = 0.0;
      ql[71] = x - 0.5; ql[72] = y; ql[73] = 0.0;

      qr[0] = sqrt(spatial_det_r) * rho_r * W_r;
      qr[1] = sqrt(spatial_det_r) * ((rho_r * h_star_r * (W_r * W_r) * cov_vel_r[0]) - (lapse_r * b0_r * cov_b_r[0]));
      qr[2] = sqrt(spatial_det_r) * ((rho_r * h_star_r * (W_r * W_r) * cov_vel_r[1]) - (lapse_r * b0_r * cov_b_r[1]));
      qr[3] = sqrt(spatial_det_r) * ((rho_r * h_star_r * (W_r * W_r) * cov_vel_r[2]) - (lapse_r * b0_r * cov_b_r[2]));
      qr[4] = sqrt(spatial_det_r) * ((rho_r * h_star_r * (W_r * W_r)) - p_star_r - ((lapse_r * lapse_r) * (b0_r * b0_r)) - (rho_r * W_r));

      qr[5] = sqrt(spatial_det_r) * mag_x_r;
      qr[6] = sqrt(spatial_det_r) * mag_y_r;
      qr[7] = sqrt(spatial_det_r) * mag_z_r;

      qr[8] = lapse_r;
      qr[9] = shift_r[0]; qr[10] = shift_r[1]; qr[11] = shift_r[2];

      qr[12] = spatial_metric_r[0][0]; qr[13] = spatial_metric_r[0][1]; qr[14] = spatial_metric_r[0][2];
      qr[15] = spatial_metric_r[1][0]; qr[16] = spatial_metric_r[1][1]; qr[17] = spatial_metric_r[1][2];
      qr[18] = spatial_metric_r[2][0]; qr[19] = spatial_metric_r[2][1]; qr[20] = spatial_metric_r[2][2];

      qr[21] = extrinsic_curvature_r[0][0]; qr[22] = extrinsic_curvature_r[0][1]; qr[23] = extrinsic_curvature_r[0][2];
      qr[24] = extrinsic_curvature_r[1][0]; qr[25] = extrinsic_curvature_r[1][1]; qr[26] = extrinsic_curvature_r[1][2];
      qr[27] = extrinsic_curvature_r[2][0]; qr[28] = extrinsic_curvature_r[2][1]; qr[29] = extrinsic_curvature_r[2][2];

      qr[30] = 1.0;

      qr[31] = lapse_der_r[0]; qr[32] = lapse_der_r[1]; qr[33] = lapse_der_r[2];
      qr[34] = shift_der_r[0][0]; qr[35] = shift_der_r[0][1]; qr[36] = shift_der_r[0][2];
      qr[37] = shift_der_r[1][0]; qr[38] = shift_der_r[1][1]; qr[39] = shift_der_r[1][2];
      qr[40] = shift_der_r[2][0]; qr[41] = shift_der_r[2][1]; qr[42] = shift_der_r[2][2];

      qr[43] = spatial_metric_der_r[0][0][0]; qr[44] = spatial_metric_der_r[0][0][1]; qr[45] = spatial_metric_der_r[0][0][2];
      qr[46] = spatial_metric_der_r[0][1][0]; qr[47] = spatial_metric_der_r[0][1][1]; qr[48] = spatial_metric_der_r[0][1][2];
      qr[49] = spatial_metric_der_r[0][2][0]; qr[50] = spatial_metric_der_r[0][2][1]; qr[51] = spatial_metric_der_r[0][2][2];

      qr[52] = spatial_metric_der_r[1][0][0]; qr[53] = spatial_metric_der_r[1][0][1]; qr[54] = spatial_metric_der_r[1][0][2];
      qr[55] = spatial_metric_der_r[1][1][0]; qr[56] = spatial_metric_der_r[1][1][1]; qr[57] = spatial_metric_der_r[1][1][2];
      qr[58] = spatial_metric_der_r[1][2][0]; qr[59] = spatial_metric_der_r[1][2][1]; qr[60] = spatial_metric_der_r[1][2][2];

      qr[61] = spatial_metric_der_r[2][0][0]; qr[62] = spatial_metric_der_r[2][0][1]; qr[63] = spatial_metric_der_r[2][0][2];
      qr[64] = spatial_metric_der_r[2][1][0]; qr[65] = spatial_metric_der_r[2][1][1]; qr[66] = spatial_metric_der_r[2][1][2];
      qr[67] = spatial_metric_der_r[2][2][0]; qr[68] = spatial_metric_der_r[2][2][1]; qr[69] = spatial_metric_der_r[2][2][2];

      qr[70] = 0.0;
      qr[71] = x + 0.5; qr[72] = y; qr[73] = 0.0;

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

      for (int d = 0; d < 3; d++) {
        double speeds[2], waves[2 * 74], waves_local[2 * 74];

        double ql_local[74], qr_local[74];
        gkyl_wv_eqn_rotate_to_local(gr_mhd, tau1[d], tau2[d], norm[d], ql, ql_local);
        gkyl_wv_eqn_rotate_to_local(gr_mhd, tau1[d], tau2[d], norm[d], qr, qr_local);

        double delta[74];
        for (int i = 0; i < 74; i++) {
          delta[i] = qr_local[i] - ql_local[i];
        }

        gkyl_wv_eqn_waves(gr_mhd, GKYL_WV_LOW_ORDER_FLUX, delta, ql_local, qr_local, waves_local, speeds);

        double apdq_local[74], amdq_local[74];
        gkyl_wv_eqn_qfluct(gr_mhd, GKYL_WV_LOW_ORDER_FLUX, ql_local, qr_local, waves_local, speeds, amdq_local, apdq_local);

        for (int i = 0; i < 2; i++) {
          gkyl_wv_eqn_rotate_to_global(gr_mhd, tau1[d], tau2[d], norm[d], &waves_local[i * 74], &waves[i * 74]);
        }

        double apdq[74], amdq[74];
        gkyl_wv_eqn_rotate_to_global(gr_mhd, tau1[d], tau2[d], norm[d], apdq_local, apdq);
        gkyl_wv_eqn_rotate_to_global(gr_mhd, tau1[d], tau2[d], norm[d], amdq_local, amdq);

        double fl_local[74], fr_local[74];
        gkyl_gr_mhd_flux(gas_gamma, ql_local, fl_local);
        gkyl_gr_mhd_flux(gas_gamma, qr_local, fr_local);

        double fl[74], fr[74];
        gkyl_wv_eqn_rotate_to_global(gr_mhd, tau1[d], tau2[d], norm[d], fl_local, fl);
        gkyl_wv_eqn_rotate_to_global(gr_mhd, tau1[d], tau2[d], norm[d], fr_local, fr);

        for (int i = 0; i < 74; i++) {
          TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-15) );
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric_l[i]);
        gkyl_free(spatial_metric_r[i]);
        gkyl_free(extrinsic_curvature_l[i]);
        gkyl_free(extrinsic_curvature_r[i]);
        gkyl_free(shift_der_l[i]);
        gkyl_free(shift_der_r[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der_l[i][j]);
          gkyl_free(spatial_metric_der_r[i][j]);
        }
        gkyl_free(spatial_metric_der_l[i]);
        gkyl_free(spatial_metric_der_r[i]);
      }
      gkyl_free(spatial_metric_l);
      gkyl_free(spatial_metric_r);
      gkyl_free(extrinsic_curvature_l);
      gkyl_free(extrinsic_curvature_r);
      gkyl_free(shift_l);
      gkyl_free(shift_r);
      gkyl_free(vel_l);
      gkyl_free(vel_r);
      gkyl_free(mag_l);
      gkyl_free(mag_r);
      gkyl_free(cov_mag_l);
      gkyl_free(cov_mag_r);
      gkyl_free(cov_vel_l);
      gkyl_free(cov_vel_r);
      gkyl_free(spacetime_vel_l);
      gkyl_free(spacetime_vel_r);
      gkyl_free(b_l);
      gkyl_free(b_r);
      gkyl_free(cov_b_l);
      gkyl_free(cov_b_r);
      gkyl_free(lapse_der_l);
      gkyl_free(lapse_der_r);
      gkyl_free(shift_der_l);
      gkyl_free(shift_der_r);
      gkyl_free(spatial_metric_der_l);
      gkyl_free(spatial_metric_der_r);
    }
  }

  gkyl_wv_eqn_release(gr_mhd);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_mhd_waves_schwarzschild()
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_mhd = gkyl_wv_gr_mhd_new(gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho_l = 1.0, u_l = 0.1, v_l = 0.2, w_l = 0.3, p_l = 1.5;
      double mag_x_l = 0.3, mag_y_l = 0.2, mag_z_l = 0.1;
      double rho_r = 0.1, u_r = 0.2, v_r = 0.3, w_r = 0.4, p_r = 0.15;
      double mag_x_r = 0.4, mag_y_r = 0.3, mag_z_r = 0.2;

      double spatial_det_l, spatial_det_r;
      double lapse_l, lapse_r;
      double *shift_l = gkyl_malloc(sizeof(double[3]));
      double *shift_r = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region_l, in_excision_region_r;

      double **spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
      double **spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature_l = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
        spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature_l[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der_l = gkyl_malloc(sizeof(double[3]));
      double *lapse_der_r = gkyl_malloc(sizeof(double[3]));
      double **shift_der_l = gkyl_malloc(sizeof(double*[3]));
      double **shift_der_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der_l[i] = gkyl_malloc(sizeof(double[3]));
        shift_der_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der_l = gkyl_malloc(sizeof(double**[3]));
      double ***spatial_metric_der_r = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der_l[i] = gkyl_malloc(sizeof(double*[3]));
        spatial_metric_der_r[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der_l[i][j] = gkyl_malloc(sizeof(double[3]));
          spatial_metric_der_r[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_det_l);
      spacetime->spatial_metric_det_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_det_r);
      spacetime->lapse_function_func(spacetime, 0.0, x - 0.1, y, 0.0, &lapse_l);
      spacetime->lapse_function_func(spacetime, 0.0, x + 0.1, y, 0.0, &lapse_r);
      spacetime->shift_vector_func(spacetime, 0.0, x - 0.1, y, 0.0, &shift_l);
      spacetime->shift_vector_func(spacetime, 0.0, x + 0.1, y, 0.0, &shift_r);
      spacetime->excision_region_func(spacetime, 0.0, x - 0.1, y, 0.0, &in_excision_region_l);
      spacetime->excision_region_func(spacetime, 0.0, x + 0.1, y, 0.0, &in_excision_region_r);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_metric_l);
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_metric_r);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_l);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_r);

      spacetime->lapse_function_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der_l);
      spacetime->lapse_function_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der_r);
      spacetime->shift_vector_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der_l);
      spacetime->shift_vector_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der_r);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der_l);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der_r);

      double *vel_l = gkyl_malloc(sizeof(double[3]));
      double *vel_r = gkyl_malloc(sizeof(double[3]));
      vel_l[0] = u_l; vel_l[1] = v_l; vel_l[2] = w_l;
      vel_r[0] = u_r; vel_r[1] = v_r; vel_r[2] = w_r;

      double v_sq_l = 0.0, v_sq_r = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq_l += spatial_metric_l[i][j] * vel_l[i] * vel_l[j];
          v_sq_r += spatial_metric_r[i][j] * vel_r[i] * vel_r[j];
        }
      }

      double W_l = 1.0 / sqrt(1.0 - v_sq_l);
      double W_r = 1.0 / sqrt(1.0 - v_sq_r);

      double *mag_l = gkyl_malloc(sizeof(double[3]));
      double *mag_r = gkyl_malloc(sizeof(double[3]));
      mag_l[0] = mag_x_l; mag_l[1] = mag_y_l; mag_l[2] = mag_z_l;
      mag_r[0] = mag_x_r; mag_r[1] = mag_y_r; mag_r[2] = mag_z_r;

      double *cov_mag_l = gkyl_malloc(sizeof(double[3]));
      double *cov_mag_r = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        cov_mag_l[i] = 0.0;
        cov_mag_r[i] = 0.0;

        for (int j = 0; j < 3; j++) {
          cov_mag_l[i] += spatial_metric_l[i][j] * mag_l[j];
          cov_mag_r[i] += spatial_metric_r[i][j] * mag_r[j];
        }
      }

      double *cov_vel_l = gkyl_malloc(sizeof(double[3]));
      double *cov_vel_r = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        cov_vel_l[i] = 0.0;
        cov_vel_r[i] = 0.0;

        for (int j = 0; j < 3; j++) {
          cov_vel_l[i] += spatial_metric_l[i][j] * vel_l[j];
          cov_vel_r[i] += spatial_metric_r[i][j] * vel_r[j];
        }
      }
      
      double b0_l = 0.0, b0_r = 0.0;
      for (int i = 0; i < 3; i++) {
        b0_l += W_l * mag_l[i] * (cov_vel_l[i] / lapse_l);
        b0_r += W_r * mag_r[i] * (cov_vel_r[i] / lapse_r);
      }

      double *spacetime_vel_l = gkyl_malloc(sizeof(double[4]));
      double *spacetime_vel_r = gkyl_malloc(sizeof(double[4]));
      spacetime_vel_l[0] = W_l / lapse_l;
      spacetime_vel_r[0] = W_r / lapse_r;
      for (int i = 0; i < 3; i++) {
        spacetime_vel_l[i + 1] = (W_l * vel_l[i]) - (shift_l[i] * (W_l / lapse_l));
        spacetime_vel_r[i + 1] = (W_r * vel_r[i]) - (shift_r[i] * (W_r / lapse_r));
      }

      double *b_l = gkyl_malloc(sizeof(double[3]));
      double *b_r = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        b_l[i] = (mag_l[i] + (lapse_l * b0_l * spacetime_vel_l[i + 1])) / W_l;
        b_r[i] = (mag_r[i] + (lapse_r * b0_r * spacetime_vel_r[i + 1])) / W_r;
      }

      double b_sq_l = 0.0, b_sq_r = 0.0;
      for (int i = 0; i < 3; i++) {
        b_sq_l += (mag_l[i] * cov_mag_l[i]) / (W_l * W_l);
        b_sq_r += (mag_r[i] * cov_mag_r[i]) / (W_r * W_r);
      }
      b_sq_l += ((lapse_l * lapse_l) * (b0_l * b0_l)) / (W_l * W_l);
      b_sq_r += ((lapse_r * lapse_r) * (b0_r * b0_r)) / (W_r * W_r);

      double *cov_b_l = gkyl_malloc(sizeof(double[3]));
      double *cov_b_r = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        cov_b_l[i] = 0.0;
        cov_b_r[i] = 0.0;

        for (int j = 0; j < 3; j++) {
          cov_b_l[i] += spatial_metric_l[i][j] * b_l[j];
          cov_b_r[i] += spatial_metric_r[i][j] * b_r[j];
        }
      }

      double h_star_l = 1.0 + ((p_l / rho_l) * (gas_gamma / (gas_gamma - 1.0))) + (b_sq_l / rho_l);
      double h_star_r = 1.0 + ((p_r / rho_r) * (gas_gamma / (gas_gamma - 1.0))) + (b_sq_r / rho_r);
      double p_star_l = p_l + (0.5 * b_sq_l);
      double p_star_r = p_r + (0.5 * b_sq_r);

      if (!in_excision_region_l && !in_excision_region_r) {
        double ql[74], qr[74];
        ql[0] = sqrt(spatial_det_l) * rho_l * W_l;
        ql[1] = sqrt(spatial_det_l) * ((rho_l * h_star_l * (W_l * W_l) * cov_vel_l[0]) - (lapse_l * b0_l * cov_b_l[0]));
        ql[2] = sqrt(spatial_det_l) * ((rho_l * h_star_l * (W_l * W_l) * cov_vel_l[1]) - (lapse_l * b0_l * cov_b_l[1]));
        ql[3] = sqrt(spatial_det_l) * ((rho_l * h_star_l * (W_l * W_l) * cov_vel_l[2]) - (lapse_l * b0_l * cov_b_l[2]));
        ql[4] = sqrt(spatial_det_l) * ((rho_l * h_star_l * (W_l * W_l)) - p_star_l - ((lapse_l * lapse_l) * (b0_l * b0_l)) - (rho_l * W_l));

        ql[5] = sqrt(spatial_det_l) * mag_x_l;
        ql[6] = sqrt(spatial_det_l) * mag_y_l;
        ql[7] = sqrt(spatial_det_l) * mag_z_l;

        ql[8] = lapse_l;
        ql[9] = shift_l[0]; ql[10] = shift_l[1]; ql[11] = shift_l[2];

        ql[12] = spatial_metric_l[0][0]; ql[13] = spatial_metric_l[0][1]; ql[14] = spatial_metric_l[0][2];
        ql[15] = spatial_metric_l[1][0]; ql[16] = spatial_metric_l[1][1]; ql[17] = spatial_metric_l[1][2];
        ql[18] = spatial_metric_l[2][0]; ql[19] = spatial_metric_l[2][1]; ql[20] = spatial_metric_l[2][2];

        ql[21] = extrinsic_curvature_l[0][0]; ql[22] = extrinsic_curvature_l[0][1]; ql[23] = extrinsic_curvature_l[0][2];
        ql[24] = extrinsic_curvature_l[1][0]; ql[25] = extrinsic_curvature_l[1][1]; ql[26] = extrinsic_curvature_l[1][2];
        ql[27] = extrinsic_curvature_l[2][0]; ql[28] = extrinsic_curvature_l[2][1]; ql[29] = extrinsic_curvature_l[2][2];

        ql[30] = 1.0;

        ql[31] = lapse_der_l[0]; ql[32] = lapse_der_l[1]; ql[33] = lapse_der_l[2];
        ql[34] = shift_der_l[0][0]; ql[35] = shift_der_l[0][1]; ql[36] = shift_der_l[0][2];
        ql[37] = shift_der_l[1][0]; ql[38] = shift_der_l[1][1]; ql[39] = shift_der_l[1][2];
        ql[40] = shift_der_l[2][0]; ql[41] = shift_der_l[2][1]; ql[42] = shift_der_l[2][2];

        ql[43] = spatial_metric_der_l[0][0][0]; ql[44] = spatial_metric_der_l[0][0][1]; ql[45] = spatial_metric_der_l[0][0][2];
        ql[46] = spatial_metric_der_l[0][1][0]; ql[47] = spatial_metric_der_l[0][1][1]; ql[48] = spatial_metric_der_l[0][1][2];
        ql[49] = spatial_metric_der_l[0][2][0]; ql[50] = spatial_metric_der_l[0][2][1]; ql[51] = spatial_metric_der_l[0][2][2];

        ql[52] = spatial_metric_der_l[1][0][0]; ql[53] = spatial_metric_der_l[1][0][1]; ql[54] = spatial_metric_der_l[1][0][2];
        ql[55] = spatial_metric_der_l[1][1][0]; ql[56] = spatial_metric_der_l[1][1][1]; ql[57] = spatial_metric_der_l[1][1][2];
        ql[58] = spatial_metric_der_l[1][2][0]; ql[59] = spatial_metric_der_l[1][2][1]; ql[60] = spatial_metric_der_l[1][2][2];

        ql[61] = spatial_metric_der_l[2][0][0]; ql[62] = spatial_metric_der_l[2][0][1]; ql[63] = spatial_metric_der_l[2][0][2];
        ql[64] = spatial_metric_der_l[2][1][0]; ql[65] = spatial_metric_der_l[2][1][1]; ql[66] = spatial_metric_der_l[2][1][2];
        ql[67] = spatial_metric_der_l[2][2][0]; ql[68] = spatial_metric_der_l[2][2][1]; ql[69] = spatial_metric_der_l[2][2][2];

        ql[70] = 0.0;
        ql[71] = x - 0.5; ql[72] = y; ql[73] = 0.0;

        qr[0] = sqrt(spatial_det_r) * rho_r * W_r;
        qr[1] = sqrt(spatial_det_r) * ((rho_r * h_star_r * (W_r * W_r) * cov_vel_r[0]) - (lapse_r * b0_r * cov_b_r[0]));
        qr[2] = sqrt(spatial_det_r) * ((rho_r * h_star_r * (W_r * W_r) * cov_vel_r[1]) - (lapse_r * b0_r * cov_b_r[1]));
        qr[3] = sqrt(spatial_det_r) * ((rho_r * h_star_r * (W_r * W_r) * cov_vel_r[2]) - (lapse_r * b0_r * cov_b_r[2]));
        qr[4] = sqrt(spatial_det_r) * ((rho_r * h_star_r * (W_r * W_r)) - p_star_r - ((lapse_r * lapse_r) * (b0_r * b0_r)) - (rho_r * W_r));

        qr[5] = sqrt(spatial_det_r) * mag_x_r;
        qr[6] = sqrt(spatial_det_r) * mag_y_r;
        qr[7] = sqrt(spatial_det_r) * mag_z_r;

        qr[8] = lapse_r;
        qr[9] = shift_r[0]; qr[10] = shift_r[1]; qr[11] = shift_r[2];

        qr[12] = spatial_metric_r[0][0]; qr[13] = spatial_metric_r[0][1]; qr[14] = spatial_metric_r[0][2];
        qr[15] = spatial_metric_r[1][0]; qr[16] = spatial_metric_r[1][1]; qr[17] = spatial_metric_r[1][2];
        qr[18] = spatial_metric_r[2][0]; qr[19] = spatial_metric_r[2][1]; qr[20] = spatial_metric_r[2][2];

        qr[21] = extrinsic_curvature_r[0][0]; qr[22] = extrinsic_curvature_r[0][1]; qr[23] = extrinsic_curvature_r[0][2];
        qr[24] = extrinsic_curvature_r[1][0]; qr[25] = extrinsic_curvature_r[1][1]; qr[26] = extrinsic_curvature_r[1][2];
        qr[27] = extrinsic_curvature_r[2][0]; qr[28] = extrinsic_curvature_r[2][1]; qr[29] = extrinsic_curvature_r[2][2];

        qr[30] = 1.0;

        qr[31] = lapse_der_r[0]; qr[32] = lapse_der_r[1]; qr[33] = lapse_der_r[2];
        qr[34] = shift_der_r[0][0]; qr[35] = shift_der_r[0][1]; qr[36] = shift_der_r[0][2];
        qr[37] = shift_der_r[1][0]; qr[38] = shift_der_r[1][1]; qr[39] = shift_der_r[1][2];
        qr[40] = shift_der_r[2][0]; qr[41] = shift_der_r[2][1]; qr[42] = shift_der_r[2][2];

        qr[43] = spatial_metric_der_r[0][0][0]; qr[44] = spatial_metric_der_r[0][0][1]; qr[45] = spatial_metric_der_r[0][0][2];
        qr[46] = spatial_metric_der_r[0][1][0]; qr[47] = spatial_metric_der_r[0][1][1]; qr[48] = spatial_metric_der_r[0][1][2];
        qr[49] = spatial_metric_der_r[0][2][0]; qr[50] = spatial_metric_der_r[0][2][1]; qr[51] = spatial_metric_der_r[0][2][2];

        qr[52] = spatial_metric_der_r[1][0][0]; qr[53] = spatial_metric_der_r[1][0][1]; qr[54] = spatial_metric_der_r[1][0][2];
        qr[55] = spatial_metric_der_r[1][1][0]; qr[56] = spatial_metric_der_r[1][1][1]; qr[57] = spatial_metric_der_r[1][1][2];
        qr[58] = spatial_metric_der_r[1][2][0]; qr[59] = spatial_metric_der_r[1][2][1]; qr[60] = spatial_metric_der_r[1][2][2];

        qr[61] = spatial_metric_der_r[2][0][0]; qr[62] = spatial_metric_der_r[2][0][1]; qr[63] = spatial_metric_der_r[2][0][2];
        qr[64] = spatial_metric_der_r[2][1][0]; qr[65] = spatial_metric_der_r[2][1][1]; qr[66] = spatial_metric_der_r[2][1][2];
        qr[67] = spatial_metric_der_r[2][2][0]; qr[68] = spatial_metric_der_r[2][2][1]; qr[69] = spatial_metric_der_r[2][2][2];

        qr[70] = 0.0;
        qr[71] = x + 0.5; qr[72] = y; qr[73] = 0.0;

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

        for (int d = 0; d < 3; d++) {
          double speeds[2], waves[2 * 74], waves_local[2 * 74];

          double ql_local[74], qr_local[74];
          gkyl_wv_eqn_rotate_to_local(gr_mhd, tau1[d], tau2[d], norm[d], ql, ql_local);
          gkyl_wv_eqn_rotate_to_local(gr_mhd, tau1[d], tau2[d], norm[d], qr, qr_local);

          double delta[74];
          for (int i = 0; i < 74; i++) {
            delta[i] = qr_local[i] - ql_local[i];
          }

          gkyl_wv_eqn_waves(gr_mhd, GKYL_WV_LOW_ORDER_FLUX, delta, ql_local, qr_local, waves_local, speeds);

          double apdq_local[74], amdq_local[74];
          gkyl_wv_eqn_qfluct(gr_mhd, GKYL_WV_LOW_ORDER_FLUX, ql_local, qr_local, waves_local, speeds, amdq_local, apdq_local);

          for (int i = 0; i < 2; i++) {
            gkyl_wv_eqn_rotate_to_global(gr_mhd, tau1[d], tau2[d], norm[d], &waves_local[i * 74], &waves[i * 74]);
          }

          double apdq[74], amdq[74];
          gkyl_wv_eqn_rotate_to_global(gr_mhd, tau1[d], tau2[d], norm[d], apdq_local, apdq);
          gkyl_wv_eqn_rotate_to_global(gr_mhd, tau1[d], tau2[d], norm[d], amdq_local, amdq);

          double fl_local[74], fr_local[74];
          gkyl_gr_mhd_flux(gas_gamma, ql_local, fl_local);
          gkyl_gr_mhd_flux(gas_gamma, qr_local, fr_local);

          double fl[74], fr[74];
          gkyl_wv_eqn_rotate_to_global(gr_mhd, tau1[d], tau2[d], norm[d], fl_local, fl);
          gkyl_wv_eqn_rotate_to_global(gr_mhd, tau1[d], tau2[d], norm[d], fr_local, fr);

          for (int i = 0; i < 74; i++) {
            TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-12) );
          }
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric_l[i]);
        gkyl_free(spatial_metric_r[i]);
        gkyl_free(extrinsic_curvature_l[i]);
        gkyl_free(extrinsic_curvature_r[i]);
        gkyl_free(shift_der_l[i]);
        gkyl_free(shift_der_r[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der_l[i][j]);
          gkyl_free(spatial_metric_der_r[i][j]);
        }
        gkyl_free(spatial_metric_der_l[i]);
        gkyl_free(spatial_metric_der_r[i]);
      }
      gkyl_free(spatial_metric_l);
      gkyl_free(spatial_metric_r);
      gkyl_free(extrinsic_curvature_l);
      gkyl_free(extrinsic_curvature_r);
      gkyl_free(shift_l);
      gkyl_free(shift_r);
      gkyl_free(vel_l);
      gkyl_free(vel_r);
      gkyl_free(mag_l);
      gkyl_free(mag_r);
      gkyl_free(cov_mag_l);
      gkyl_free(cov_mag_r);
      gkyl_free(cov_vel_l);
      gkyl_free(cov_vel_r);
      gkyl_free(spacetime_vel_l);
      gkyl_free(spacetime_vel_r);
      gkyl_free(b_l);
      gkyl_free(b_r);
      gkyl_free(cov_b_l);
      gkyl_free(cov_b_r);
      gkyl_free(lapse_der_l);
      gkyl_free(lapse_der_r);
      gkyl_free(shift_der_l);
      gkyl_free(shift_der_r);
      gkyl_free(spatial_metric_der_l);
      gkyl_free(spatial_metric_der_r);
    }
  }

  gkyl_wv_eqn_release(gr_mhd);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_mhd_waves_kerr()
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_mhd = gkyl_wv_gr_mhd_new(gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho_l = 1.0, u_l = 0.1, v_l = 0.2, w_l = 0.3, p_l = 1.5;
      double mag_x_l = 0.3, mag_y_l = 0.2, mag_z_l = 0.1;
      double rho_r = 0.1, u_r = 0.2, v_r = 0.3, w_r = 0.4, p_r = 0.15;
      double mag_x_r = 0.4, mag_y_r = 0.3, mag_z_r = 0.2;

      double spatial_det_l, spatial_det_r;
      double lapse_l, lapse_r;
      double *shift_l = gkyl_malloc(sizeof(double[3]));
      double *shift_r = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region_l, in_excision_region_r;

      double **spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
      double **spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature_l = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
        spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature_l[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der_l = gkyl_malloc(sizeof(double[3]));
      double *lapse_der_r = gkyl_malloc(sizeof(double[3]));
      double **shift_der_l = gkyl_malloc(sizeof(double*[3]));
      double **shift_der_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der_l[i] = gkyl_malloc(sizeof(double[3]));
        shift_der_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der_l = gkyl_malloc(sizeof(double**[3]));
      double ***spatial_metric_der_r = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der_l[i] = gkyl_malloc(sizeof(double*[3]));
        spatial_metric_der_r[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der_l[i][j] = gkyl_malloc(sizeof(double[3]));
          spatial_metric_der_r[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_det_l);
      spacetime->spatial_metric_det_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_det_r);
      spacetime->lapse_function_func(spacetime, 0.0, x - 0.1, y, 0.0, &lapse_l);
      spacetime->lapse_function_func(spacetime, 0.0, x + 0.1, y, 0.0, &lapse_r);
      spacetime->shift_vector_func(spacetime, 0.0, x - 0.1, y, 0.0, &shift_l);
      spacetime->shift_vector_func(spacetime, 0.0, x + 0.1, y, 0.0, &shift_r);
      spacetime->excision_region_func(spacetime, 0.0, x - 0.1, y, 0.0, &in_excision_region_l);
      spacetime->excision_region_func(spacetime, 0.0, x + 0.1, y, 0.0, &in_excision_region_r);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_metric_l);
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_metric_r);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_l);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_r);

      spacetime->lapse_function_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der_l);
      spacetime->lapse_function_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der_r);
      spacetime->shift_vector_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der_l);
      spacetime->shift_vector_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der_r);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der_l);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der_r);

      double *vel_l = gkyl_malloc(sizeof(double[3]));
      double *vel_r = gkyl_malloc(sizeof(double[3]));
      vel_l[0] = u_l; vel_l[1] = v_l; vel_l[2] = w_l;
      vel_r[0] = u_r; vel_r[1] = v_r; vel_r[2] = w_r;

      double v_sq_l = 0.0, v_sq_r = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq_l += spatial_metric_l[i][j] * vel_l[i] * vel_l[j];
          v_sq_r += spatial_metric_r[i][j] * vel_r[i] * vel_r[j];
        }
      }

      double W_l = 1.0 / sqrt(1.0 - v_sq_l);
      double W_r = 1.0 / sqrt(1.0 - v_sq_r);

      double *mag_l = gkyl_malloc(sizeof(double[3]));
      double *mag_r = gkyl_malloc(sizeof(double[3]));
      mag_l[0] = mag_x_l; mag_l[1] = mag_y_l; mag_l[2] = mag_z_l;
      mag_r[0] = mag_x_r; mag_r[1] = mag_y_r; mag_r[2] = mag_z_r;

      double *cov_mag_l = gkyl_malloc(sizeof(double[3]));
      double *cov_mag_r = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        cov_mag_l[i] = 0.0;
        cov_mag_r[i] = 0.0;

        for (int j = 0; j < 3; j++) {
          cov_mag_l[i] += spatial_metric_l[i][j] * mag_l[j];
          cov_mag_r[i] += spatial_metric_r[i][j] * mag_r[j];
        }
      }

      double *cov_vel_l = gkyl_malloc(sizeof(double[3]));
      double *cov_vel_r = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        cov_vel_l[i] = 0.0;
        cov_vel_r[i] = 0.0;

        for (int j = 0; j < 3; j++) {
          cov_vel_l[i] += spatial_metric_l[i][j] * vel_l[j];
          cov_vel_r[i] += spatial_metric_r[i][j] * vel_r[j];
        }
      }
      
      double b0_l = 0.0, b0_r = 0.0;
      for (int i = 0; i < 3; i++) {
        b0_l += W_l * mag_l[i] * (cov_vel_l[i] / lapse_l);
        b0_r += W_r * mag_r[i] * (cov_vel_r[i] / lapse_r);
      }

      double *spacetime_vel_l = gkyl_malloc(sizeof(double[4]));
      double *spacetime_vel_r = gkyl_malloc(sizeof(double[4]));
      spacetime_vel_l[0] = W_l / lapse_l;
      spacetime_vel_r[0] = W_r / lapse_r;
      for (int i = 0; i < 3; i++) {
        spacetime_vel_l[i + 1] = (W_l * vel_l[i]) - (shift_l[i] * (W_l / lapse_l));
        spacetime_vel_r[i + 1] = (W_r * vel_r[i]) - (shift_r[i] * (W_r / lapse_r));
      }

      double *b_l = gkyl_malloc(sizeof(double[3]));
      double *b_r = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        b_l[i] = (mag_l[i] + (lapse_l * b0_l * spacetime_vel_l[i + 1])) / W_l;
        b_r[i] = (mag_r[i] + (lapse_r * b0_r * spacetime_vel_r[i + 1])) / W_r;
      }

      double b_sq_l = 0.0, b_sq_r = 0.0;
      for (int i = 0; i < 3; i++) {
        b_sq_l += (mag_l[i] * cov_mag_l[i]) / (W_l * W_l);
        b_sq_r += (mag_r[i] * cov_mag_r[i]) / (W_r * W_r);
      }
      b_sq_l += ((lapse_l * lapse_l) * (b0_l * b0_l)) / (W_l * W_l);
      b_sq_r += ((lapse_r * lapse_r) * (b0_r * b0_r)) / (W_r * W_r);

      double *cov_b_l = gkyl_malloc(sizeof(double[3]));
      double *cov_b_r = gkyl_malloc(sizeof(double[3]));
      for (int i = 0; i < 3; i++) {
        cov_b_l[i] = 0.0;
        cov_b_r[i] = 0.0;

        for (int j = 0; j < 3; j++) {
          cov_b_l[i] += spatial_metric_l[i][j] * b_l[j];
          cov_b_r[i] += spatial_metric_r[i][j] * b_r[j];
        }
      }

      double h_star_l = 1.0 + ((p_l / rho_l) * (gas_gamma / (gas_gamma - 1.0))) + (b_sq_l / rho_l);
      double h_star_r = 1.0 + ((p_r / rho_r) * (gas_gamma / (gas_gamma - 1.0))) + (b_sq_r / rho_r);
      double p_star_l = p_l + (0.5 * b_sq_l);
      double p_star_r = p_r + (0.5 * b_sq_r);

      if (!in_excision_region_l && !in_excision_region_r) {
        double ql[74], qr[74];
        ql[0] = sqrt(spatial_det_l) * rho_l * W_l;
        ql[1] = sqrt(spatial_det_l) * ((rho_l * h_star_l * (W_l * W_l) * cov_vel_l[0]) - (lapse_l * b0_l * cov_b_l[0]));
        ql[2] = sqrt(spatial_det_l) * ((rho_l * h_star_l * (W_l * W_l) * cov_vel_l[1]) - (lapse_l * b0_l * cov_b_l[1]));
        ql[3] = sqrt(spatial_det_l) * ((rho_l * h_star_l * (W_l * W_l) * cov_vel_l[2]) - (lapse_l * b0_l * cov_b_l[2]));
        ql[4] = sqrt(spatial_det_l) * ((rho_l * h_star_l * (W_l * W_l)) - p_star_l - ((lapse_l * lapse_l) * (b0_l * b0_l)) - (rho_l * W_l));

        ql[5] = sqrt(spatial_det_l) * mag_x_l;
        ql[6] = sqrt(spatial_det_l) * mag_y_l;
        ql[7] = sqrt(spatial_det_l) * mag_z_l;

        ql[8] = lapse_l;
        ql[9] = shift_l[0]; ql[10] = shift_l[1]; ql[11] = shift_l[2];

        ql[12] = spatial_metric_l[0][0]; ql[13] = spatial_metric_l[0][1]; ql[14] = spatial_metric_l[0][2];
        ql[15] = spatial_metric_l[1][0]; ql[16] = spatial_metric_l[1][1]; ql[17] = spatial_metric_l[1][2];
        ql[18] = spatial_metric_l[2][0]; ql[19] = spatial_metric_l[2][1]; ql[20] = spatial_metric_l[2][2];

        ql[21] = extrinsic_curvature_l[0][0]; ql[22] = extrinsic_curvature_l[0][1]; ql[23] = extrinsic_curvature_l[0][2];
        ql[24] = extrinsic_curvature_l[1][0]; ql[25] = extrinsic_curvature_l[1][1]; ql[26] = extrinsic_curvature_l[1][2];
        ql[27] = extrinsic_curvature_l[2][0]; ql[28] = extrinsic_curvature_l[2][1]; ql[29] = extrinsic_curvature_l[2][2];

        ql[30] = 1.0;

        ql[31] = lapse_der_l[0]; ql[32] = lapse_der_l[1]; ql[33] = lapse_der_l[2];
        ql[34] = shift_der_l[0][0]; ql[35] = shift_der_l[0][1]; ql[36] = shift_der_l[0][2];
        ql[37] = shift_der_l[1][0]; ql[38] = shift_der_l[1][1]; ql[39] = shift_der_l[1][2];
        ql[40] = shift_der_l[2][0]; ql[41] = shift_der_l[2][1]; ql[42] = shift_der_l[2][2];

        ql[43] = spatial_metric_der_l[0][0][0]; ql[44] = spatial_metric_der_l[0][0][1]; ql[45] = spatial_metric_der_l[0][0][2];
        ql[46] = spatial_metric_der_l[0][1][0]; ql[47] = spatial_metric_der_l[0][1][1]; ql[48] = spatial_metric_der_l[0][1][2];
        ql[49] = spatial_metric_der_l[0][2][0]; ql[50] = spatial_metric_der_l[0][2][1]; ql[51] = spatial_metric_der_l[0][2][2];

        ql[52] = spatial_metric_der_l[1][0][0]; ql[53] = spatial_metric_der_l[1][0][1]; ql[54] = spatial_metric_der_l[1][0][2];
        ql[55] = spatial_metric_der_l[1][1][0]; ql[56] = spatial_metric_der_l[1][1][1]; ql[57] = spatial_metric_der_l[1][1][2];
        ql[58] = spatial_metric_der_l[1][2][0]; ql[59] = spatial_metric_der_l[1][2][1]; ql[60] = spatial_metric_der_l[1][2][2];

        ql[61] = spatial_metric_der_l[2][0][0]; ql[62] = spatial_metric_der_l[2][0][1]; ql[63] = spatial_metric_der_l[2][0][2];
        ql[64] = spatial_metric_der_l[2][1][0]; ql[65] = spatial_metric_der_l[2][1][1]; ql[66] = spatial_metric_der_l[2][1][2];
        ql[67] = spatial_metric_der_l[2][2][0]; ql[68] = spatial_metric_der_l[2][2][1]; ql[69] = spatial_metric_der_l[2][2][2];

        ql[70] = 0.0;
        ql[71] = x - 0.5; ql[72] = y; ql[73] = 0.0;

        qr[0] = sqrt(spatial_det_r) * rho_r * W_r;
        qr[1] = sqrt(spatial_det_r) * ((rho_r * h_star_r * (W_r * W_r) * cov_vel_r[0]) - (lapse_r * b0_r * cov_b_r[0]));
        qr[2] = sqrt(spatial_det_r) * ((rho_r * h_star_r * (W_r * W_r) * cov_vel_r[1]) - (lapse_r * b0_r * cov_b_r[1]));
        qr[3] = sqrt(spatial_det_r) * ((rho_r * h_star_r * (W_r * W_r) * cov_vel_r[2]) - (lapse_r * b0_r * cov_b_r[2]));
        qr[4] = sqrt(spatial_det_r) * ((rho_r * h_star_r * (W_r * W_r)) - p_star_r - ((lapse_r * lapse_r) * (b0_r * b0_r)) - (rho_r * W_r));

        qr[5] = sqrt(spatial_det_r) * mag_x_r;
        qr[6] = sqrt(spatial_det_r) * mag_y_r;
        qr[7] = sqrt(spatial_det_r) * mag_z_r;

        qr[8] = lapse_r;
        qr[9] = shift_r[0]; qr[10] = shift_r[1]; qr[11] = shift_r[2];

        qr[12] = spatial_metric_r[0][0]; qr[13] = spatial_metric_r[0][1]; qr[14] = spatial_metric_r[0][2];
        qr[15] = spatial_metric_r[1][0]; qr[16] = spatial_metric_r[1][1]; qr[17] = spatial_metric_r[1][2];
        qr[18] = spatial_metric_r[2][0]; qr[19] = spatial_metric_r[2][1]; qr[20] = spatial_metric_r[2][2];

        qr[21] = extrinsic_curvature_r[0][0]; qr[22] = extrinsic_curvature_r[0][1]; qr[23] = extrinsic_curvature_r[0][2];
        qr[24] = extrinsic_curvature_r[1][0]; qr[25] = extrinsic_curvature_r[1][1]; qr[26] = extrinsic_curvature_r[1][2];
        qr[27] = extrinsic_curvature_r[2][0]; qr[28] = extrinsic_curvature_r[2][1]; qr[29] = extrinsic_curvature_r[2][2];

        qr[30] = 1.0;

        qr[31] = lapse_der_r[0]; qr[32] = lapse_der_r[1]; qr[33] = lapse_der_r[2];
        qr[34] = shift_der_r[0][0]; qr[35] = shift_der_r[0][1]; qr[36] = shift_der_r[0][2];
        qr[37] = shift_der_r[1][0]; qr[38] = shift_der_r[1][1]; qr[39] = shift_der_r[1][2];
        qr[40] = shift_der_r[2][0]; qr[41] = shift_der_r[2][1]; qr[42] = shift_der_r[2][2];

        qr[43] = spatial_metric_der_r[0][0][0]; qr[44] = spatial_metric_der_r[0][0][1]; qr[45] = spatial_metric_der_r[0][0][2];
        qr[46] = spatial_metric_der_r[0][1][0]; qr[47] = spatial_metric_der_r[0][1][1]; qr[48] = spatial_metric_der_r[0][1][2];
        qr[49] = spatial_metric_der_r[0][2][0]; qr[50] = spatial_metric_der_r[0][2][1]; qr[51] = spatial_metric_der_r[0][2][2];

        qr[52] = spatial_metric_der_r[1][0][0]; qr[53] = spatial_metric_der_r[1][0][1]; qr[54] = spatial_metric_der_r[1][0][2];
        qr[55] = spatial_metric_der_r[1][1][0]; qr[56] = spatial_metric_der_r[1][1][1]; qr[57] = spatial_metric_der_r[1][1][2];
        qr[58] = spatial_metric_der_r[1][2][0]; qr[59] = spatial_metric_der_r[1][2][1]; qr[60] = spatial_metric_der_r[1][2][2];

        qr[61] = spatial_metric_der_r[2][0][0]; qr[62] = spatial_metric_der_r[2][0][1]; qr[63] = spatial_metric_der_r[2][0][2];
        qr[64] = spatial_metric_der_r[2][1][0]; qr[65] = spatial_metric_der_r[2][1][1]; qr[66] = spatial_metric_der_r[2][1][2];
        qr[67] = spatial_metric_der_r[2][2][0]; qr[68] = spatial_metric_der_r[2][2][1]; qr[69] = spatial_metric_der_r[2][2][2];

        qr[70] = 0.0;
        qr[71] = x + 0.5; qr[72] = y; qr[73] = 0.0;

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

        for (int d = 0; d < 3; d++) {
          double speeds[2], waves[2 * 74], waves_local[2 * 74];

          double ql_local[74], qr_local[74];
          gkyl_wv_eqn_rotate_to_local(gr_mhd, tau1[d], tau2[d], norm[d], ql, ql_local);
          gkyl_wv_eqn_rotate_to_local(gr_mhd, tau1[d], tau2[d], norm[d], qr, qr_local);

          double delta[74];
          for (int i = 0; i < 74; i++) {
            delta[i] = qr_local[i] - ql_local[i];
          }

          gkyl_wv_eqn_waves(gr_mhd, GKYL_WV_LOW_ORDER_FLUX, delta, ql_local, qr_local, waves_local, speeds);

          double apdq_local[74], amdq_local[74];
          gkyl_wv_eqn_qfluct(gr_mhd, GKYL_WV_LOW_ORDER_FLUX, ql_local, qr_local, waves_local, speeds, amdq_local, apdq_local);

          for (int i = 0; i < 2; i++) {
            gkyl_wv_eqn_rotate_to_global(gr_mhd, tau1[d], tau2[d], norm[d], &waves_local[i * 74], &waves[i * 74]);
          }

          double apdq[74], amdq[74];
          gkyl_wv_eqn_rotate_to_global(gr_mhd, tau1[d], tau2[d], norm[d], apdq_local, apdq);
          gkyl_wv_eqn_rotate_to_global(gr_mhd, tau1[d], tau2[d], norm[d], amdq_local, amdq);

          double fl_local[74], fr_local[74];
          gkyl_gr_mhd_flux(gas_gamma, ql_local, fl_local);
          gkyl_gr_mhd_flux(gas_gamma, qr_local, fr_local);

          double fl[74], fr[74];
          gkyl_wv_eqn_rotate_to_global(gr_mhd, tau1[d], tau2[d], norm[d], fl_local, fl);
          gkyl_wv_eqn_rotate_to_global(gr_mhd, tau1[d], tau2[d], norm[d], fr_local, fr);

          for (int i = 0; i < 74; i++) {
            TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-12) );
          }
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric_l[i]);
        gkyl_free(spatial_metric_r[i]);
        gkyl_free(extrinsic_curvature_l[i]);
        gkyl_free(extrinsic_curvature_r[i]);
        gkyl_free(shift_der_l[i]);
        gkyl_free(shift_der_r[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der_l[i][j]);
          gkyl_free(spatial_metric_der_r[i][j]);
        }
        gkyl_free(spatial_metric_der_l[i]);
        gkyl_free(spatial_metric_der_r[i]);
      }
      gkyl_free(spatial_metric_l);
      gkyl_free(spatial_metric_r);
      gkyl_free(extrinsic_curvature_l);
      gkyl_free(extrinsic_curvature_r);
      gkyl_free(shift_l);
      gkyl_free(shift_r);
      gkyl_free(vel_l);
      gkyl_free(vel_r);
      gkyl_free(mag_l);
      gkyl_free(mag_r);
      gkyl_free(cov_mag_l);
      gkyl_free(cov_mag_r);
      gkyl_free(cov_vel_l);
      gkyl_free(cov_vel_r);
      gkyl_free(spacetime_vel_l);
      gkyl_free(spacetime_vel_r);
      gkyl_free(b_l);
      gkyl_free(b_r);
      gkyl_free(cov_b_l);
      gkyl_free(cov_b_r);
      gkyl_free(lapse_der_l);
      gkyl_free(lapse_der_r);
      gkyl_free(shift_der_l);
      gkyl_free(shift_der_r);
      gkyl_free(spatial_metric_der_l);
      gkyl_free(spatial_metric_der_r);
    }
  }

  gkyl_wv_eqn_release(gr_mhd);
  gkyl_gr_spacetime_release(spacetime);
}

TEST_LIST = {
  { "gr_mhd_basic_minkowski" , test_gr_mhd_basic_minkowski},
  { "gr_mhd_basic_schwarzschild", test_gr_mhd_basic_schwarzschild },
  { "gr_mhd_basic_kerr", test_gr_mhd_basic_kerr },
  { "gr_mhd_waves_minkowski", test_gr_mhd_waves_minkowski },
  { "gr_mhd_waves_schwarzschild", test_gr_mhd_waves_schwarzschild },
  { "gr_mhd_waves_kerr", test_gr_mhd_waves_kerr },
  { NULL, NULL },
};