#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_mhd.h>
#include <gkyl_wv_gr_mhd_priv.h>

void
gkyl_gr_mhd_flux(double gas_gamma, const double q[74], double flux[74])
{
  double v[74] = { 0.0 };
  gkyl_gr_mhd_prim_vars(gas_gamma, q, v);
  double rho =  v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double mag_x = v[5];
  double mag_y = v[6];
  double mag_z = v[7];

  double lapse = v[8];
  double shift_x = v[9];
  double shift_y = v[10];
  double shift_z = v[11];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[12]; spatial_metric[0][1] = v[13]; spatial_metric[0][2] = v[14];
  spatial_metric[1][0] = v[15]; spatial_metric[1][1] = v[16]; spatial_metric[1][2] = v[17];
  spatial_metric[2][0] = v[18]; spatial_metric[2][1] = v[19]; spatial_metric[2][2] = v[20];

  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

  bool in_excision_region = false;
  if (v[30] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double vel[3];
    double v_sq = 0.0;
    vel[0] = vx; vel[1] = vy; vel[2] = vz;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq += spatial_metric[i][j] * vel[i] * vel[j];
      }
    }

    double W = 1.0 / (sqrt(1.0 - v_sq));
    if (v_sq > 1.0 - pow(10.0, -8.0)) {
      W = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double mag[3];
    mag[0] = mag_x; mag[1] = mag_y; mag[2] = mag_z;

    double cov_mag[3];
    for (int i = 0; i < 3; i++) {
      cov_mag[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_mag[i] += spatial_metric[i][j] * mag[j];
      }
    }

    double cov_vel[3];
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

    double shift[3];
    shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

    double spacetime_vel[4];
    spacetime_vel[0] = W / lapse;
    for (int i = 0; i < 3; i++) {
      spacetime_vel[i + 1] = (W * vel[i]) - (shift[i] * (W / lapse));
    }

    double b[3];
    for (int i = 0; i < 3; i++) {
      b[i] = (mag[i] + (lapse * b0 * spacetime_vel[i + 1])) / W;
    }

    double b_sq = 0.0;
    for (int i = 0; i < 3; i++) {
      b_sq += (mag[i] * cov_mag[i]) / (W * W);
    }
    b_sq += ((lapse * lapse) * (b0 * b0)) / (W * W);

    double cov_b[3];
    for (int i = 0; i < 3; i++) {
      cov_b[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_b[i] += spatial_metric[i][j] * b[j];
      }
    }

    double h_star = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0))) + (b_sq / rho);
    double p_star = p + (0.5 * b_sq);

    double D = rho * W;
    double Sx = (rho * h_star * (W * W) * cov_vel[0]) - (lapse * b0 * cov_b[0]);
    double Sy = (rho * h_star * (W * W) * cov_vel[1]) - (lapse * b0 * cov_b[1]);
    double Sz = (rho * h_star * (W * W) * cov_vel[2]) - (lapse * b0 * cov_b[2]);
    double Etot = (rho * h_star * (W * W)) - p_star - ((lapse * lapse) * (b0 * b0)) - (rho * W);

    flux[0] = (lapse * sqrt(spatial_det)) * (D * (vx - (shift_x / lapse)));
    flux[1] = (lapse * sqrt(spatial_det)) * ((Sx * (vx - (shift_x / lapse))) + p_star - ((cov_b[0] * mag_x) / W));
    flux[2] = (lapse * sqrt(spatial_det)) * ((Sy * (vx - (shift_x / lapse))) - ((cov_b[1] * mag_x) / W));
    flux[3] = (lapse * sqrt(spatial_det)) * ((Sz * (vx - (shift_x / lapse))) - ((cov_b[2] * mag_x) / W));
    flux[4] = (lapse * sqrt(spatial_det)) * ((Etot * (vx - (shift_x / lapse))) + (p_star * vx) - ((lapse * b0 * mag_x) / W));

    flux[5] = (lapse * sqrt(spatial_det)) * (((vx - (shift_x / lapse)) * mag_x) - ((vx - (shift_x / lapse)) * mag_x));
    flux[6] = (lapse * sqrt(spatial_det)) * (((vx - (shift_x / lapse)) * mag_y) - ((vy - (shift_y / lapse)) * mag_x));
    flux[7] = (lapse * sqrt(spatial_det)) * (((vx - (shift_x / lapse)) * mag_z) - ((vz - (shift_z / lapse)) * mag_x));

    for (int i = 8; i < 74; i++) {
      flux[i] = 0.0;
    }
  }
  else {
    for (int i = 0; i < 74; i++) {
      flux[i] = 0.0;
    }
  }
}

void
gkyl_gr_mhd_prim_vars(double gas_gamma, const double q[74], double v[74])
{
  double lapse = q[8];
  double shift_x = q[9];
  double shift_y = q[10];
  double shift_z = q[11];

  double spatial_metric[3][3];
  spatial_metric[0][0] = q[12]; spatial_metric[0][1] = q[13]; spatial_metric[0][2] = q[14];
  spatial_metric[1][0] = q[15]; spatial_metric[1][1] = q[16]; spatial_metric[1][2] = q[17];
  spatial_metric[2][0] = q[18]; spatial_metric[2][1] = q[19]; spatial_metric[2][2] = q[20];

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_mhd_inv_spatial_metric(q, &inv_spatial_metric);
  
  double extrinsic_curvature[3][3];
  extrinsic_curvature[0][0] = q[21]; extrinsic_curvature[0][1] = q[22]; extrinsic_curvature[0][2] = q[23];
  extrinsic_curvature[1][0] = q[24]; extrinsic_curvature[1][1] = q[25]; extrinsic_curvature[1][2] = q[26];
  extrinsic_curvature[2][0] = q[27]; extrinsic_curvature[2][1] = q[28]; extrinsic_curvature[2][2] = q[29];

  double lapse_der[3];
  lapse_der[0] = q[31];
  lapse_der[1] = q[32];
  lapse_der[2] = q[33];

  double shift_der[3][3];
  shift_der[0][0] = q[34]; shift_der[0][1] = q[35]; shift_der[0][2] = q[36];
  shift_der[1][0] = q[37]; shift_der[1][1] = q[38]; shift_der[1][2] = q[39];
  shift_der[2][0] = q[40]; shift_der[2][1] = q[41]; shift_der[2][2] = q[42];

  double spatial_metric_der[3][3][3];
  spatial_metric_der[0][0][0] = q[43]; spatial_metric_der[0][0][1] = q[44]; spatial_metric_der[0][0][2] = q[45];
  spatial_metric_der[0][1][0] = q[46]; spatial_metric_der[0][1][1] = q[47]; spatial_metric_der[0][1][2] = q[48];
  spatial_metric_der[0][2][0] = q[49]; spatial_metric_der[0][2][1] = q[50]; spatial_metric_der[0][2][2] = q[51];

  spatial_metric_der[1][0][0] = q[52]; spatial_metric_der[1][0][1] = q[53]; spatial_metric_der[1][0][2] = q[54];
  spatial_metric_der[1][1][0] = q[55]; spatial_metric_der[1][1][1] = q[56]; spatial_metric_der[1][1][2] = q[57];
  spatial_metric_der[1][2][0] = q[58]; spatial_metric_der[1][2][1] = q[59]; spatial_metric_der[1][2][2] = q[60];

  spatial_metric_der[0][0][0] = q[61]; spatial_metric_der[0][0][1] = q[62]; spatial_metric_der[0][0][2] = q[63];
  spatial_metric_der[0][1][0] = q[64]; spatial_metric_der[0][1][1] = q[65]; spatial_metric_der[0][1][2] = q[66];
  spatial_metric_der[0][2][0] = q[67]; spatial_metric_der[0][2][1] = q[68]; spatial_metric_der[0][2][2] = q[69];

  double evol_param = q[70];
  double x = q[71];
  double y = q[72];
  double z = q[73];

  bool in_excision_region = false;
  if (q[30] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
      (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
      (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

    double D = q[0] / sqrt(spatial_det);
    double momx = q[1] / sqrt(spatial_det);
    double momy = q[2] / sqrt(spatial_det);
    double momz = q[3] / sqrt(spatial_det);
    double Etot = q[4] / sqrt(spatial_det);

    double mag_x = q[5] / sqrt(spatial_det);
    double mag_y = q[6] / sqrt(spatial_det);
    double mag_z = q[7] / sqrt(spatial_det);

    double cov_mom[3];
    cov_mom[0] = momx; cov_mom[1] = momy; cov_mom[2] = momz;

    double mom[3];
    for (int i = 0; i < 3; i++) {
      mom[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        mom[i] += inv_spatial_metric[i][j] * cov_mom[j];
      }
    }

    double M_sq = 0.0;
    for (int i = 0; i < 3; i++) {
      M_sq += (mom[i] * cov_mom[i]);
    }

    double mag[3];
    mag[0] = mag_x; mag[1] = mag_y; mag[2] = mag_z;

    double cov_mag[3];
    for (int i = 0; i < 3; i++) {
      cov_mag[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_mag[i] += spatial_metric[i][j] * mag[j];
      }
    }
    
    double mag_sq = 0.0;
    for (int i = 0; i < 3; i++) {
      mag_sq += mag[i] * cov_mag[i];
    }

    double tau_star = 0.0;
    for (int i = 0; i < 3; i++) {
      tau_star += (mag[i] * cov_mom[i]);
    }

    double p_guess = 0.0, p_new = 0.0, rho = 0.0, z = 0.0;
    int iter = 0;

    while (iter < 100) {
      double a = Etot + D + p_guess + (0.5 * mag_sq);
      double d = 0.5 * ((M_sq * mag_sq) - (tau_star * tau_star));

      double phi = acos((1.0 / a) * sqrt((27.0 * d) / (4.0 * a)));
      double epsilon1 = ((1.0 / 3.0) * a) - ((2.0 / 3.0) * a * cos(((2.0 / 3.0) * phi) + ((2.0 / 3.0) * M_PI)));
      z = epsilon1 - mag_sq;

      double v_sq = ((M_sq * (z * z)) + ((tau_star * tau_star) * (mag_sq + (2.0 * z)))) / ((z * z) * (mag_sq + z) * (mag_sq + z));
      double W = 1.0 / sqrt(1.0 - v_sq);
      rho = D / W;
      double h = z / (W * W * rho);

      p_new = rho * (h - 1.0) * ((gas_gamma - 1.0) / gas_gamma);

      if (fabs(p_guess - p_new) < pow(10.0, -15.0)) {
        iter = 100;
      }
      else {
        iter += 1;
        p_guess = p_new;
      }
    }

    if (p_new < pow(10.0, -8.0)) {
      p_new = pow(10.0, -8.0);
    }
    if (rho < pow(10.0, -8.0)) {
      rho = pow(10.0, -8.0);
    }

    double vel[3];
    for (int i = 0; i < 3; i++) {
      vel[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        vel[i] += (inv_spatial_metric[i][j] * cov_mom[j]) / (z + mag_sq);
        vel[i] += ((mag[j] * cov_mom[j]) * mag[i]) / (z * (z + mag_sq));
      }
    }

    v[0] = rho;
    v[1] = vel[0];
    v[2] = vel[1];
    v[3] = vel[2];
    v[4] = p_new;

    v[5] = mag_x;
    v[6] = mag_y;
    v[7] = mag_z;

    v[8] = lapse;
    v[9] = shift_x;
    v[10] = shift_y;
    v[11] = shift_z;

    v[12] = spatial_metric[0][0]; v[13] = spatial_metric[0][1]; v[14] = spatial_metric[0][2];
    v[15] = spatial_metric[1][0]; v[16] = spatial_metric[1][1]; v[17] = spatial_metric[1][2];
    v[18] = spatial_metric[2][0]; v[19] = spatial_metric[2][1]; v[20] = spatial_metric[2][2];

    v[21] = extrinsic_curvature[0][0]; v[22] = extrinsic_curvature[0][1]; v[23] = extrinsic_curvature[0][2];
    v[24] = extrinsic_curvature[1][0]; v[25] = extrinsic_curvature[1][1]; v[26] = extrinsic_curvature[1][2];
    v[27] = extrinsic_curvature[2][0]; v[28] = extrinsic_curvature[2][1]; v[29] = extrinsic_curvature[2][2];

    v[30] = 1.0;

    v[31] = lapse_der[0];
    v[32] = lapse_der[1];
    v[33] = lapse_der[2];

    v[34] = shift_der[0][0]; v[35] = shift_der[0][1]; v[36] = shift_der[0][2];
    v[37] = shift_der[1][0]; v[38] = shift_der[1][1]; v[39] = shift_der[1][2];
    v[40] = shift_der[2][0]; v[41] = shift_der[2][1]; v[42] = shift_der[2][2];

    v[43] = spatial_metric_der[0][0][0]; v[44] = spatial_metric_der[0][0][1]; v[45] = spatial_metric_der[0][0][2];
    v[46] = spatial_metric_der[0][1][0]; v[47] = spatial_metric_der[0][1][1]; v[48] = spatial_metric_der[0][1][2];
    v[49] = spatial_metric_der[0][2][0]; v[50] = spatial_metric_der[0][2][1]; v[51] = spatial_metric_der[0][2][2];

    v[52] = spatial_metric_der[1][0][0]; v[53] = spatial_metric_der[1][0][1]; v[54] = spatial_metric_der[1][0][2];
    v[55] = spatial_metric_der[1][1][0]; v[56] = spatial_metric_der[1][1][1]; v[57] = spatial_metric_der[1][1][2];
    v[58] = spatial_metric_der[1][2][0]; v[59] = spatial_metric_der[1][2][1]; v[60] = spatial_metric_der[1][2][2];

    v[61] = spatial_metric_der[2][0][0]; v[62] = spatial_metric_der[2][0][1]; v[63] = spatial_metric_der[2][0][2];
    v[64] = spatial_metric_der[2][1][0]; v[65] = spatial_metric_der[2][1][1]; v[66] = spatial_metric_der[2][1][2];
    v[67] = spatial_metric_der[2][2][0]; v[68] = spatial_metric_der[2][2][1]; v[69] = spatial_metric_der[2][2][2];

    v[70] = evol_param;
    v[71] = x;
    v[72] = y;
    v[73] = z;
  }
  else {
    for (int i = 0; i < 74; i++) {
      v[i] = 0.0;
    }
    
    v[30] = -1.0;
  }
}

void 
gkyl_gr_mhd_inv_spatial_metric(const double q[74], double ***inv_spatial_metric)
{
  double spatial_metric[3][3];
  spatial_metric[0][0] = q[12]; spatial_metric[0][1] = q[13]; spatial_metric[0][2] = q[14];
  spatial_metric[1][0] = q[15]; spatial_metric[1][1] = q[16]; spatial_metric[1][2] = q[17];
  spatial_metric[2][0] = q[18]; spatial_metric[2][1] = q[19]; spatial_metric[2][2] = q[20];

  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  double trace = 0.0;
  for (int i = 0; i < 3; i++) {
    trace += spatial_metric[i][i];
  }

  double spatial_metric_sq[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      spatial_metric_sq[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        spatial_metric_sq[i][j] += spatial_metric[i][k] * spatial_metric[k][j];
      }
    }
  }

  double sq_trace = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_trace += spatial_metric_sq[i][i];
  }

  double euclidean_metric[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        euclidean_metric[i][j] = 1.0;
      }
      else {
        euclidean_metric[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*inv_spatial_metric)[i][j] = (1.0 / spatial_det) *
        ((0.5 * ((trace * trace) - sq_trace) * euclidean_metric[i][j]) - (trace * spatial_metric[i][j]) + spatial_metric_sq[i][j]);
    }
  }
}

void
gkyl_gr_mhd_stress_energy_tensor(double gas_gamma, const double q[74], double ***stress_energy)
{
  double v[74] = { 0.0 };
  gkyl_gr_mhd_prim_vars(gas_gamma, q, v);
  double rho = v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double mag_x = v[5];
  double mag_y = v[6];
  double mag_z = v[7];

  double lapse = v[8];
  double shift_x = v[9];
  double shift_y = v[10];
  double shift_z = v[11];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[12]; spatial_metric[0][1] = v[13]; spatial_metric[0][2] = v[14];
  spatial_metric[1][0] = v[15]; spatial_metric[1][1] = v[16]; spatial_metric[1][2] = v[17];
  spatial_metric[2][0] = v[18]; spatial_metric[2][1] = v[19]; spatial_metric[2][2] = v[20];

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_mhd_inv_spatial_metric(q, &inv_spatial_metric);

  bool in_excision_region = false;
  if (v[30] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double vel[3];
    double v_sq = 0.0;
    vel[0] = vx; vel[1] = vy; vel[2] = vz;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq += spatial_metric[i][j] * vel[i] * vel[j];
      }
    }

    double W = 1.0 / sqrt(1.0 - v_sq);
    if (v_sq > 1.0 - pow(10.0, -8.0)) {
      W = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double mag[3];
    mag[0] = mag_x; mag[1] = mag_y; mag[2] = mag_z;

    double cov_mag[3];
    for (int i = 0; i < 3; i++) {
      cov_mag[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_mag[i] += spatial_metric[i][j] * mag[j];
      }
    }

    double cov_vel[3];
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

    double shift[3];
    shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

    double spacetime_vel[4];
    spacetime_vel[0] = W / lapse;
    for (int i = 0; i < 3; i++) {
      spacetime_vel[i + 1] = (W * vel[i]) - (shift[i] * (W / lapse));
    }

    double b[3];
    for (int i = 0; i < 3; i++) {
      b[i] = (mag[i] + (lapse * b0 * spacetime_vel[i + 1])) / W;
    }

    double b_sq = 0.0;
    for (int i = 0; i < 3; i++) {
      b_sq += (mag[i] * cov_mag[i]) / (W * W);
    }
    b_sq += ((lapse * lapse) * (b0 * b0)) / (W * W);

    double cov_b[3];
    for (int i = 0; i < 3; i++) {
      cov_b[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_b[i] += spatial_metric[i][j] * b[j];
      }
    }

    double h_star = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0))) + (b_sq / rho);
    double p_star = p + (0.5 * b_sq);

    double inv_spacetime_metric[4][4];
    inv_spacetime_metric[0][0] = - (1.0 / (lapse * lapse));
    for (int i = 0; i < 3; i++) {
      inv_spacetime_metric[0][i] = (1.0 / (lapse * lapse)) * shift[i];
      inv_spacetime_metric[i][0] = (1.0 / (lapse * lapse)) * shift[i];
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        inv_spacetime_metric[i][j] = inv_spatial_metric[i][j] - ((1.0 / (lapse * lapse)) * shift[i] * shift[j]);
      }
    }

    double spacetime_b[4];
    spacetime_b[0] = b0;
    for (int i = 0; i < 3; i++) {
      spacetime_b[i + 1] = b[i];
    }

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        (*stress_energy)[i][j] = (rho * h_star * spacetime_vel[i] * spacetime_vel[j]) + (p_star * inv_spacetime_metric[i][j]) - (spacetime_b[i] * spacetime_b[j]);
      }
    }
  }
  else {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        (*stress_energy)[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(inv_spatial_metric[i]);
  }
  gkyl_free(inv_spatial_metric);
}

static inline double
gkyl_gr_mhd_max_abs_speed(double gas_gamma, const double q[74])
{
  double v[74] = { 0.0 };
  gkyl_gr_mhd_prim_vars(gas_gamma, q, v);
  double rho =  v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double mag_x = v[5];
  double mag_y = v[6];
  double mag_z = v[7];

  double lapse = v[8];
  double shift_x = v[9];
  double shift_y = v[10];
  double shift_z = v[11];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[12]; spatial_metric[0][1] = v[13]; spatial_metric[0][2] = v[14];
  spatial_metric[1][0] = v[15]; spatial_metric[1][1] = v[16]; spatial_metric[1][2] = v[17];
  spatial_metric[2][0] = v[18]; spatial_metric[2][1] = v[19]; spatial_metric[2][2] = v[20];

  bool in_excision_region = false;
  if (v[30] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double vel[3];
    double v_sq = 0.0;
    vel[0] = vx; vel[1] = vy; vel[2] = vz;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq += spatial_metric[i][j] * vel[i] * vel[j];
      }
    }

    double W = 1.0 / sqrt(1.0 - v_sq);
    if (v_sq > 1.0 - pow(10.0, -8.0)) {
      W = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double mag[3];
    mag[0] = mag_x; mag[1] = mag_y; mag[2] = mag_z;

    double cov_mag[3];
    for (int i = 0; i < 3; i++) {
      cov_mag[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_mag[i] += spatial_metric[i][j] * mag[j];
      }
    }

    double cov_vel[3];
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

    double shift[3];
    shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

    double spacetime_vel[4];
    spacetime_vel[0] = W / lapse;
    for (int i = 0; i < 3; i++) {
      spacetime_vel[i + 1] = (W * vel[i]) - (shift[i] * (W / lapse));
    }

    double b[3];
    for (int i = 0; i < 3; i++) {
      b[i] = (mag[i] + (lapse * b0 * spacetime_vel[i + 1])) / W;
    }

    double b_sq = 0.0;
    for (int i = 0; i < 3; i++) {
      b_sq += (mag[i] * cov_mag[i]) / (W * W);
    }
    b_sq += ((lapse * lapse) * (b0 * b0)) / (W * W);

    double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));
    double C = (rho * h) + b_sq;

    double entropy_eigs[3];
    double fast_alfven_eigs[3];
    double slow_alfven_eigs[3];

    for (int i = 0; i < 3; i++) {
      entropy_eigs[i] = (lapse * vel[i] - shift[i]);

      fast_alfven_eigs[i] = (b[i] + (sqrt(C) * spacetime_vel[i + 1])) / (b0 + (sqrt(C) * spacetime_vel[0]));
      slow_alfven_eigs[i] = (b[i] - (sqrt(C) * spacetime_vel[i + 1])) / (b0 - (sqrt(C) * spacetime_vel[0]));
    }

    double max_eig = 0.0;
    for (int i = 0; i < 3; i++) {
      if (fabs(entropy_eigs[i]) > max_eig) {
        max_eig = fabs(entropy_eigs[i]);
      }
      if (fabs(fast_alfven_eigs[i]) > max_eig) {
        max_eig = fabs(fast_alfven_eigs[i]);
      }
      if (fabs(slow_alfven_eigs[i]) > max_eig) {
        max_eig = fabs(slow_alfven_eigs[i]);
      }
    }

    double num = (gas_gamma * p) / rho;
    double den = 1.0 + ((p / rho) * (gas_gamma) / (gas_gamma - 1.0));
    double c_s = sqrt(num / den);

    return fabs(v_sq) + c_s;
  }
  else {
    return pow(10.0, -8.0);
  }
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 74; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 74; i++) {
    qout[i] = win[i];
  }
}

static void
gr_mhd_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 74; i++) {
    ghost[i] = skin[i];
  }

  ghost[1] = -ghost[1];
}

static void
gr_mhd_no_slip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 1; i < 4; i++) {
    ghost[i] = -skin[i];
  }

  ghost[0] = skin[0];
  ghost[4] = skin[4];

  for (int i = 5; i < 74; i++) {
    ghost[i] = skin[i];
  }
}

static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal,
  double* GKYL_RESTRICT qlocal)
{
  qlocal[0] = qglobal[0];
  qlocal[1] = (qglobal[1] * norm[0]) + (qglobal[2] * norm[1]) + (qglobal[3] * norm[2]);
  qlocal[2] = (qglobal[1] * tau1[0]) + (qglobal[2] * tau1[1]) + (qglobal[3] * tau1[2]);
  qlocal[3] = (qglobal[1] * tau2[0]) + (qglobal[2] * tau2[1]) + (qglobal[3] * tau2[2]);
  qlocal[4] = qglobal[4];

  qlocal[5] = (qglobal[5] * norm[0]) + (qglobal[6] * norm[1]) + (qglobal[7] * norm[2]);
  qlocal[6] = (qglobal[5] * tau1[0]) + (qglobal[6] * tau1[1]) + (qglobal[7] * tau1[2]);
  qlocal[7] = (qglobal[5] * tau2[0]) + (qglobal[6] * tau2[1]) + (qglobal[7] * tau2[2]);

  qlocal[8] = qglobal[8];
  qlocal[9] = (qglobal[9] * norm[0]) + (qglobal[10] * norm[1]) + (qglobal[11] * norm[2]);
  qlocal[10] = (qglobal[9] * tau1[0]) + (qglobal[10] * tau1[1]) + (qglobal[11] * tau1[2]);
  qlocal[11] = (qglobal[9] * tau2[0]) + (qglobal[10] * tau2[1]) + (qglobal[11] * tau2[2]);

  // Temporary arrays to store rotated column vectors.
  double r1[3], r2[3], r3[3];
  r1[0] = (qglobal[12] * norm[0]) + (qglobal[13] * norm[1]) + (qglobal[14] * norm[2]);
  r1[1] = (qglobal[12] * tau1[0]) + (qglobal[13] * tau1[1]) + (qglobal[14] * tau1[2]);
  r1[2] = (qglobal[12] * tau2[0]) + (qglobal[13] * tau2[1]) + (qglobal[14] * tau2[2]);

  r2[0] = (qglobal[15] * norm[0]) + (qglobal[16] * norm[1]) + (qglobal[17] * norm[2]);
  r2[1] = (qglobal[15] * tau1[0]) + (qglobal[16] * tau1[1]) + (qglobal[17] * tau1[2]);
  r2[2] = (qglobal[15] * tau2[0]) + (qglobal[16] * tau2[1]) + (qglobal[17] * tau2[2]);

  r3[0] = (qglobal[18] * norm[0]) + (qglobal[19] * norm[1]) + (qglobal[20] * norm[2]);
  r3[1] = (qglobal[18] * tau1[0]) + (qglobal[19] * tau1[1]) + (qglobal[20] * tau1[2]);
  r3[2] = (qglobal[18] * tau2[0]) + (qglobal[19] * tau2[1]) + (qglobal[20] * tau2[2]);

  // Temporary arrays to store rotated row vectors.
  double v1[3], v2[3], v3[3];
  v1[0] = (r1[0] * norm[0]) + (r2[0] * norm[1]) + (r3[0] * norm[2]);
  v1[1] = (r1[0] * tau1[0]) + (r2[0] * tau1[1]) + (r3[0] * tau1[2]);
  v1[2] = (r1[0] * tau2[0]) + (r2[0] * tau2[1]) + (r3[0] * tau2[2]);

  v2[0] = (r1[1] * norm[0]) + (r2[1] * norm[1]) + (r3[1] * norm[2]);
  v2[1] = (r1[1] * tau1[0]) + (r2[1] * tau1[1]) + (r3[1] * tau1[2]);
  v2[2] = (r1[1] * tau2[0]) + (r2[1] * tau2[1]) + (r3[1] * tau2[2]);

  v3[0] = (r1[2] * norm[0]) + (r2[2] * norm[1]) + (r3[2] * norm[2]);
  v3[1] = (r1[2] * tau1[0]) + (r2[2] * tau1[1]) + (r3[2] * tau1[2]);
  v3[2] = (r1[2] * tau2[0]) + (r2[2] * tau2[1]) + (r3[2] * tau2[2]);

  // Rotate spatial metric tensor to local coordinate frame.
  qlocal[12] = v1[0]; qlocal[13] = v1[1]; qlocal[14] = v1[2];
  qlocal[15] = v2[0]; qlocal[16] = v2[1]; qlocal[17] = v2[2];
  qlocal[18] = v3[0]; qlocal[19] = v3[1]; qlocal[20] = v3[2];

  // Temporary arrays to store rotated extrinsic column vectors.
  double extr_r1[3], extr_r2[3], extr_r3[3];
  extr_r1[0] = (qglobal[21] * norm[0]) + (qglobal[22] * norm[1]) + (qglobal[23] * norm[2]);
  extr_r1[1] = (qglobal[21] * tau1[0]) + (qglobal[22] * tau1[1]) + (qglobal[23] * tau1[2]);
  extr_r1[2] = (qglobal[21] * tau2[0]) + (qglobal[22] * tau2[1]) + (qglobal[23] * tau2[2]);

  extr_r2[0] = (qglobal[24] * norm[0]) + (qglobal[25] * norm[1]) + (qglobal[26] * norm[2]);
  extr_r2[1] = (qglobal[24] * tau1[0]) + (qglobal[25] * tau1[1]) + (qglobal[26] * tau1[2]);
  extr_r2[2] = (qglobal[24] * tau2[0]) + (qglobal[25] * tau2[1]) + (qglobal[26] * tau2[2]);

  extr_r3[0] = (qglobal[27] * norm[0]) + (qglobal[28] * norm[1]) + (qglobal[29] * norm[2]);
  extr_r3[1] = (qglobal[27] * tau1[0]) + (qglobal[28] * tau1[1]) + (qglobal[29] * tau1[2]);
  extr_r3[2] = (qglobal[27] * tau2[0]) + (qglobal[28] * tau2[1]) + (qglobal[29] * tau2[2]);

  // Temporary arrays to store rotated extrinsic row vectors.
  double inv_v1[3], inv_v2[3], inv_v3[3];
  inv_v1[0] = (extr_r1[0] * norm[0]) + (extr_r2[0] * norm[1]) + (extr_r3[0] * norm[2]);
  inv_v1[1] = (extr_r1[0] * tau1[0]) + (extr_r2[0] * tau1[1]) + (extr_r3[0] * tau1[2]);
  inv_v1[2] = (extr_r1[0] * tau2[0]) + (extr_r2[0] * tau2[1]) + (extr_r3[0] * tau2[2]);

  inv_v2[0] = (extr_r1[1] * norm[0]) + (extr_r2[1] * norm[1]) + (extr_r3[1] * norm[2]);
  inv_v2[1] = (extr_r1[1] * tau1[0]) + (extr_r2[1] * tau1[1]) + (extr_r3[1] * tau1[2]);
  inv_v2[2] = (extr_r1[1] * tau2[0]) + (extr_r2[1] * tau2[1]) + (extr_r3[1] * tau2[2]);

  inv_v3[0] = (extr_r1[2] * norm[0]) + (extr_r2[2] * norm[1]) + (extr_r3[2] * norm[2]);
  inv_v3[1] = (extr_r1[2] * tau1[0]) + (extr_r2[2] * tau1[1]) + (extr_r3[2] * tau1[2]);
  inv_v3[2] = (extr_r1[2] * tau2[0]) + (extr_r2[2] * tau2[1]) + (extr_r3[2] * tau2[2]);

  // Rotate extrinsic curvature tensor to local coordinate frame.
  qlocal[21] = inv_v1[0]; qlocal[22] = inv_v1[1]; qlocal[23] = inv_v1[2];
  qlocal[24] = inv_v2[0]; qlocal[25] = inv_v2[1]; qlocal[26] = inv_v2[2];
  qlocal[27] = inv_v3[0]; qlocal[28] = inv_v3[1]; qlocal[29] = inv_v3[2];

  qlocal[30] = qglobal[30];

  qlocal[31] = (qglobal[31] * norm[0]) + (qglobal[32] * norm[1]) + (qglobal[33] * norm[2]);
  qlocal[32] = (qglobal[31] * tau1[0]) + (qglobal[32] * tau1[1]) + (qglobal[33] * tau1[2]);
  qlocal[33] = (qglobal[31] * tau2[0]) + (qglobal[32] * tau2[1]) + (qglobal[33] * tau2[2]);

  // Temporary arrays to store rotated shift derivative column vectors.
  double shiftder_r1[3], shiftder_r2[3], shiftder_r3[3];
  shiftder_r1[0] = (qglobal[34] * norm[0]) + (qglobal[35] * norm[1]) + (qglobal[36] * norm[2]);
  shiftder_r1[1] = (qglobal[34] * tau1[0]) + (qglobal[35] * tau1[1]) + (qglobal[36] * tau1[2]);
  shiftder_r1[2] = (qglobal[34] * tau2[0]) + (qglobal[35] * tau2[1]) + (qglobal[36] * tau2[2]);

  shiftder_r2[0] = (qglobal[37] * norm[0]) + (qglobal[38] * norm[1]) + (qglobal[39] * norm[2]);
  shiftder_r2[1] = (qglobal[37] * tau1[0]) + (qglobal[38] * tau1[1]) + (qglobal[39] * tau1[2]);
  shiftder_r2[2] = (qglobal[37] * tau2[0]) + (qglobal[38] * tau2[1]) + (qglobal[39] * tau2[2]);

  shiftder_r3[0] = (qglobal[40] * norm[0]) + (qglobal[41] * norm[1]) + (qglobal[42] * norm[2]);
  shiftder_r3[1] = (qglobal[40] * tau1[0]) + (qglobal[41] * tau1[1]) + (qglobal[42] * tau1[2]);
  shiftder_r3[2] = (qglobal[40] * tau2[0]) + (qglobal[41] * tau2[1]) + (qglobal[42] * tau2[2]);

  // Temporary arrays to store rotated shift derivative row vectors.
  double shiftder_v1[3], shiftder_v2[3], shiftder_v3[3];
  shiftder_v1[0] = (shiftder_r1[0] * norm[0]) + (shiftder_r2[0] * norm[1]) + (shiftder_r3[0] * norm[2]);
  shiftder_v1[1] = (shiftder_r1[0] * tau1[0]) + (shiftder_r2[0] * tau1[1]) + (shiftder_r3[0] * tau1[2]);
  shiftder_v1[2] = (shiftder_r1[0] * tau2[0]) + (shiftder_r2[0] * tau2[1]) + (shiftder_r3[0] * tau2[2]);

  shiftder_v2[0] = (shiftder_r1[1] * norm[0]) + (shiftder_r2[1] * norm[1]) + (shiftder_r3[1] * norm[2]);
  shiftder_v2[1] = (shiftder_r1[1] * tau1[0]) + (shiftder_r2[1] * tau1[1]) + (shiftder_r3[1] * tau1[2]);
  shiftder_v2[2] = (shiftder_r1[1] * tau2[0]) + (shiftder_r2[1] * tau2[1]) + (shiftder_r3[1] * tau2[2]);

  shiftder_v3[0] = (shiftder_r1[2] * norm[0]) + (shiftder_r2[2] * norm[1]) + (shiftder_r3[2] * norm[2]);
  shiftder_v3[1] = (shiftder_r1[2] * tau1[0]) + (shiftder_r2[2] * tau1[1]) + (shiftder_r3[2] * tau1[2]);
  shiftder_v3[2] = (shiftder_r1[2] * tau2[0]) + (shiftder_r2[2] * tau2[1]) + (shiftder_r3[2] * tau2[2]);

  // Rotate shift vector derivative to local coordinate frame.
  qlocal[34] = shiftder_v1[0]; qlocal[35] = shiftder_v1[1]; qlocal[36] = shiftder_v1[2];
  qlocal[37] = shiftder_v2[0]; qlocal[38] = shiftder_v2[1]; qlocal[39] = shiftder_v2[2];
  qlocal[40] = shiftder_v3[0]; qlocal[41] = shiftder_v3[1]; qlocal[42] = shiftder_v3[2];

  // Temporary arrays to store rotated column vectors.
  double r11[3], r12[3], r13[3];
  double r21[3], r22[3], r23[3];
  double r31[3], r32[3], r33[3];

  r11[0] = (qglobal[43] * norm[0]) + (qglobal[44] * norm[1]) + (qglobal[45] * norm[2]);
  r11[1] = (qglobal[43] * tau1[0]) + (qglobal[44] * tau1[1]) + (qglobal[45] * tau1[2]);
  r11[2] = (qglobal[43] * tau2[0]) + (qglobal[44] * tau2[1]) + (qglobal[45] * tau2[2]);

  r12[0] = (qglobal[46] * norm[0]) + (qglobal[47] * norm[1]) + (qglobal[48] * norm[2]);
  r12[1] = (qglobal[46] * tau1[0]) + (qglobal[47] * tau1[1]) + (qglobal[48] * tau1[2]);
  r12[2] = (qglobal[46] * tau2[0]) + (qglobal[47] * tau2[1]) + (qglobal[48] * tau2[2]);

  r13[0] = (qglobal[49] * norm[0]) + (qglobal[50] * norm[1]) + (qglobal[51] * norm[2]);
  r13[1] = (qglobal[49] * tau1[0]) + (qglobal[50] * tau1[1]) + (qglobal[51] * tau1[2]);
  r13[2] = (qglobal[49] * tau2[0]) + (qglobal[50] * tau2[1]) + (qglobal[51] * tau2[2]);

  r21[0] = (qglobal[52] * norm[0]) + (qglobal[53] * norm[1]) + (qglobal[54] * norm[2]);
  r21[1] = (qglobal[52] * tau1[0]) + (qglobal[53] * tau1[1]) + (qglobal[54] * tau1[2]);
  r21[2] = (qglobal[52] * tau2[0]) + (qglobal[53] * tau2[1]) + (qglobal[54] * tau2[2]);

  r22[0] = (qglobal[55] * norm[0]) + (qglobal[56] * norm[1]) + (qglobal[57] * norm[2]);
  r22[1] = (qglobal[55] * tau1[0]) + (qglobal[56] * tau1[1]) + (qglobal[57] * tau1[2]);
  r22[2] = (qglobal[55] * tau2[0]) + (qglobal[56] * tau2[1]) + (qglobal[57] * tau2[2]);

  r23[0] = (qglobal[58] * norm[0]) + (qglobal[59] * norm[1]) + (qglobal[60] * norm[2]);
  r23[1] = (qglobal[58] * tau1[0]) + (qglobal[59] * tau1[1]) + (qglobal[60] * tau1[2]);
  r23[2] = (qglobal[58] * tau2[0]) + (qglobal[59] * tau2[1]) + (qglobal[60] * tau2[2]);

  r31[0] = (qglobal[61] * norm[0]) + (qglobal[62] * norm[1]) + (qglobal[63] * norm[2]);
  r31[1] = (qglobal[61] * tau1[0]) + (qglobal[62] * tau1[1]) + (qglobal[63] * tau1[2]);
  r31[2] = (qglobal[61] * tau2[0]) + (qglobal[62] * tau2[1]) + (qglobal[63] * tau2[2]);

  r32[0] = (qglobal[64] * norm[0]) + (qglobal[65] * norm[1]) + (qglobal[66] * norm[2]);
  r32[1] = (qglobal[64] * tau1[0]) + (qglobal[65] * tau1[1]) + (qglobal[66] * tau1[2]);
  r32[2] = (qglobal[64] * tau2[0]) + (qglobal[65] * tau2[1]) + (qglobal[66] * tau2[2]);

  r33[0] = (qglobal[67] * norm[0]) + (qglobal[68] * norm[1]) + (qglobal[69] * norm[2]);
  r33[1] = (qglobal[67] * tau1[0]) + (qglobal[68] * tau1[1]) + (qglobal[69] * tau1[2]);
  r33[2] = (qglobal[67] * tau2[0]) + (qglobal[68] * tau2[1]) + (qglobal[69] * tau2[2]);

  // Temporary arrays to store rotated row vectors.
  double s11[3], s12[3], s13[3];
  double s21[3], s22[3], s23[3];
  double s31[3], s32[3], s33[3];

  s11[0] = (r11[0] * norm[0]) + (r12[0] * norm[1]) + (r13[0] * norm[2]);
  s11[1] = (r11[1] * norm[0]) + (r12[1] * norm[1]) + (r13[1] * norm[2]);
  s11[2] = (r11[2] * norm[0]) + (r12[2] * norm[1]) + (r13[2] * norm[2]);

  s12[0] = (r11[0] * tau1[0]) + (r12[0] * tau1[1]) + (r13[0] * tau1[2]);
  s12[1] = (r11[1] * tau1[0]) + (r12[1] * tau1[1]) + (r13[1] * tau1[2]);
  s12[2] = (r11[2] * tau1[0]) + (r12[2] * tau1[1]) + (r13[2] * tau1[2]);

  s13[0] = (r11[0] * tau2[0]) + (r12[0] * tau2[1]) + (r13[0] * tau2[2]);
  s13[1] = (r11[1] * tau2[0]) + (r12[1] * tau2[1]) + (r13[1] * tau2[2]);
  s13[2] = (r11[2] * tau2[0]) + (r12[2] * tau2[1]) + (r13[2] * tau2[2]);

  s21[0] = (r21[0] * norm[0]) + (r22[0] * norm[1]) + (r23[0] * norm[2]);
  s21[1] = (r21[1] * norm[0]) + (r22[1] * norm[1]) + (r23[1] * norm[2]);
  s21[2] = (r21[2] * norm[0]) + (r22[2] * norm[1]) + (r23[2] * norm[2]);

  s22[0] = (r21[0] * tau1[0]) + (r22[0] * tau1[1]) + (r23[0] * tau1[2]);
  s22[1] = (r21[1] * tau1[0]) + (r22[1] * tau1[1]) + (r23[1] * tau1[2]);
  s22[2] = (r21[2] * tau1[0]) + (r22[2] * tau1[1]) + (r23[2] * tau1[2]);

  s23[0] = (r21[0] * tau2[0]) + (r22[0] * tau2[1]) + (r23[0] * tau2[2]);
  s23[1] = (r21[1] * tau2[0]) + (r22[1] * tau2[1]) + (r23[1] * tau2[2]);
  s23[2] = (r21[2] * tau2[0]) + (r22[2] * tau2[1]) + (r23[2] * tau2[2]);

  s31[0] = (r31[0] * norm[0]) + (r32[0] * norm[1]) + (r33[0] * norm[2]);
  s31[1] = (r31[1] * norm[0]) + (r32[1] * norm[1]) + (r33[1] * norm[2]);
  s31[2] = (r31[2] * norm[0]) + (r32[2] * norm[1]) + (r33[2] * norm[2]);

  s32[0] = (r31[0] * tau1[0]) + (r32[0] * tau1[1]) + (r33[0] * tau1[2]);
  s32[1] = (r31[1] * tau1[0]) + (r32[1] * tau1[1]) + (r33[1] * tau1[2]);
  s32[2] = (r31[2] * tau1[0]) + (r32[2] * tau1[1]) + (r33[2] * tau1[2]);

  s33[0] = (r31[0] * tau2[0]) + (r32[0] * tau2[1]) + (r33[0] * tau2[2]);
  s33[1] = (r31[1] * tau2[0]) + (r32[1] * tau2[1]) + (r33[1] * tau2[2]);
  s33[2] = (r31[2] * tau2[0]) + (r32[2] * tau2[1]) + (r33[2] * tau2[2]);
  
  // Rotate spatial metric tensor derivative to local coordinate frame.
  qlocal[43] = (s11[0] * norm[0]) + (s21[0] * norm[1]) + (s31[0] * norm[2]);
  qlocal[44] = (s11[1] * norm[0]) + (s21[1] * norm[1]) + (s31[1] * norm[2]);
  qlocal[45] = (s11[2] * norm[0]) + (s21[2] * norm[1]) + (s31[2] * norm[2]);

  qlocal[46] = (s12[0] * norm[0]) + (s22[0] * norm[1]) + (s32[0] * norm[2]);
  qlocal[47] = (s12[1] * norm[0]) + (s22[1] * norm[1]) + (s32[1] * norm[2]);
  qlocal[48] = (s12[2] * norm[0]) + (s22[2] * norm[1]) + (s32[2] * norm[2]);

  qlocal[49] = (s13[0] * norm[0]) + (s23[0] * norm[1]) + (s33[0] * norm[2]);
  qlocal[50] = (s13[1] * norm[0]) + (s23[1] * norm[1]) + (s33[1] * norm[2]);
  qlocal[51] = (s13[2] * norm[0]) + (s23[2] * norm[1]) + (s33[2] * norm[2]);

  qlocal[52] = (s11[0] * tau1[0]) + (s21[0] * tau1[1]) + (s31[0] * tau1[2]);
  qlocal[53] = (s11[1] * tau1[0]) + (s21[1] * tau1[1]) + (s31[1] * tau1[2]);
  qlocal[54] = (s11[2] * tau1[0]) + (s21[2] * tau1[1]) + (s31[2] * tau1[2]);

  qlocal[55] = (s12[0] * tau1[0]) + (s22[0] * tau1[1]) + (s32[0] * tau1[2]);
  qlocal[56] = (s12[1] * tau1[0]) + (s22[1] * tau1[1]) + (s32[1] * tau1[2]);
  qlocal[57] = (s12[2] * tau1[0]) + (s22[2] * tau1[1]) + (s32[2] * tau1[2]);

  qlocal[58] = (s13[0] * tau1[0]) + (s23[0] * tau1[1]) + (s33[0] * tau1[2]);
  qlocal[59] = (s13[1] * tau1[0]) + (s23[1] * tau1[1]) + (s33[1] * tau1[2]);
  qlocal[60] = (s13[2] * tau1[0]) + (s23[2] * tau1[1]) + (s33[2] * tau1[2]);

  qlocal[61] = (s11[0] * tau2[0]) + (s21[0] * tau2[1]) + (s31[0] * tau2[2]);
  qlocal[62] = (s11[1] * tau2[0]) + (s21[1] * tau2[1]) + (s31[1] * tau2[2]);
  qlocal[63] = (s11[2] * tau2[0]) + (s21[2] * tau2[1]) + (s31[2] * tau2[2]);

  qlocal[64] = (s12[0] * tau2[0]) + (s22[0] * tau2[1]) + (s32[0] * tau2[2]);
  qlocal[65] = (s12[1] * tau2[0]) + (s22[1] * tau2[1]) + (s32[1] * tau2[2]);
  qlocal[66] = (s12[2] * tau2[0]) + (s22[2] * tau2[1]) + (s32[2] * tau2[2]);

  qlocal[67] = (s13[0] * tau2[0]) + (s23[0] * tau2[1]) + (s33[0] * tau2[2]);
  qlocal[68] = (s13[1] * tau2[0]) + (s23[1] * tau2[1]) + (s33[1] * tau2[2]);
  qlocal[69] = (s13[2] * tau2[0]) + (s23[2] * tau2[1]) + (s33[2] * tau2[2]);

  qlocal[70] = qglobal[70];
  qlocal[71] = (qglobal[71] * norm[0]) + (qglobal[72] * norm[1]) + (qglobal[73] * norm[2]);
  qlocal[72] = (qglobal[71] * tau1[0]) + (qglobal[72] * tau1[1]) + (qglobal[73] * tau1[2]);
  qlocal[73] = (qglobal[71] * tau2[0]) + (qglobal[72] * tau2[1]) + (qglobal[73] * tau2[2]);
}

static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal,
  double* GKYL_RESTRICT qglobal)
{
  qglobal[0] = qlocal[0];
  qglobal[1] = (qlocal[1] * norm[0]) + (qlocal[2] * tau1[0]) + (qlocal[3] * tau2[0]);
  qglobal[2] = (qlocal[1] * norm[1]) + (qlocal[2] * tau1[1]) + (qlocal[3] * tau2[1]);
  qglobal[3] = (qlocal[1] * norm[2]) + (qlocal[2] * tau1[2]) + (qlocal[3] * tau2[2]);
  qglobal[4] = qlocal[4];

  qglobal[5] = (qlocal[5] * norm[0]) + (qlocal[6] * tau1[0]) + (qlocal[7] * tau2[0]);
  qglobal[6] = (qlocal[5] * norm[1]) + (qlocal[6] * tau1[1]) + (qlocal[7] * tau2[1]);
  qglobal[7] = (qlocal[5] * norm[2]) + (qlocal[6] * tau1[2]) + (qlocal[7] * tau2[2]);

  qglobal[8] = qlocal[8];
  qglobal[9] = (qlocal[9] * norm[0]) + (qlocal[10] * tau1[0]) + (qlocal[11] * tau2[0]);
  qglobal[10] = (qlocal[9] * norm[1]) + (qlocal[10] * tau1[1]) + (qlocal[11] * tau2[1]);
  qglobal[11] = (qlocal[9] * norm[2]) + (qlocal[10] * tau1[2]) + (qlocal[11] * tau2[2]);

  // Temporary arrays to store rotated column vectors.
  double r1[3], r2[3], r3[3];
  r1[0] = (qlocal[12] * norm[0]) + (qlocal[13] * tau1[0]) + (qlocal[14] * tau2[0]);
  r1[1] = (qlocal[12] * norm[1]) + (qlocal[13] * tau1[1]) + (qlocal[14] * tau2[1]);
  r1[2] = (qlocal[12] * norm[2]) + (qlocal[13] * tau1[2]) + (qlocal[14] * tau2[2]);

  r2[0] = (qlocal[15] * norm[0]) + (qlocal[16] * tau1[0]) + (qlocal[17] * tau2[0]);
  r2[1] = (qlocal[15] * norm[1]) + (qlocal[16] * tau1[1]) + (qlocal[17] * tau2[1]);
  r2[2] = (qlocal[15] * norm[2]) + (qlocal[16] * tau1[2]) + (qlocal[17] * tau2[2]);

  r3[0] = (qlocal[18] * norm[0]) + (qlocal[19] * tau1[0]) + (qlocal[20] * tau2[0]);
  r3[1] = (qlocal[18] * norm[1]) + (qlocal[19] * tau1[1]) + (qlocal[20] * tau2[1]);
  r3[2] = (qlocal[18] * norm[2]) + (qlocal[19] * tau1[2]) + (qlocal[20] * tau2[2]);

  // Temporary arrays to store rotated row vectors.
  double v1[3], v2[3], v3[3];
  v1[0] = (r1[0] * norm[0]) + (r2[0] * tau1[0]) + (r3[0] * tau2[0]);
  v1[1] = (r1[0] * norm[1]) + (r2[0] * tau1[1]) + (r3[0] * tau2[1]);
  v1[2] = (r1[0] * norm[2]) + (r2[0] * tau1[2]) + (r3[0] * tau2[2]);

  v2[0] = (r1[1] * norm[0]) + (r2[1] * tau1[0]) + (r3[1] * tau2[0]);
  v2[1] = (r1[1] * norm[1]) + (r2[1] * tau1[1]) + (r3[1] * tau2[1]);
  v2[2] = (r1[1] * norm[2]) + (r2[1] * tau1[2]) + (r3[1] * tau2[2]);

  v3[0] = (r1[2] * norm[0]) + (r2[2] * tau1[0]) + (r3[2] * tau2[0]);
  v3[1] = (r1[2] * norm[1]) + (r2[2] * tau1[1]) + (r3[2] * tau2[1]);
  v3[2] = (r1[2] * norm[2]) + (r2[2] * tau1[2]) + (r3[2] * tau2[2]);

  // Rotate spatial metric tensor back to global coordinate frame.
  qglobal[12] = v1[0]; qglobal[13] = v1[1]; qglobal[14] = v1[2];
  qglobal[15] = v2[0]; qglobal[16] = v2[1]; qglobal[17] = v2[2];
  qglobal[18] = v3[0]; qglobal[19] = v3[1]; qglobal[20] = v3[2];

  // Temporary arrays to store rotated extrinsic column vectors.
  double extr_r1[3], extr_r2[3], extr_r3[3];
  extr_r1[0] = (qlocal[21] * norm[0]) + (qlocal[22] * tau1[0]) + (qlocal[23] * tau2[0]);
  extr_r1[1] = (qlocal[21] * norm[1]) + (qlocal[22] * tau1[1]) + (qlocal[23] * tau2[1]);
  extr_r1[2] = (qlocal[21] * norm[2]) + (qlocal[22] * tau1[2]) + (qlocal[23] * tau2[2]);

  extr_r2[0] = (qlocal[24] * norm[0]) + (qlocal[25] * tau1[0]) + (qlocal[26] * tau2[0]);
  extr_r2[1] = (qlocal[24] * norm[1]) + (qlocal[25] * tau1[1]) + (qlocal[26] * tau2[1]);
  extr_r2[2] = (qlocal[24] * norm[2]) + (qlocal[25] * tau1[2]) + (qlocal[26] * tau2[2]);

  extr_r3[0] = (qlocal[27] * norm[0]) + (qlocal[28] * tau1[0]) + (qlocal[29] * tau2[0]);
  extr_r3[1] = (qlocal[27] * norm[1]) + (qlocal[28] * tau1[1]) + (qlocal[29] * tau2[1]);
  extr_r3[2] = (qlocal[27] * norm[2]) + (qlocal[28] * tau1[2]) + (qlocal[29] * tau2[2]);

  // Temporary arrays to store rotated extrinsic row vectors.
  double inv_v1[3], inv_v2[3], inv_v3[3];
  inv_v1[0] = (extr_r1[0] * norm[0]) + (extr_r2[0] * tau1[0]) + (extr_r3[0] * tau2[0]);
  inv_v1[1] = (extr_r1[0] * norm[1]) + (extr_r2[0] * tau1[1]) + (extr_r3[0] * tau2[1]);
  inv_v1[2] = (extr_r1[0] * norm[2]) + (extr_r2[0] * tau1[2]) + (extr_r3[0] * tau2[2]);

  inv_v2[0] = (extr_r1[1] * norm[0]) + (extr_r2[1] * tau1[0]) + (extr_r3[1] * tau2[0]);
  inv_v2[1] = (extr_r1[1] * norm[1]) + (extr_r2[1] * tau1[1]) + (extr_r3[1] * tau2[1]);
  inv_v2[2] = (extr_r1[1] * norm[2]) + (extr_r2[1] * tau1[2]) + (extr_r3[1] * tau2[2]);

  inv_v3[0] = (extr_r1[2] * norm[0]) + (extr_r2[2] * tau1[0]) + (extr_r3[2] * tau2[0]);
  inv_v3[1] = (extr_r1[2] * norm[1]) + (extr_r2[2] * tau1[1]) + (extr_r3[2] * tau2[1]);
  inv_v3[2] = (extr_r1[2] * norm[2]) + (extr_r2[2] * tau1[2]) + (extr_r3[2] * tau2[2]);

  // Rotate extrinsic curvature tensor back to global coordinate frame.
  qglobal[21] = inv_v1[0]; qglobal[22] = inv_v1[1]; qglobal[23] = inv_v1[2];
  qglobal[24] = inv_v2[0]; qglobal[25] = inv_v2[1]; qglobal[26] = inv_v2[2];
  qglobal[27] = inv_v3[0]; qglobal[28] = inv_v3[1]; qglobal[29] = inv_v3[2];

  qglobal[30] = qlocal[30];

  qglobal[31] = (qlocal[31] * norm[0]) + (qlocal[32] * tau1[0]) + (qlocal[33] * tau2[0]);
  qglobal[32] = (qlocal[31] * norm[1]) + (qlocal[32] * tau1[1]) + (qlocal[33] * tau2[1]);
  qglobal[33] = (qlocal[31] * norm[2]) + (qlocal[32] * tau1[2]) + (qlocal[33] * tau2[2]);

  // Temporary arrays to store rotated shift derivative column vectors.
  double shiftder_r1[3], shiftder_r2[3], shiftder_r3[3];
  shiftder_r1[0] = (qlocal[34] * norm[0]) + (qlocal[35] * tau1[0]) + (qlocal[36] * tau2[0]);
  shiftder_r1[1] = (qlocal[34] * norm[1]) + (qlocal[35] * tau1[1]) + (qlocal[36] * tau2[1]);
  shiftder_r1[2] = (qlocal[34] * norm[2]) + (qlocal[35] * tau1[2]) + (qlocal[36] * tau2[2]);

  shiftder_r2[0] = (qlocal[37] * norm[0]) + (qlocal[38] * tau1[0]) + (qlocal[39] * tau2[0]);
  shiftder_r2[1] = (qlocal[37] * norm[1]) + (qlocal[38] * tau1[1]) + (qlocal[39] * tau2[1]);
  shiftder_r2[2] = (qlocal[37] * norm[2]) + (qlocal[38] * tau1[2]) + (qlocal[39] * tau2[2]);

  shiftder_r3[0] = (qlocal[40] * norm[0]) + (qlocal[41] * tau1[0]) + (qlocal[42] * tau2[0]);
  shiftder_r3[1] = (qlocal[40] * norm[1]) + (qlocal[41] * tau1[1]) + (qlocal[42] * tau2[1]);
  shiftder_r3[2] = (qlocal[40] * norm[2]) + (qlocal[41] * tau1[2]) + (qlocal[42] * tau2[2]);

  // Temporary arrays to store rotated shift derivative row vectors.
  double shiftder_v1[3], shiftder_v2[3], shiftder_v3[3];
  shiftder_v1[0] = (shiftder_r1[0] * norm[0]) + (shiftder_r2[0] * tau1[0]) + (shiftder_r3[0] * tau2[0]);
  shiftder_v1[1] = (shiftder_r1[0] * norm[1]) + (shiftder_r2[0] * tau1[1]) + (shiftder_r3[0] * tau2[1]);
  shiftder_v1[2] = (shiftder_r1[0] * norm[2]) + (shiftder_r2[0] * tau1[2]) + (shiftder_r3[0] * tau2[2]);

  shiftder_v2[0] = (shiftder_r1[1] * norm[0]) + (shiftder_r2[1] * tau1[0]) + (shiftder_r3[1] * tau2[0]);
  shiftder_v2[1] = (shiftder_r1[1] * norm[1]) + (shiftder_r2[1] * tau1[1]) + (shiftder_r3[1] * tau2[1]);
  shiftder_v2[2] = (shiftder_r1[1] * norm[2]) + (shiftder_r2[1] * tau1[2]) + (shiftder_r3[1] * tau2[2]);

  shiftder_v3[0] = (shiftder_r1[2] * norm[0]) + (shiftder_r2[2] * tau1[0]) + (shiftder_r3[2] * tau2[0]);
  shiftder_v3[1] = (shiftder_r1[2] * norm[1]) + (shiftder_r2[2] * tau1[1]) + (shiftder_r3[2] * tau2[1]);
  shiftder_v3[2] = (shiftder_r1[2] * norm[2]) + (shiftder_r2[2] * tau1[2]) + (shiftder_r3[2] * tau2[2]);

  // Rotate shift vector derivative back to global coordinate frame.
  qglobal[34] = shiftder_v1[0]; qglobal[35] = shiftder_v1[1]; qglobal[36] = shiftder_v1[2];
  qglobal[37] = shiftder_v2[0]; qglobal[38] = shiftder_v2[1]; qglobal[39] = shiftder_v2[2];
  qglobal[40] = shiftder_v3[0]; qglobal[41] = shiftder_v3[1]; qglobal[42] = shiftder_v3[2];

  // Temporary arrays to store rotated column vectors.
  double r11[3], r12[3], r13[3];
  double r21[3], r22[3], r23[3];
  double r31[3], r32[3], r33[3];

  r11[0] = (qlocal[43] * norm[0]) + (qlocal[52] * tau1[0]) + (qlocal[61] * tau2[0]);
  r11[1] = (qlocal[43] * norm[1]) + (qlocal[52] * tau1[1]) + (qlocal[61] * tau2[1]);
  r11[2] = (qlocal[43] * norm[2]) + (qlocal[52] * tau1[2]) + (qlocal[61] * tau2[2]);

  r12[0] = (qlocal[44] * norm[0]) + (qlocal[53] * tau1[0]) + (qlocal[62] * tau2[0]);
  r12[1] = (qlocal[44] * norm[1]) + (qlocal[53] * tau1[1]) + (qlocal[62] * tau2[1]);
  r12[2] = (qlocal[44] * norm[2]) + (qlocal[53] * tau1[2]) + (qlocal[62] * tau2[2]);

  r13[0] = (qlocal[45] * norm[0]) + (qlocal[54] * tau1[0]) + (qlocal[63] * tau2[0]);
  r13[1] = (qlocal[45] * norm[1]) + (qlocal[54] * tau1[1]) + (qlocal[63] * tau2[1]);
  r13[2] = (qlocal[45] * norm[2]) + (qlocal[54] * tau1[2]) + (qlocal[63] * tau2[2]);

  r21[0] = (qlocal[46] * norm[0]) + (qlocal[55] * tau1[0]) + (qlocal[64] * tau2[0]);
  r21[1] = (qlocal[46] * norm[1]) + (qlocal[55] * tau1[1]) + (qlocal[64] * tau2[1]);
  r21[2] = (qlocal[46] * norm[2]) + (qlocal[55] * tau1[2]) + (qlocal[64] * tau2[2]);

  r22[0] = (qlocal[47] * norm[0]) + (qlocal[56] * tau1[0]) + (qlocal[65] * tau2[0]);
  r22[1] = (qlocal[47] * norm[1]) + (qlocal[56] * tau1[1]) + (qlocal[65] * tau2[1]);
  r22[2] = (qlocal[47] * norm[2]) + (qlocal[56] * tau1[2]) + (qlocal[65] * tau2[2]);

  r23[0] = (qlocal[48] * norm[0]) + (qlocal[57] * tau1[0]) + (qlocal[66] * tau2[0]);
  r23[1] = (qlocal[48] * norm[1]) + (qlocal[57] * tau1[1]) + (qlocal[66] * tau2[1]);
  r23[2] = (qlocal[48] * norm[2]) + (qlocal[57] * tau1[2]) + (qlocal[66] * tau2[2]);

  r31[0] = (qlocal[49] * norm[0]) + (qlocal[58] * tau1[0]) + (qlocal[67] * tau2[0]);
  r31[1] = (qlocal[49] * norm[1]) + (qlocal[58] * tau1[1]) + (qlocal[67] * tau2[1]);
  r31[2] = (qlocal[49] * norm[2]) + (qlocal[58] * tau1[2]) + (qlocal[67] * tau2[2]);

  r32[0] = (qlocal[50] * norm[0]) + (qlocal[59] * tau1[0]) + (qlocal[68] * tau2[0]);
  r32[1] = (qlocal[50] * norm[1]) + (qlocal[59] * tau1[1]) + (qlocal[68] * tau2[1]);
  r32[2] = (qlocal[50] * norm[2]) + (qlocal[59] * tau1[2]) + (qlocal[68] * tau2[2]);

  r33[0] = (qlocal[51] * norm[0]) + (qlocal[60] * tau1[0]) + (qlocal[69] * tau2[0]);
  r33[1] = (qlocal[51] * norm[1]) + (qlocal[60] * tau1[1]) + (qlocal[69] * tau2[1]);
  r33[2] = (qlocal[51] * norm[2]) + (qlocal[60] * tau1[2]) + (qlocal[69] * tau2[2]);

  // Temporary arrays to store rotated row vectors.
  double s11[3], s12[3], s13[3];
  double s21[3], s22[3], s23[3];
  double s31[3], s32[3], s33[3];

  s11[0] = (r11[0] * norm[0]) + (r21[0] * tau1[0]) + (r31[0] * tau2[0]);
  s11[1] = (r11[1] * norm[0]) + (r21[1] * tau1[0]) + (r31[1] * tau2[0]);
  s11[2] = (r11[2] * norm[0]) + (r21[2] * tau1[0]) + (r31[2] * tau2[0]);

  s12[0] = (r11[0] * norm[1]) + (r21[0] * tau1[1]) + (r31[0] * tau2[1]);
  s12[1] = (r11[1] * norm[1]) + (r21[1] * tau1[1]) + (r31[1] * tau2[1]);
  s12[2] = (r11[2] * norm[1]) + (r21[2] * tau1[1]) + (r31[2] * tau2[1]);

  s13[0] = (r11[0] * norm[2]) + (r21[0] * tau1[2]) + (r31[0] * tau2[2]);
  s13[1] = (r11[1] * norm[2]) + (r21[1] * tau1[2]) + (r31[1] * tau2[2]);
  s13[2] = (r11[2] * norm[2]) + (r21[2] * tau1[2]) + (r31[2] * tau2[2]);

  s21[0] = (r12[0] * norm[0]) + (r22[0] * tau1[0]) + (r32[0] * tau2[0]);
  s21[1] = (r12[1] * norm[0]) + (r22[1] * tau1[0]) + (r32[1] * tau2[0]);
  s21[2] = (r12[2] * norm[0]) + (r22[2] * tau1[0]) + (r32[2] * tau2[0]);

  s22[0] = (r12[0] * norm[1]) + (r22[0] * tau1[1]) + (r32[0] * tau2[1]);
  s22[1] = (r12[1] * norm[1]) + (r22[1] * tau1[1]) + (r32[1] * tau2[1]);
  s22[2] = (r12[2] * norm[1]) + (r22[2] * tau1[1]) + (r32[2] * tau2[1]);

  s23[0] = (r12[0] * norm[2]) + (r22[0] * tau1[2]) + (r32[0] * tau2[2]);
  s23[1] = (r12[1] * norm[2]) + (r22[1] * tau1[2]) + (r32[1] * tau2[2]);
  s23[2] = (r12[2] * norm[2]) + (r22[2] * tau1[2]) + (r32[2] * tau2[2]);

  s31[0] = (r13[0] * norm[0]) + (r23[0] * tau1[0]) + (r33[0] * tau2[0]);
  s31[1] = (r13[1] * norm[0]) + (r23[1] * tau1[0]) + (r33[1] * tau2[0]);
  s31[2] = (r13[2] * norm[0]) + (r23[2] * tau1[0]) + (r33[2] * tau2[0]);

  s32[0] = (r13[0] * norm[1]) + (r23[0] * tau1[1]) + (r33[0] * tau2[1]);
  s32[1] = (r13[1] * norm[1]) + (r23[1] * tau1[1]) + (r33[1] * tau2[1]);
  s32[2] = (r13[2] * norm[1]) + (r23[2] * tau1[1]) + (r33[2] * tau2[1]);

  s33[0] = (r13[0] * norm[2]) + (r23[0] * tau1[2]) + (r33[0] * tau2[2]);
  s33[1] = (r13[1] * norm[2]) + (r23[1] * tau1[2]) + (r33[1] * tau2[2]);
  s33[2] = (r13[2] * norm[2]) + (r23[2] * tau1[2]) + (r33[2] * tau2[2]);

  // Rotate spatial metric tensor derivative back to global coordinate frame.
  qglobal[43] = (s11[0] * norm[0]) + (s12[0] * tau1[0]) + (s13[0] * tau2[0]);
  qglobal[44] = (s11[1] * norm[0]) + (s12[1] * tau1[0]) + (s13[1] * tau2[0]);
  qglobal[45] = (s11[2] * norm[0]) + (s12[2] * tau1[0]) + (s13[2] * tau2[0]);

  qglobal[46] = (s11[0] * norm[1]) + (s12[0] * tau1[1]) + (s13[0] * tau2[1]);
  qglobal[47] = (s11[1] * norm[1]) + (s12[1] * tau1[1]) + (s13[1] * tau2[1]);
  qglobal[48] = (s11[2] * norm[1]) + (s12[2] * tau1[1]) + (s13[2] * tau2[1]);

  qglobal[49] = (s11[0] * norm[2]) + (s12[0] * tau1[2]) + (s13[0] * tau2[2]);
  qglobal[50] = (s11[1] * norm[2]) + (s12[1] * tau1[2]) + (s13[1] * tau2[2]);
  qglobal[51] = (s11[2] * norm[2]) + (s12[2] * tau1[2]) + (s13[2] * tau2[2]);

  qglobal[52] = (s21[0] * norm[0]) + (s22[0] * tau1[0]) + (s23[0] * tau2[0]);
  qglobal[53] = (s21[1] * norm[0]) + (s22[1] * tau1[0]) + (s23[1] * tau2[0]);
  qglobal[54] = (s21[2] * norm[0]) + (s22[2] * tau1[0]) + (s23[2] * tau2[0]);

  qglobal[55] = (s21[0] * norm[1]) + (s22[0] * tau1[1]) + (s23[0] * tau2[1]);
  qglobal[56] = (s21[1] * norm[1]) + (s22[1] * tau1[1]) + (s23[1] * tau2[1]);
  qglobal[57] = (s21[2] * norm[1]) + (s22[2] * tau1[1]) + (s23[2] * tau2[1]);

  qglobal[58] = (s21[0] * norm[2]) + (s22[0] * tau1[2]) + (s23[0] * tau2[2]);
  qglobal[59] = (s21[1] * norm[2]) + (s22[1] * tau1[2]) + (s23[1] * tau2[2]);
  qglobal[60] = (s21[2] * norm[2]) + (s22[2] * tau1[2]) + (s23[2] * tau2[2]);

  qglobal[61] = (s31[0] * norm[0]) + (s32[0] * tau1[0]) + (s33[0] * tau2[0]);
  qglobal[62] = (s31[1] * norm[0]) + (s32[1] * tau1[0]) + (s33[1] * tau2[0]);
  qglobal[63] = (s31[2] * norm[0]) + (s32[2] * tau1[0]) + (s33[2] * tau2[0]);

  qglobal[64] = (s31[0] * norm[1]) + (s32[0] * tau1[1]) + (s33[0] * tau2[1]);
  qglobal[65] = (s31[1] * norm[1]) + (s32[1] * tau1[1]) + (s33[1] * tau2[1]);
  qglobal[66] = (s31[2] * norm[1]) + (s32[2] * tau1[1]) + (s33[2] * tau2[1]);

  qglobal[67] = (s31[0] * norm[2]) + (s32[0] * tau1[2]) + (s33[0] * tau2[2]);
  qglobal[68] = (s31[1] * norm[2]) + (s32[1] * tau1[2]) + (s33[1] * tau2[2]);
  qglobal[69] = (s31[2] * norm[2]) + (s32[2] * tau1[2]) + (s33[2] * tau2[2]);

  qglobal[70] = qlocal[70];
  qglobal[71] = (qlocal[71] * norm[0]) + (qlocal[72] * tau1[0]) + (qlocal[73] * tau2[0]);
  qglobal[72] = (qlocal[71] * norm[1]) + (qlocal[72] * tau1[1]) + (qlocal[73] * tau2[1]);
  qglobal[73] = (qlocal[71] * norm[2]) + (qlocal[72] * tau1[2]) + (qlocal[73] * tau2[2]);
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_mhd *gr_mhd = container_of(eqn, struct wv_gr_mhd, eqn);
  double gas_gamma = gr_mhd->gas_gamma;

  double sl = gkyl_gr_mhd_max_abs_speed(gas_gamma, ql);
  double sr = gkyl_gr_mhd_max_abs_speed(gas_gamma, qr);
  double amax = fmax(sl, sr);

  double fl[74], fr[74];
  gkyl_gr_mhd_flux(gas_gamma, ql, fl);
  gkyl_gr_mhd_flux(gas_gamma, qr, fr);

  bool in_excision_region_l = false;
  if (ql[30] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  bool in_excision_region_r = false;
  if (qr[30] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  double *w0 = &waves[0], *w1 = &waves[74];
  if (!in_excision_region_l && !in_excision_region_r) {
    for (int i = 0; i < 74; i++) {
      w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr[i] - fl[i]) / amax);
      w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr[i] - fl[i]) / amax);
    }
  }
  else {
    for (int i = 0; i < 74; i++) {
      w0[i] = 0.0;
      w1[i] = 0.0;
    }
  }

  s[0] = -amax;
  s[1] = amax;

  return s[1];
}

static void
qfluct_lax(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[74];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 74; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}

static double
wave_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  return wave_lax(eqn, delta, ql, qr, waves, s);
}

static void
qfluct_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double* waves, const double* s,
  double* amdq, double* apdq)
{
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static double
flux_jump(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, double* flux_jump)
{
  const struct wv_gr_mhd *gr_mhd = container_of(eqn, struct wv_gr_mhd, eqn);
  double gas_gamma = gr_mhd->gas_gamma;

  double fr[74], fl[74];
  gkyl_gr_mhd_flux(gas_gamma, ql, fl);
  gkyl_gr_mhd_flux(gas_gamma, qr, fr);

  bool in_excision_region_l = false;
  if (ql[30] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  bool in_excision_region_r = false;
  if (qr[30] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  if (!in_excision_region_l && !in_excision_region_r) {
    for (int m = 0; m < 74; m++) {
      flux_jump[m] = fr[m] - fl[m];
    }
  }
  else {
    for (int m = 0; m < 74; m++) {
      flux_jump[m] = 0.0;
    }
  }

  double amaxl = gkyl_gr_mhd_max_abs_speed(gas_gamma, ql);
  double amaxr = gkyl_gr_mhd_max_abs_speed(gas_gamma, qr);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_gr_mhd *gr_mhd = container_of(eqn, struct wv_gr_mhd, eqn);
  double gas_gamma = gr_mhd->gas_gamma;

  double v[74] = { 0.0 };
  gkyl_gr_mhd_prim_vars(gas_gamma, q, v);

  if (v[0] < 0.0 || v[4] < 0.0) {
    return false;
  }
  else {
    return true;
  }
}

static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_gr_mhd *gr_mhd = container_of(eqn, struct wv_gr_mhd, eqn);
  double gas_gamma = gr_mhd->gas_gamma;

  return gkyl_gr_mhd_max_abs_speed(gas_gamma, q);
}

static inline void
gr_mhd_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  for (int i = 0; i < 8; i++) {
    diag[i] = qin[i];
  }
}

void
gkyl_gr_mhd_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_gr_mhd *gr_mhd = container_of(base->on_dev, struct wv_gr_mhd, eqn);
    gkyl_cu_free(gr_mhd);
  }

  struct wv_gr_mhd *gr_mhd = container_of(base, struct wv_gr_mhd, eqn);
  gkyl_free(gr_mhd);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_mhd_new(double gas_gamma, enum gkyl_spacetime_gauge spacetime_gauge, int reinit_freq, struct gkyl_gr_spacetime* spacetime, bool use_gpu)
{
  return gkyl_wv_gr_mhd_inew(&(struct gkyl_wv_gr_mhd_inp) {
      .gas_gamma = gas_gamma,
      .spacetime_gauge = spacetime_gauge,
      .reinit_freq = reinit_freq,
      .spacetime = spacetime,
      .rp_type = WV_GR_MHD_RP_LAX,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_gr_mhd_inew(const struct gkyl_wv_gr_mhd_inp* inp)
{
  struct wv_gr_mhd *gr_mhd = gkyl_malloc(sizeof(struct wv_gr_mhd));

  gr_mhd->eqn.type = GKYL_EQN_GR_MHD;
  gr_mhd->eqn.num_equations = 74;
  gr_mhd->eqn.num_diag = 8;

  gr_mhd->gas_gamma = inp->gas_gamma;
  gr_mhd->spacetime_gauge = inp->spacetime_gauge;
  gr_mhd->reinit_freq = inp->reinit_freq;
  gr_mhd->spacetime = inp->spacetime;

  if (inp->rp_type == WV_GR_MHD_RP_LAX) {
    gr_mhd->eqn.num_waves = 2;
    gr_mhd->eqn.waves_func = wave_lax_l;
    gr_mhd->eqn.qfluct_func = qfluct_lax_l;
  }

  gr_mhd->eqn.flux_jump = flux_jump;
  gr_mhd->eqn.check_inv_func = check_inv;
  gr_mhd->eqn.max_speed_func = max_speed;
  gr_mhd->eqn.rotate_to_local_func = rot_to_local;
  gr_mhd->eqn.rotate_to_global_func = rot_to_global;

  gr_mhd->eqn.wall_bc_func = gr_mhd_wall;
  gr_mhd->eqn.no_slip_bc_func = gr_mhd_no_slip;

  gr_mhd->eqn.cons_to_riem = cons_to_riem;
  gr_mhd->eqn.riem_to_cons = riem_to_cons;

  gr_mhd->eqn.cons_to_diag = gr_mhd_cons_to_diag;

  gr_mhd->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_mhd->eqn.flags);
  gr_mhd->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_mhd_free);
  gr_mhd->eqn.on_dev = &gr_mhd->eqn; // On the CPU, the equation object points to itself.

  return &gr_mhd->eqn;
}

double
gkyl_wv_gr_mhd_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_mhd *gr_mhd = container_of(eqn, struct wv_gr_mhd, eqn);
  double gas_gamma = gr_mhd->gas_gamma;

  return gas_gamma;
}

enum gkyl_spacetime_gauge
gkyl_wv_gr_mhd_spacetime_gauge(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_mhd *gr_mhd = container_of(eqn, struct wv_gr_mhd, eqn);
  enum gkyl_spacetime_gauge spacetime_gauge = gr_mhd->spacetime_gauge;

  return spacetime_gauge;
}

int
gkyl_wv_gr_mhd_reinit_freq(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_mhd *gr_mhd = container_of(eqn, struct wv_gr_mhd, eqn);
  int reinit_freq = gr_mhd->reinit_freq;

  return reinit_freq;
}

struct gkyl_gr_spacetime*
gkyl_wv_gr_mhd_spacetime(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_mhd *gr_mhd = container_of(eqn, struct wv_gr_mhd, eqn);
  struct gkyl_gr_spacetime *spacetime = gr_mhd->spacetime;

  return spacetime;
}