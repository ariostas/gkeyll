#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_mhd.h>
#include <gkyl_wv_gr_mhd_priv.h>

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

    double magx = q[5] / sqrt(spatial_det);
    double magy = q[6] / sqrt(spatial_det);
    double magz = q[7] / sqrt(spatial_det);

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
      M_sq += (mom[i] * cov_mom[i]) / sqrt(spatial_det);
    }

    double mag[3];
    mag[0] = magx; mag[1] = magy; mag[2] = magz;

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
      tau_star += (mag[i] * cov_mom[i]) / sqrt(spatial_det);
    }

    double p_guess = 0.0;
    int iter = 0;

    while (iter < 100) {
      double a = Etot + D + p_guess + (0.5 * mag_sq);
      double d = 0.5 * ((M_sq * mag_sq) - (tau_star * tau_star));

      double phi = acos((1.0 / a) * sqrt((27.0 * d) / (4.0 * a)));
      double epsilon1 = (1.0 / 3.0) - ((2.0 / 3.0) * a * cos(((2.0 / 3.0) * phi) + ((2.0 / 3.0) * M_PI)));
      double z = epsilon1 - mag_sq;

      double v_sq = ((M_sq * (z * z)) + ((tau_star * tau_star) * (mag_sq + (2.0 * z)))) / ((z * z) * (mag_sq + z) * (mag_sq + z));
      double W = sqrt(1.0 / sqrt(1 - v_sq));
      double rho = D / W;
      double h = z / (W * W * rho);

      double p_new = rho * (h - 1.0) * ((gas_gamma - 1.0) / gas_gamma);

      if (fabs(p_guess - p_new) < pow(10.0, -8.0)) {
        iter = 100;
      }
      else {
        iter += 1;
        p_guess = p_new;
      }
    }

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
  gr_mhd->eqn.num_diag = 5;

  gr_mhd->gas_gamma = inp->gas_gamma;
  gr_mhd->spacetime_gauge = inp->spacetime_gauge;
  gr_mhd->reinit_freq = inp->reinit_freq;
  gr_mhd->spacetime = inp->spacetime;

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