#include <stdio.h>

#include <check.h>

#include "baseline.h"
#include "amb_kf.h"
#include "observation.h"
#include "linear_algebra.h"
#include "check_utils.h"

START_TEST(test_lsq)
{
  sdiff_t sdiffs[5];
  sdiffs[0].sat_pos[0] = 1;
  sdiffs[0].sat_pos[1] = 0;
  sdiffs[0].sat_pos[2] = 0;

  sdiffs[1].sat_pos[0] = 1;
  sdiffs[1].sat_pos[1] = 1;
  sdiffs[1].sat_pos[2] = 0;

  sdiffs[2].sat_pos[0] = 0;
  sdiffs[2].sat_pos[1] = 1;
  sdiffs[2].sat_pos[2] = 0;

  sdiffs[3].sat_pos[0] = 0;
  sdiffs[3].sat_pos[1] = 0;
  sdiffs[3].sat_pos[2] = 1;

  sdiffs[4].sat_pos[0] = 0;
  sdiffs[4].sat_pos[1] = 1;
  sdiffs[4].sat_pos[2] = 1;

  nkf_t kf;
  kf.state_mean[0] = 0;
  kf.state_mean[1] = 0;
  kf.state_mean[2] = 0;
  kf.state_mean[3] = 0;
  kf.state_dim = 4;

  double meas[4];
  meas[0] = 0;
  meas[1] = 0;
  meas[2] = 0;
  meas[3] = 0;

  double b[3];

  double ref_ecef[3];
  ref_ecef[0] = 0;
  ref_ecef[1] = 0;
  ref_ecef[2] = 0;

  least_squares_solve_b_external_ambs(kf.state_dim, kf.state_mean,
      sdiffs, meas, ref_ecef, b, false);
  fail_unless(within_epsilon(b[0], 0));
  fail_unless(within_epsilon(b[1], 0));
  fail_unless(within_epsilon(b[2], 0));

  meas[1] = 1;
  least_squares_solve_b_external_ambs(kf.state_dim, kf.state_mean,
      sdiffs, meas, ref_ecef, b, false);
  fail_unless(within_epsilon(b[0], -0.324757)); /* check that it matches computation made elsewhere */
  fail_unless(within_epsilon(b[1], -0.134519));
  fail_unless(within_epsilon(b[2], -0.324757));

  meas[1] = 2;
  least_squares_solve_b_external_ambs(kf.state_dim, kf.state_mean,
      sdiffs, meas, ref_ecef, b, false);
  fail_unless(within_epsilon(b[0], -0.324757 * 2)); /* check that it's linear */
  fail_unless(within_epsilon(b[1], -0.134519 * 2));
  fail_unless(within_epsilon(b[2], -0.324757 * 2));
}
END_TEST

START_TEST(test_sos_innov)
{
  /* Test with H, U, D, R = I */
  nkf_t kf;
  kf.state_dim = 2;
  kf.obs_dim = 2;
  matrix_eye(2, kf.decor_obs_mtx);
  kf.decor_obs_cov[0] = 1;
  kf.decor_obs_cov[1] = 1;
  matrix_eye(2, kf.state_cov_U);
  kf.state_cov_D[0] = 2;
  kf.state_cov_D[1] = 3;
  kf.state_mean[0] = 1;
  kf.state_mean[1] = 2;
  double obs[2] = {1,2};
  /* Test with the predicted and actual observations the same. */
  fail_unless(within_epsilon(get_sos_innov(&kf, obs), 0));
  fail_unless(get_sos_innov(&kf, obs) >= 0);
  kf.state_mean[0] = 0;
  kf.state_mean[1] = 0;
  /* Test with the predicted and actual observations slightly different.
   * S = R + H * U * D * U^T * H^T = diag(3,4).
   * y = z - H*x = z = (1, 2).
   * sos = sum_i (y_i / S_ii) = 1/3 + 2^2 / 4. */
  fail_unless(within_epsilon(get_sos_innov(&kf, obs), 1.0f/3 + 4.0f/4));
  /* Test it with a singular matrix.
     kf.state_cov = {{1,1},{1,1}} */
  matrix_eye(2, kf.state_cov_U);
  kf.state_cov_U[1] = 1;
  kf.state_cov_D[0] = 0;
  kf.state_cov_D[1] = 1;
  memset(kf.decor_obs_cov, 0, 2*sizeof(double));
  fail_unless(isfinite(get_sos_innov(&kf, obs)));
  fail_unless(get_sos_innov(&kf, obs) >= 0);
}
END_TEST

START_TEST(test_sos_innov_dims)
{
  /* Make sure that the SOS innovation calculation works for a few different
   * combinations of state_dim and obs_dim, including zeros.*/
  nkf_t kf;
  kf.decor_obs_cov[0] = 1;
  kf.decor_obs_cov[1] = 1;
  matrix_eye(2, kf.state_cov_U);
  kf.state_cov_D[0] = 2;
  kf.state_cov_D[1] = 3;
  kf.state_mean[0] = -1;
  kf.state_mean[1] = -2;
  double obs[2] = {1,2};
  kf.state_dim = 1;
  kf.obs_dim = 1;
  matrix_eye(2, kf.decor_obs_mtx);
  fail_unless(within_epsilon(get_sos_innov(&kf, obs), 4.0f/3));
  kf.state_dim = 0;
  kf.obs_dim = 1;
  fail_unless(get_sos_innov(&kf, obs) == 0);
  kf.state_dim = 1;
  kf.obs_dim = 0;
  fail_unless(get_sos_innov(&kf, obs) == 0);
  kf.state_dim = 0;
  kf.obs_dim = 0;
  fail_unless(get_sos_innov(&kf, obs) == 0);
  kf.state_dim = 1;
  kf.obs_dim = 2;
  kf.decor_obs_mtx[0] = 1;
  kf.decor_obs_mtx[1] = 1;
  fail_unless(within_epsilon(get_sos_innov(&kf, obs), 4.0f/3 + 9.0f/3));
  kf.state_dim = 2;
  kf.obs_dim = 1;
  fail_unless(within_epsilon(get_sos_innov(&kf, obs), 16.0f/6));
}
END_TEST

START_TEST(test_outlier)
{
  /* Test with H, U, D, R = I */
  nkf_t kf;
  kf.state_dim = 2;
  kf.obs_dim = 2;
  matrix_eye(2, kf.decor_obs_mtx);
  kf.decor_obs_cov[0] = 1;
  kf.decor_obs_cov[1] = 1;
  matrix_eye(2, kf.state_cov_U);
  kf.state_cov_D[0] = 2;
  kf.state_cov_D[1] = 3;
  kf.state_mean[0] = 1;
  kf.state_mean[1] = 2;
  double obs[2] = {1,2};
  double k_scalar;
  /* Test with the predicted and actual observations the same. */
  bool is_outlier = outlier_check(&kf, obs, &k_scalar);
  /* A perfect match should never be an outlier, so k_scalar must be 1. */
  fail_unless(!is_outlier);
  fail_unless(within_epsilon(k_scalar, 1));
  kf.state_mean[0] = 1e20;
  kf.state_mean[1] = 1e20;

  /* Test with the predicted and actual observations wildly different. */
  is_outlier = outlier_check(&kf, obs, &k_scalar);
  /* This horrible of a match should always be an outlier,
   * so k_scalar must be < 1. */
  fail_unless(is_outlier);
  fail_unless(k_scalar < 1);
}
END_TEST

START_TEST(test_outlier_dims)
{
  /* Test that the outlier detection says that a null measurement
   * (either dim = 0) is a good measurement. Tests different combos
   * of which dim is zero. */
  nkf_t kf = {.state_dim = 1,
              .obs_dim = 0};
  double obs;
  double k_scalar;
  bool bad = outlier_check(&kf, &obs, &k_scalar);
  fail_unless(bad == false);
  fail_unless(within_epsilon(k_scalar, 1));
  kf.state_dim = 0;
  kf.obs_dim = 1;
  bad = outlier_check(&kf, &obs, &k_scalar);
  fail_unless(bad == false);
  fail_unless(within_epsilon(k_scalar, 1));
  kf.state_dim = 0;
  kf.obs_dim = 0;
  bad = outlier_check(&kf, &obs, &k_scalar);
  fail_unless(bad == false);
  fail_unless(within_epsilon(k_scalar, 1));
}
END_TEST

START_TEST(test_kf_update_noop)
{
  /* Make sure that k_scalar = 0 results in a noop in the KF measurement
   * update. */
  nkf_t kf;
  memset(&kf, 1, sizeof(nkf_t));
  nkf_t kf2;
  memcpy(&kf2, &kf, sizeof(nkf_t));
  double R = 0;
  double innov = 0;
  double f[1] = {0};
  double g[1] = {0};
  double alpha = 0;
  double k_scalar = 0;
  update_kf_state(&kf, R, f, g, alpha, k_scalar, innov);
  fail_unless(memcmp(&kf, &kf2, sizeof(nkf_t)) == 0);
}
END_TEST

START_TEST(test_kf_update)
{
  /* Test that random full rank KFs coded the slow, but naive way match up with
   * those done the factored way. */
  u8 dim = 10;
  for (u32 i=0; i < 1000; i++) {
    /* All the allocations up here because they get in the way. */
    nkf_t kf = {.state_dim = dim};
    double f[dim];
    double g[dim];
    double m[dim * dim];
    double mt[dim * dim];
    double p[dim * dim];
    double p2[dim * dim];
    double h[dim];
    double ph[dim];
    double state_mean[dim];
    double k[dim];
    double x[dim];
    double kh[dim * dim];
    double eye[dim * dim];
    double sc[dim * dim];
    double p3[dim * dim];
    /* A bunch of random values to try. The only special random array is P, which
     * is also (almost surely) full rank pos def (M*M' for random M). */
    arr_frand(dim, -1, 1, state_mean);
    double R = frand(1e-12, 1);
    double k_scalar = frand(0, 1);
    double innov = frand(-1,1);
    arr_frand(dim, -1, 1, h);
    arr_frand(dim * dim, -1, 1, m);
    matrix_transpose(dim, dim, m, mt);
    matrix_multiply(dim, dim, dim, m, mt, p);
    matrix_transpose(dim, dim, p, p2);
    fail_unless(arr_within_epsilon(dim * dim, p, p2));
    matrix_udu(dim, p, kf.state_cov_U, kf.state_cov_D);
    matrix_transpose(dim, dim, p2, p);
    memcpy(kf.state_mean, state_mean, dim * sizeof(double));
    /* Compute the factored KF update */
    double alpha = compute_innovation_terms(dim, h, R,
                                     kf.state_cov_U, kf.state_cov_D,
                                     f, g);
    update_kf_state(&kf, R, f, g,
                     alpha, k_scalar,
                     innov);
    matrix_reconstruct_udu(dim, kf.state_cov_U, kf.state_cov_D, p2);
    /* Compute the simple KF update: */
    /* S = H * P * H' + R; */
    matrix_multiply(dim, dim, 1, p, h, ph);
    double s = vector_dot(dim, h, ph) + R;
    /* Make sure the innovation variances match */
    fail_unless(within_epsilon(s, alpha));
    /* K = P* H' * k_scalar * S^-1 */
    matrix_multiply(dim, dim, 1, p, h, ph);
    for (u8 j=0; j < dim; j++) {
      k[j] = ph[j] * k_scalar / s;
    }
    /* x_next = x_prev + k * innov */
    vector_add_sc(dim, state_mean, k, innov, x);
    /* p_next = (I - K * H) * p_prev */
    matrix_multiply(dim, 1, dim, k, h, kh);
    matrix_eye(dim, eye);
    matrix_add_sc(dim, dim, eye, kh, -1, sc);
    matrix_multiply(dim, dim, dim, sc, p, p3);
    /* Make sure the mean and covariances match */
    fail_unless(arr_within_epsilon(dim, x, kf.state_mean));
    fail_unless(arr_within_epsilon(dim * dim, p2, p3));
  }
}
END_TEST

void assign_state_rebase_mtx(const u8 num_sats, const u8 *old_prns,
                             const u8 *new_prns, double *rebase_mtx);

START_TEST(test_rebase_state)
{
  int num_sats = 6;
  int dim = num_sats - 1;

  double m1[dim*dim];
  double m2[dim*dim];

  double m3[dim*dim];
  double m4[dim*dim];
  double id[dim*dim];

  u8 prns1[] = {2,1,3,4,5,6};
  u8 prns2[] = {5,1,2,3,4,6};

  assign_state_rebase_mtx(num_sats, prns1, prns2, m1);
  assign_state_rebase_mtx(num_sats, prns2, prns1, m2);

  matrix_multiply(dim, dim, dim, m1, m2, m3);
  matrix_multiply(dim, dim, dim, m2, m1, m4);

  matrix_eye(dim, id);

  fail_unless(arr_within_epsilon(dim*dim, m3, id));
  fail_unless(arr_within_epsilon(dim*dim, m4, id));
}
END_TEST

Suite* amb_kf_test_suite(void)
{
  Suite *s = suite_create("Ambiguity Kalman Filter");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_lsq);
  tcase_add_test(tc_core, test_sos_innov);
  tcase_add_test(tc_core, test_outlier);
  tcase_add_test(tc_core, test_sos_innov_dims);
  tcase_add_test(tc_core, test_outlier_dims);
  tcase_add_test(tc_core, test_kf_update_noop);
  tcase_add_test(tc_core, test_kf_update);
  tcase_add_test(tc_core, test_rebase_state);
  suite_add_tcase(s, tc_core);

  return s;
}

