
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
      sdiffs, meas, ref_ecef, b);
  fail_unless(within_epsilon(b[0], 0));
  fail_unless(within_epsilon(b[1], 0));
  fail_unless(within_epsilon(b[2], 0));

  meas[1] = 1;
  least_squares_solve_b_external_ambs(kf.state_dim, kf.state_mean,
      sdiffs, meas, ref_ecef, b);
  fail_unless(within_epsilon(b[0], -0.324757)); /* check that it matches computation made elsewhere */
  fail_unless(within_epsilon(b[1], -0.134519));
  fail_unless(within_epsilon(b[2], -0.324757));

  meas[1] = 2;
  least_squares_solve_b_external_ambs(kf.state_dim, kf.state_mean,
      sdiffs, meas, ref_ecef, b);
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
  nkf_t kf;
  kf.state_dim = 1;
  kf.obs_dim = 1;
  matrix_eye(2, kf.decor_obs_mtx);
  kf.decor_obs_cov[0] = 1;
  kf.decor_obs_cov[1] = 1;
  matrix_eye(2, kf.state_cov_U);
  kf.state_cov_D[0] = 2;
  kf.state_cov_D[1] = 3;
  kf.state_mean[0] = -1;
  kf.state_mean[1] = -2;
  double obs[2] = {1,2};
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
  /* Test that the outlier stuff says it's a good measurement
     When the dimension is 0. */
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

START_TEST(test_kf_update)
{
  nkf_t kf;
  memset(&kf, 1, sizeof(nkf_t));
  kf.state_dim = 0;
  nkf_t kf2;
  memcpy(&kf2, &kf, sizeof(nkf_t));
  double R = 0;
  double innov = 0;
  double f[10] = {0,0,0,0,0,0,0,0,0,0};
  double g[10] = {0,0,0,0,0,0,0,0,0,0};
  double alpha = 0;
  double k_scalar = 0;
  update_kf_state(&kf, R, f, g, alpha, k_scalar, innov);
  fail_unless(memcmp(&kf, &kf2, sizeof(nkf_t)) == 0);
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
  tcase_add_test(tc_core, test_kf_update);
  suite_add_tcase(s, tc_core);

  return s;
}

