#include <stdio.h>

#include <check.h>

#include "baseline.h"
#include "amb_kf.h"
#include "observation.h"
#include "linear_algebra.h"
#include "check_utils.h"

/* Need static method assign_state_rebase_mtx */
#include "amb_kf.c"

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
      sdiffs, meas, ref_ecef, b, false, DEFAULT_RAIM_THRESHOLD);
  fail_unless(within_epsilon(b[0], 0));
  fail_unless(within_epsilon(b[1], 0));
  fail_unless(within_epsilon(b[2], 0));

  meas[1] = 1;
  least_squares_solve_b_external_ambs(kf.state_dim, kf.state_mean,
      sdiffs, meas, ref_ecef, b, false, DEFAULT_RAIM_THRESHOLD);
  fail_unless(within_epsilon(b[0], -0.324757)); /* check that it matches computation made elsewhere */
  fail_unless(within_epsilon(b[1], -0.134519));
  fail_unless(within_epsilon(b[2], -0.324757));

  meas[1] = 2;
  least_squares_solve_b_external_ambs(kf.state_dim, kf.state_mean,
      sdiffs, meas, ref_ecef, b, false, DEFAULT_RAIM_THRESHOLD);
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


static void kf_test_setup(nkf_t *kf, nkf_t *kf_old, double sigma[4*4])
{
  kf->state_dim = 4;
  memcpy(kf->state_mean, (double [4]){1,2,3,4}, 4*sizeof(double));
  memcpy(sigma, (double [16]){9, 4, 8, 6,
                              4, 6, 4, 7,
                              8, 4, 8, 5,
                              6, 7, 5, 10}, 16*sizeof(double));
  double sig_test[4*4];
  memcpy(sig_test, sigma, 4*4*sizeof(double));
  matrix_udu(4, sig_test, kf->state_cov_U, kf->state_cov_D);
  /* Make sure we have a properly decomposable input matrix */
  matrix_reconstruct_udu(4, kf->state_cov_U, kf->state_cov_D, sig_test);
  fail_unless(arr_within_epsilon(16, sigma, sig_test));
  memcpy(kf_old, kf, sizeof(nkf_t));
}

static bool is_upper_unit(u8 n, double *u)
{
  for (u8 i=0; i < n; i++) {
    for (u8 j=0; j <= n; j++) {
      double uij = u[i * n + j];
      if (j < i && !within_epsilon(uij, 0)) {
        return false;
      }
      if (j == i && !within_epsilon(uij, 1)) {
        return false;
      }
      if (!isfinite(uij)) {
        return false;
      }
    }
  }
  return true;
}

static bool is_nonneg_arr(u8 n, double *a)
{
  for (u8 i=0; i < n; i++) {
    if (a[i] < 0) {
      return false;
    }
  }
  return true;
}

static bool kf_state_valid(nkf_t *kf)
{
  if (!is_upper_unit(kf->state_dim, kf->state_cov_U)) return false;
  if (!is_nonneg_arr(kf->state_dim, kf->state_cov_D)) return false;
  for (u8 i=0; i<kf->state_dim; i++) {
    if (!isfinite(kf->state_mean[i])) return false;
    if (!isfinite(kf->state_cov_D[i])) return false;
    for (u8 j=0; j<kf->state_dim; j++) {
      if (!isfinite(kf->state_cov_U[i*kf->state_dim + j])) return false;
    }
  }
  return true;
}

static bool kf_state_match(nkf_t *kf1, nkf_t *kf2)
{
  if (kf1->state_dim != kf2->state_dim) return false;
  if (!arr_within_epsilon(kf1->state_dim,
                          kf1->state_mean, kf2->state_mean)) return false;
  if (!arr_within_epsilon(kf1->state_dim,
                          kf1->state_cov_D, kf2->state_cov_D)) return false;
  if (!arr_within_epsilon(kf1->state_dim*kf1->state_dim,
                          kf1->state_cov_U, kf2->state_cov_U)) return false;
  return true;
}

static void check_inverse(u8 num_sats, const u8 *old_prns, const u8 *new_prns,
                          nkf_t *kf, nkf_t *kf_old)
{
  rebase_nkf(kf, num_sats, old_prns, new_prns);
  rebase_nkf(kf, num_sats, new_prns, old_prns);
  u8 dim = CLAMP_DIFF(num_sats, 1);
  fail_unless(kf->state_dim == dim);
  fail_unless(kf_state_match(kf, kf_old));
  fail_unless(kf_state_valid(kf));
}

START_TEST(test_rebase_nkf_inverse)
{
  nkf_t kf, kf_old;
  double sig_old[4*4];
  kf_test_setup(&kf, &kf_old, sig_old);

  /* Tests with old ref at the beginning */
  u8 old_prns[5] = {1,2,3,4,5};
  /* Rebase back and forth with new ref at beginning. */
  check_inverse(5, old_prns, (u8 [5]){2,1,3,4,5}, &kf, &kf_old);
  /* Reset and test again with new ref in middle. */
  memcpy(&kf, &kf_old, sizeof(nkf_t));
  check_inverse(5, old_prns, (u8 [5]){3,1,2,4,5}, &kf, &kf_old);
  /* Reset and test again with new ref in end. */
  memcpy(&kf, &kf_old, sizeof(nkf_t));
  check_inverse(5, old_prns, (u8 [5]){3,1,2,4,5}, &kf, &kf_old);

  /* Tests with old ref at the middle */
  memcpy(old_prns, (u8 [5]){3,1,2,4,5}, 5*sizeof(u8));
  /* Rebase back and forth with new ref at beginning. */
  memcpy(&kf, &kf_old, sizeof(nkf_t));
  check_inverse(5, old_prns, (u8 [5]){1,2,3,4,5}, &kf, &kf_old);
  /* Reset and test again with new ref in middle. */
  memcpy(&kf, &kf_old, sizeof(nkf_t));
  check_inverse(5, old_prns, (u8 [5]){2,1,3,4,5}, &kf, &kf_old);
  /* Reset and test again with new ref in end. */
  memcpy(&kf, &kf_old, sizeof(nkf_t));
  check_inverse(5, old_prns, (u8 [5]){5,1,2,3,4}, &kf, &kf_old);

  /* Tests with old ref at the end */
  memcpy(old_prns, (u8 [5]){5,1,2,3,4}, 5*sizeof(u8));
  /* Rebase back and forth with new ref at beginning. */
  memcpy(&kf, &kf_old, sizeof(nkf_t));
  check_inverse(5, old_prns, (u8 [5]){1,2,3,4,5}, &kf, &kf_old);
  /* Reset and test again with new ref in middle. */
  memcpy(&kf, &kf_old, sizeof(nkf_t));
  check_inverse(5, old_prns, (u8 [5]){2,1,3,4,5}, &kf, &kf_old);
  /* Reset and test again with new ref in end. */
  memcpy(&kf, &kf_old, sizeof(nkf_t));
  check_inverse(5, old_prns, (u8 [5]){4,1,2,3,5}, &kf, &kf_old);
}
END_TEST

static void check_identity(u8 num_sats, const u8 *old_prns, const u8 *new_prns,
                          nkf_t *kf, nkf_t *kf_old)
{
  rebase_nkf(kf, num_sats, new_prns, old_prns);
  u8 dim = CLAMP_DIFF(num_sats, 1);
  fail_unless(kf->state_dim == dim);
  fail_unless(kf_state_match(kf, kf_old));
  fail_unless(kf_state_valid(kf));
}

START_TEST(test_rebase_nkf_identity)
{
  nkf_t kf, kf_old;
  double sig_old[4 * 4];
  kf_test_setup(&kf, &kf_old, sig_old);

  /* Tests with refs at the beginning */
  u8 old_prns[5] = {1,2,3,4,5};
  u8 new_prns[5] = {1,2,3,4,5};
  check_identity(5, old_prns, new_prns, &kf, &kf_old);
  /* Reset and test again with refs in middle. */
  memcpy(&kf, &kf_old, sizeof(nkf_t));
  memcpy(old_prns, (u8 [5]){3,1,2,4,5}, 5*sizeof(u8));
  memcpy(new_prns, (u8 [5]){3,1,2,4,5}, 5*sizeof(u8));
  check_identity(5, old_prns, new_prns, &kf, &kf_old);
  /* Reset and test again with refs at end. */
  memcpy(&kf, &kf_old, sizeof(nkf_t));
  memcpy(old_prns, (u8 [5]){5,1,2,3,4}, 5*sizeof(u8));
  memcpy(new_prns, (u8 [5]){5,1,2,3,4}, 5*sizeof(u8));
  check_identity(5, old_prns, new_prns, &kf, &kf_old);
}
END_TEST

START_TEST(test_rebase_nkf)
{
  nkf_t kf, kf_old;
  double sig_old[4*4];
  kf_test_setup(&kf, &kf_old, sig_old);
  double sig_new[4*4];
  u8 old_prns[5];
  u8 new_prns[5];
  /* sig_old should now be {9, 4, 8, 6,
                            4, 6, 4, 7,
                            8, 4, 8, 5,
                            6, 7, 5, 10} */

  memcpy(old_prns, (u8 [5]){3,1,2,4,5}, 5*sizeof(u8));
  /* Rebase back and forth with new ref at beginning. */
  memcpy(new_prns, (u8 [5]){1,2,3,4,5}, 5*sizeof(u8));
  rebase_nkf(&kf, 5, new_prns, old_prns);
  matrix_reconstruct_udu(4, kf.state_cov_U, kf.state_cov_D, sig_new);
  fail_unless(arr_within_epsilon(16, sig_new,
    (double[16]) {6, 2, 2, -1,
                  2, 7, 6,  1,
                  2, 6, 6,  0,
                 -1, 1, 0,  2}));
  fail_unless(arr_within_epsilon(4, kf.state_mean,
    (double[4]) {-2, -1, 1, 2}));
  fail_unless(kf_state_valid(&kf));
  /* Reset and test again with new ref in middle. */
  memcpy(&kf, &kf_old, sizeof(nkf_t));
  memcpy(new_prns, (u8 [5]){2,1,3,4,5}, 5*sizeof(u8));
  rebase_nkf(&kf, 5, new_prns, old_prns);
  matrix_reconstruct_udu(4, kf.state_cov_U, kf.state_cov_D, sig_new);
  fail_unless(arr_within_epsilon(16, sig_new,
    (double[16]) {7,  2, 6,  1,
                  2,  6, 2, -1,
                  6,  2, 6,  0,
                  1, -1, 0,  2}));
  fail_unless(arr_within_epsilon(4, kf.state_mean,
    (double[4]) {-1, -2, 1, 2}));
  fail_unless(kf_state_valid(&kf));
  /* Reset and test again with new ref in end. */
  memcpy(&kf, &kf_old, sizeof(nkf_t));
  memcpy(new_prns, (u8 [5]){5,1,2,3,4}, 5*sizeof(u8));
  rebase_nkf(&kf, 5, new_prns, old_prns);
  matrix_reconstruct_udu(4, kf.state_cov_U, kf.state_cov_D, sig_new);
  fail_unless(arr_within_epsilon(16, sig_new,
    (double[16]) {1, 0, 1, 0,
                  0, 6, 6, 4,
                  1, 6, 8, 3,
                  0, 4, 3, 8}));
  fail_unless(arr_within_epsilon(4, kf.state_mean,
    (double[4]) {-2, -1, 1, -3}));
  fail_unless(kf_state_valid(&kf));
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


START_TEST(test_projection_identity)
{
  nkf_t kf, kf_old;
  double sig_old[4*4];
  kf_test_setup(&kf, &kf_old, sig_old);

  /* Test projecting onto itself. (shouldn't change anything) */
  nkf_state_projection(&kf, 4, 4, (u8 [4]){0,1,2,3});
  fail_unless(kf.state_dim == 4);
  fail_unless(kf_state_match(&kf, &kf_old));
}
END_TEST


START_TEST(test_projection_commutativity)
{
  nkf_t kf1, kf2, kf3;
  double sig[4*4];
  kf_test_setup(&kf1, &kf2, sig);
  memcpy(&kf3, &kf1, sizeof(nkf_t));
  fail_unless(kf_state_valid(&kf1));

  /* Project out 3 then 2 */
  nkf_state_projection(&kf1, 4, 3, (u8 [3]) {0,1,2}); /* 0,1,2,3 -> 0,1,2 */
  nkf_state_projection(&kf1, 3, 2, (u8 [2]) {0,1}); /* 0,1,2 -> 0,1 */

  /* Project out 2 then 3 */
  nkf_state_projection(&kf2, 4, 3, (u8 [3]) {0,1,3}); /* 0,1,2,3 -> 0,1,3 */
  nkf_state_projection(&kf2, 3, 2, (u8 [2]) {0,1}); /* 0,1,3 -> 0,1 */

  /* Project out 2 and 3 together */
  nkf_state_projection(&kf3, 4, 2, (u8 [2]) {0,1}); /* 0,1,2,3 -> 0,1 */

  /* Make sure they all match and are valid */
  fail_unless(kf_state_valid(&kf1));
  fail_unless(kf_state_match(&kf1, &kf2));
  fail_unless(kf_state_match(&kf2, &kf3));
}
END_TEST

START_TEST(test_projection)
{
  nkf_t kf_old, kf_new;
  double sig_old[4*4];
  kf_test_setup(&kf_old, &kf_new, sig_old);
  /* kf_new.state_mean should now be {1,2,3,4}
     sig_old should now be {9, 4, 8, 6,
                            4, 6, 4, 7,
                            8, 4, 8, 5,
                            6, 7, 5, 10} */
  fail_unless(kf_state_valid(&kf_old));
  double sig_new[3*3];
  /* Test dropping last sat. */
  nkf_state_projection(&kf_new, 4, 3, (u8 [3]){0,1,2});
  assert(kf_state_valid(&kf_new));
  matrix_reconstruct_udu(3, kf_new.state_cov_U, kf_new.state_cov_D, sig_new);
  fail_unless(arr_within_epsilon(9, sig_new, (double [3*3]) {9, 4, 8,
                                                             4, 6, 4,
                                                             8, 4, 8}));
  fail_unless(arr_within_epsilon(3, kf_new.state_mean, (double [3]) {1,2,3}));
  /* Test dropping middle sat. */
  memcpy(&kf_new, &kf_old, sizeof(nkf_t));
  nkf_state_projection(&kf_new, 4, 3, (u8 [3]) {0,2,3});
  assert(kf_state_valid(&kf_new));
  matrix_reconstruct_udu(3, kf_new.state_cov_U, kf_new.state_cov_D, sig_new);
  fail_unless(arr_within_epsilon(9, sig_new, (double [3*3]) {9, 8, 6,
                                                             8, 8, 5,
                                                             6, 5, 10}));
  fail_unless(arr_within_epsilon(3, kf_new.state_mean, (double [3]) {1,3,4}));
  /*Test dropping first sat. */
  memcpy(&kf_new, &kf_old, sizeof(nkf_t));
  nkf_state_projection(&kf_new, 4, 3, (u8 [3])  {1,2,3});
  assert(kf_state_valid(&kf_new));
  matrix_reconstruct_udu(3, kf_new.state_cov_U, kf_new.state_cov_D, sig_new);
  fail_unless(arr_within_epsilon(9, sig_new, (double [3*3]) {6, 4, 7,
                                                             4, 8, 5,
                                                             7, 5, 10}));
  fail_unless(arr_within_epsilon(3, kf_new.state_mean, (double [3]) {2,3,4}));
}
END_TEST

bool kf_inclusion_projection_match(nkf_t *kf_old, nkf_t *kf_new,
                                   u8 dim1, u8 dim2, u8 *index)
{
  double est[dim2];
  memset(est, 0, dim2*sizeof(double));
  memcpy(kf_new, kf_old, sizeof(nkf_t));
  nkf_state_inclusion(kf_new, dim1, dim2, index, est, 10);
  nkf_state_projection(kf_new, dim2, dim1, index);
  return kf_state_match(kf_old, kf_new);
}

START_TEST(test_inclusion_inversion)
{
  nkf_t kf_old, kf_new;
  double sig_old[4*4];
  kf_test_setup(&kf_old, &kf_new, sig_old);
  /* kf_new.state_mean should now be {1,2,3,4}
     sig_old should now be {9, 4, 8, 6,
                            4, 6, 4, 7,
                            8, 4, 8, 5,
                            6, 7, 5, 10} */
  fail_unless(kf_state_valid(&kf_old));
  /* Test with adding/dropping at the end. */
  fail_unless(kf_inclusion_projection_match(&kf_old, &kf_new, 4, 5, (u8 [4]){0,1,2,3}));
  /* Test with adding/dropping at the middle. */
  fail_unless(kf_inclusion_projection_match(&kf_old, &kf_new, 4, 5, (u8 [4]){0,1,3,4}));
  /* Test with adding/dropping at the beginning. */
  fail_unless(kf_inclusion_projection_match(&kf_old, &kf_new, 4, 5, (u8 [4]){1,2,3,4}));
}
END_TEST

START_TEST(test_inclusion)
{
  nkf_t kf_old, kf_new;
  double sig_old[4*4], sig_new[5*5];
  kf_test_setup(&kf_old, &kf_new, sig_old);
  /* kf_new.state_mean should now be {1,2,3,4}
     sig_old should now be {9, 4, 8, 6,
                            4, 6, 4, 7,
                            8, 4, 8, 5,
                            6, 7, 5, 10} */
  fail_unless(kf_state_valid(&kf_old));
  double est[5] = {0,0,0,0,0};
  /* Test with adding at the end. */
  nkf_state_inclusion(&kf_new, 4, 5, (u8 [4]){0,1,2,3}, est, 10);
  matrix_reconstruct_udu(5, kf_new.state_cov_U, kf_new.state_cov_D, sig_new);
  fail_unless(kf_state_valid(&kf_new));
  fail_unless(arr_within_epsilon(5*5, sig_new, (double [5*5])
    {9, 4, 8, 6,  0,
     4, 6, 4, 7,  0,
     8, 4, 8, 5,  0,
     6, 7, 5, 10, 0,
     0, 0, 0, 0,  10}));
  fail_unless(arr_within_epsilon(5, kf_new.state_mean, (double [5]){1,2,3,4,0}));
  /* Test with adding at the middle */
  memcpy(&kf_new, &kf_old, sizeof(nkf_t));
  nkf_state_inclusion(&kf_new, 4, 5, (u8 [4]){0,1,3,4}, est, 10);
  matrix_reconstruct_udu(5, kf_new.state_cov_U, kf_new.state_cov_D, sig_new);
  fail_unless(kf_state_valid(&kf_new));
  fail_unless(arr_within_epsilon(5*5, sig_new, (double [5*5])
    {9, 4, 0,  8, 6,
     4, 6, 0,  4, 7,
     0, 0, 10, 0, 0,
     8, 4, 0,  8, 5,
     6, 7, 0,  5, 10}));
  fail_unless(arr_within_epsilon(5, kf_new.state_mean, (double [5]){1,2,0,3,4}));
  /* Test with adding at the beginning */
  memcpy(&kf_new, &kf_old, sizeof(nkf_t));
  nkf_state_inclusion(&kf_new, 4, 5, (u8 [4]){1,2,3,4}, est, 10);
  matrix_reconstruct_udu(5, kf_new.state_cov_U, kf_new.state_cov_D, sig_new);
  fail_unless(kf_state_valid(&kf_new));
  fail_unless(arr_within_epsilon(5*5, sig_new, (double [5*5])
    {10, 0, 0, 0, 0,
     0,  9, 4, 8, 6,
     0,  4, 6, 4, 7,
     0,  8, 4, 8, 5,
     0,  6, 7, 5, 10}));
  fail_unless(arr_within_epsilon(5, kf_new.state_mean, (double [5]){0,1,2,3,4}));
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
  tcase_add_test(tc_core, test_rebase_nkf_inverse);
  tcase_add_test(tc_core, test_rebase_nkf_identity);
  tcase_add_test(tc_core, test_rebase_nkf);
  tcase_add_test(tc_core, test_projection_identity);
  tcase_add_test(tc_core, test_projection_commutativity);
  tcase_add_test(tc_core, test_projection);
  tcase_add_test(tc_core, test_inclusion_inversion);
  tcase_add_test(tc_core, test_inclusion);
  suite_add_tcase(s, tc_core);
  return s;
}

