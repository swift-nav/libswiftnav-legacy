
#include <check.h>
#include <stdio.h>
#include "linear_algebra.h"
#include "check_utils.h"
#include "dgnss_management.h"
#include "ambiguity_test.h"

dgnss_state_t dgs;

sdiff_t sdiffs[6];
double ref_ecef[3];

void check_dgnss_management_setup()
{
  memset(ref_ecef, 0, sizeof(double) * 3);

  sdiffs[0].prn = 1;
  sdiffs[0].sat_pos[0] = 1;
  sdiffs[0].sat_pos[1] = 1;
  sdiffs[0].sat_pos[2] = 0;
  sdiffs[0].carrier_phase = 1;

  sdiffs[1].prn = 2;
  sdiffs[1].sat_pos[0] = 1;
  sdiffs[1].sat_pos[1] = 0;
  sdiffs[1].sat_pos[2] = 0;
  sdiffs[1].carrier_phase = 2;

  sdiffs[2].prn = 3;
  sdiffs[2].sat_pos[0] = 0;
  sdiffs[2].sat_pos[1] = 1;
  sdiffs[2].sat_pos[2] = 0;
  sdiffs[2].carrier_phase = 3;

  sdiffs[3].prn = 4;
  sdiffs[3].sat_pos[0] = 0;
  sdiffs[3].sat_pos[1] = 1;
  sdiffs[3].sat_pos[2] = 1;
  sdiffs[3].carrier_phase = 4;

  sdiffs[4].prn = 5;
  sdiffs[4].sat_pos[0] = 0;
  sdiffs[4].sat_pos[1] = 0;
  sdiffs[4].sat_pos[2] = 1;
  sdiffs[4].carrier_phase = 5;

  sdiffs[5].prn = 99;

  memset(dgs.nkf.state_mean, 0, sizeof(double) * 5);
  dgs.nkf.state_dim = 4;
  dgs.nkf.obs_dim = 8;

  create_ambiguity_test(&dgs.ambiguity_test);
}

void check_dgnss_management_teardown()
{}

/** Initialise an `n` x `n` identity matrix of s32's.
 *
 * \f$ M \f$ is a matrix on \f$\mathbb{R}^{n \times n}\f$
 *
 * \param n The size of the matrix.
 * \param M Pointer to the matrix.
 */
void matrix_eye_s64(u32 n, s64 *M)
{
  /* NOTE: This function has been bounds checked. Please check again if
   * modifying. */
  memset(M, 0, n * n * sizeof(s64));
  for (u32 i=0; i<n; i++) {
    M[i*n + i] = 1;
  }
}

/* Check that it works with the first sdiff as the reference sat.
 * This should verify that the loop can start correctly.*/
START_TEST(test_dgnss_low_latency_float_baseline_ref_first) {
  dgs.sats_management.num_sats = 5;
  dgs.sats_management.prns[0] = 1;
  dgs.sats_management.prns[1] = 2;
  dgs.sats_management.prns[2] = 3;
  dgs.sats_management.prns[3] = 4;
  dgs.sats_management.prns[4] = 5;

  double b[3];
  u8 num_used;
  u8 num_sdiffs = 6;

  s8 valid = _dgnss_low_latency_float_baseline(&dgs, num_sdiffs, sdiffs,
                                 ref_ecef, &num_used, b);

  fail_unless(valid == 0);
  fail_unless(within_epsilon(b[0], -0.742242));
  fail_unless(within_epsilon(b[1], -0.492905));
  fail_unless(within_epsilon(b[2], -0.0533294));
}
END_TEST

/* Check that it works with a middle sdiff as the reference sat.
 * This should verify that the induction works. */
START_TEST(test_dgnss_low_latency_float_baseline_ref_middle) {
  dgs.sats_management.num_sats = 5;
  dgs.sats_management.prns[0] = 2;
  dgs.sats_management.prns[1] = 1;
  dgs.sats_management.prns[2] = 3;
  dgs.sats_management.prns[3] = 4;
  dgs.sats_management.prns[4] = 5;

  double b[3];
  u8 num_used;
  u8 num_sdiffs = 6;

  s8 valid = _dgnss_low_latency_float_baseline(&dgs, num_sdiffs, sdiffs,
                                 ref_ecef, &num_used, b);

  fail_unless(valid == 0);
  fail_unless(within_epsilon(b[0], -0.622609));
  fail_unless(within_epsilon(b[1], -0.432371));
  fail_unless(within_epsilon(b[2], -0.00461595));
}
END_TEST

/* Check that it works with the last sdiff as the reference sat.
 * This should verify that the loop can terminate correctly.*/
START_TEST(test_dgnss_low_latency_float_baseline_ref_end) {
  dgs.sats_management.num_sats = 5;
  dgs.sats_management.prns[0] = 5;
  dgs.sats_management.prns[1] = 1;
  dgs.sats_management.prns[2] = 2;
  dgs.sats_management.prns[3] = 3;
  dgs.sats_management.prns[4] = 4;

  double b[3];
  u8 num_used;
  u8 num_sdiffs = 5;

  s8 valid = _dgnss_low_latency_float_baseline(&dgs, num_sdiffs, sdiffs,
                                 ref_ecef, &num_used, b);

  fail_unless(valid == 0);
  fail_unless(within_epsilon(b[0], -0.589178));
  fail_unless(within_epsilon(b[1], -0.35166));
  fail_unless(within_epsilon(b[2], 0.0288157));
}
END_TEST

/* Check that measurements generated from a baseline result in an estimate
 * matching the baseline. */
START_TEST(test_dgnss_low_latency_float_baseline_fixed_point) {
  dgs.sats_management.num_sats = 5;
  dgs.sats_management.prns[0] = 5;
  dgs.sats_management.prns[1] = 1;
  dgs.sats_management.prns[2] = 2;
  dgs.sats_management.prns[3] = 3;
  dgs.sats_management.prns[4] = 4;

  double b_orig[3];
  b_orig[0] = 1;
  b_orig[1] = 1;
  b_orig[2] = 1;

  ref_ecef[0] = 0; /* Done so that we can just do the vector operations on  */
  ref_ecef[1] = 0; /*  the sat_pos vectors themselves, instead of computing */
  ref_ecef[2] = 0; /*  the line of sight vectors for each sdiff.            */

  sdiffs[0].carrier_phase = vector_dot(3, b_orig, sdiffs[0].sat_pos) /
                            vector_norm(3, sdiffs[0].sat_pos) /
                            GPS_L1_LAMBDA_NO_VAC;
  sdiffs[1].carrier_phase = vector_dot(3, b_orig, sdiffs[1].sat_pos) /
                            vector_norm(3, sdiffs[1].sat_pos) /
                            GPS_L1_LAMBDA_NO_VAC;
  sdiffs[2].carrier_phase = vector_dot(3, b_orig, sdiffs[2].sat_pos) /
                            vector_norm(3, sdiffs[2].sat_pos) /
                            GPS_L1_LAMBDA_NO_VAC;
  sdiffs[3].carrier_phase = vector_dot(3, b_orig, sdiffs[3].sat_pos) /
                            vector_norm(3, sdiffs[3].sat_pos) /
                            GPS_L1_LAMBDA_NO_VAC;
  sdiffs[4].carrier_phase = vector_dot(3, b_orig, sdiffs[4].sat_pos) /
                            vector_norm(3, sdiffs[4].sat_pos) /
                            GPS_L1_LAMBDA_NO_VAC;

  double b[3];
  u8 num_used;
  u8 num_sdiffs = 6;

  s8 valid = _dgnss_low_latency_float_baseline(&dgs, num_sdiffs, sdiffs,
                                 ref_ecef, &num_used, b);

  fail_unless(valid == 0);
  fail_unless(within_epsilon(b[0], b_orig[0]));
  fail_unless(within_epsilon(b[1], b_orig[1]));
  fail_unless(within_epsilon(b[2], b_orig[2]));
}
END_TEST

START_TEST(test_dgnss_low_latency_float_baseline_few_sats) {
  dgs.sats_management.prns[0] = 5;
  dgs.sats_management.num_sats = 1;

  double b[3];
  u8 num_used;
  u8 num_sdiffs = 6;

  s8 valid = _dgnss_low_latency_float_baseline(&dgs, num_sdiffs, sdiffs,
                                              ref_ecef, &num_used, b);

  fail_unless(valid == -1);

  dgs.sats_management.num_sats = 0;

  _dgnss_low_latency_float_baseline(&dgs, num_sdiffs, sdiffs,
                                   ref_ecef, &num_used, b);

  fail_unless(valid == -1);
}
END_TEST


// TODO(dsk) update these tests
/* Check that it works with the first sdiff as the reference sat.
 * This should verify that the loop can start correctly.*/
START_TEST(test_dgnss_low_latency_IAR_baseline_ref_first) {
  u8 prns[4];
  prns[0] = 2;
  prns[1] = 3;
  prns[2] = 4;
  prns[3] = 5;
  z_t lower[4];
  z_t upper[4];
  memset(lower, 0, sizeof(z_t) * 4);
  memset(upper, 0, sizeof(z_t) * 4);

  z_t Z_inv[16];
  matrix_eye_s64(4, Z_inv);

  add_sats_old(&dgs.ambiguity_test,
           1,
           4, prns,
           lower, upper,
           Z_inv);

  dgs.ambiguity_test.amb_check.initialized = 1;
  dgs.ambiguity_test.amb_check.num_matching_ndxs = 4;
  dgs.ambiguity_test.amb_check.matching_ndxs[0] = 0;
  dgs.ambiguity_test.amb_check.matching_ndxs[1] = 1;
  dgs.ambiguity_test.amb_check.matching_ndxs[2] = 2;
  dgs.ambiguity_test.amb_check.matching_ndxs[3] = 3;
  memset(dgs.ambiguity_test.amb_check.ambs, 0, sizeof(z_t) * 5);

  double b[3];
  u8 num_used;
  u8 num_sdiffs = 6;

  s8 valid = _dgnss_low_latency_IAR_baseline(&dgs, num_sdiffs, sdiffs,
                                 ref_ecef, &num_used, b);

  fail_unless(valid == 0);
  fail_unless(within_epsilon(b[0], -0.742242));
  fail_unless(within_epsilon(b[1], -0.492905));
  fail_unless(within_epsilon(b[2], -0.0533294));
}
END_TEST

/* Check that it works with a middle sdiff as the reference sat.
 * This should verify that the induction works. */
START_TEST(test_dgnss_low_latency_IAR_baseline_ref_middle) {
  u8 ref_prn = 2;
  u8 prns[4];
  prns[0] = 1;
  prns[1] = 3;
  prns[2] = 4;
  prns[3] = 5;
  z_t lower[4];
  z_t upper[4];
  memset(lower, 0, sizeof(z_t) * 4);
  memset(upper, 0, sizeof(z_t) * 4);

  z_t Z_inv[16];
  matrix_eye_s64(4, Z_inv);

  add_sats_old(&dgs.ambiguity_test,
           ref_prn,
           4, prns,
           lower, upper,
           Z_inv);

  dgs.ambiguity_test.amb_check.initialized = 1;
  dgs.ambiguity_test.amb_check.num_matching_ndxs = 4;
  dgs.ambiguity_test.amb_check.matching_ndxs[0] = 0;
  dgs.ambiguity_test.amb_check.matching_ndxs[1] = 1;
  dgs.ambiguity_test.amb_check.matching_ndxs[2] = 2;
  dgs.ambiguity_test.amb_check.matching_ndxs[3] = 3;
  memset(dgs.ambiguity_test.amb_check.ambs, 0, sizeof(z_t) * 5);

  double b[3];
  u8 num_used;
  u8 num_sdiffs = 6;

  s8 valid = _dgnss_low_latency_IAR_baseline(&dgs, num_sdiffs, sdiffs,
                                 ref_ecef, &num_used, b);

  fail_unless(valid == 0);
  fail_unless(within_epsilon(b[0], -0.622609));
  fail_unless(within_epsilon(b[1], -0.432371));
  fail_unless(within_epsilon(b[2], -0.00461595));
}
END_TEST

/* Check that it works with the last sdiff as the reference sat.
 * This should verify that the loop can terminate correctly.*/
START_TEST(test_dgnss_low_latency_IAR_baseline_ref_end) {
  u8 ref_prn = 5;
  u8 prns[4];
  prns[0] = 1;
  prns[1] = 2;
  prns[2] = 3;
  prns[3] = 4;
  z_t lower[4];
  z_t upper[4];
  memset(lower, 0, sizeof(z_t) * 4);
  memset(upper, 0, sizeof(z_t) * 4);

  z_t Z_inv[16];
  matrix_eye_s64(4, Z_inv);

  add_sats_old(&dgs.ambiguity_test,
           ref_prn,
           4, prns,
           lower, upper,
           Z_inv);

  dgs.ambiguity_test.amb_check.initialized = 1;
  dgs.ambiguity_test.amb_check.num_matching_ndxs = 4;
  dgs.ambiguity_test.amb_check.matching_ndxs[0] = 0;
  dgs.ambiguity_test.amb_check.matching_ndxs[1] = 1;
  dgs.ambiguity_test.amb_check.matching_ndxs[2] = 2;
  dgs.ambiguity_test.amb_check.matching_ndxs[3] = 3;
  memset(dgs.ambiguity_test.amb_check.ambs, 0, sizeof(z_t) * 5);

  double b[3];
  u8 num_used;
  u8 num_sdiffs = 5;

  s8 valid = _dgnss_low_latency_IAR_baseline(&dgs, num_sdiffs, sdiffs,
                                 ref_ecef, &num_used, b);

  fail_unless(valid == 0);
  fail_unless(within_epsilon(b[0], -0.589178));
  fail_unless(within_epsilon(b[1], -0.35166));
  fail_unless(within_epsilon(b[2], 0.0288157));
}
END_TEST

/* Check that measurements generated from a baseline result in an estimate
 * matching the baseline. */
START_TEST(test_dgnss_low_latency_IAR_baseline_fixed_point) {
  u8 ref_prn = 5;
  u8 prns[4];
  prns[0] = 1;
  prns[1] = 2;
  prns[2] = 3;
  prns[3] = 4;
  z_t lower[4];
  z_t upper[4];
  memset(lower, 0, sizeof(z_t) * 4);
  memset(upper, 0, sizeof(z_t) * 4);

  z_t Z_inv[16];
  matrix_eye_s64(4, Z_inv);

  add_sats_old(&dgs.ambiguity_test,
           ref_prn,
           4, prns,
           lower, upper,
           Z_inv);

  dgs.ambiguity_test.amb_check.initialized = 1;
  dgs.ambiguity_test.amb_check.num_matching_ndxs = 4;
  dgs.ambiguity_test.amb_check.matching_ndxs[0] = 0;
  dgs.ambiguity_test.amb_check.matching_ndxs[1] = 1;
  dgs.ambiguity_test.amb_check.matching_ndxs[2] = 2;
  dgs.ambiguity_test.amb_check.matching_ndxs[3] = 3;
  memset(dgs.ambiguity_test.amb_check.ambs, 0, sizeof(s32) * 5);

  double b_orig[3];
  b_orig[0] = 1;
  b_orig[1] = 1;
  b_orig[2] = 1;

  ref_ecef[0] = 0; /* Done so that we can just do the vector operations on  */
  ref_ecef[1] = 0; /*  the sat_pos vectors themselves, instead of computing */
  ref_ecef[2] = 0; /*  the line of sight vectors for each sdiff.            */

  sdiffs[0].carrier_phase = vector_dot(3, b_orig, sdiffs[0].sat_pos) /
                            vector_norm(3, sdiffs[0].sat_pos) /
                            GPS_L1_LAMBDA_NO_VAC;
  sdiffs[1].carrier_phase = vector_dot(3, b_orig, sdiffs[1].sat_pos) /
                            vector_norm(3, sdiffs[1].sat_pos) /
                            GPS_L1_LAMBDA_NO_VAC;
  sdiffs[2].carrier_phase = vector_dot(3, b_orig, sdiffs[2].sat_pos) /
                            vector_norm(3, sdiffs[2].sat_pos) /
                            GPS_L1_LAMBDA_NO_VAC;
  sdiffs[3].carrier_phase = vector_dot(3, b_orig, sdiffs[3].sat_pos) /
                            vector_norm(3, sdiffs[3].sat_pos) /
                            GPS_L1_LAMBDA_NO_VAC;
  sdiffs[4].carrier_phase = vector_dot(3, b_orig, sdiffs[4].sat_pos) /
                            vector_norm(3, sdiffs[4].sat_pos) /
                            GPS_L1_LAMBDA_NO_VAC;

  double b[3];
  u8 num_used;
  u8 num_sdiffs = 6;

  s8 valid = _dgnss_low_latency_IAR_baseline(&dgs, num_sdiffs, sdiffs,
                                 ref_ecef, &num_used, b);

  fail_unless(valid == 0);
  fail_unless(within_epsilon(b[0], b_orig[0]));
  fail_unless(within_epsilon(b[1], b_orig[1]));
  fail_unless(within_epsilon(b[2], b_orig[2]));
}
END_TEST

START_TEST(test_dgnss_low_latency_IAR_baseline_few_sats) {
  dgs.ambiguity_test.amb_check.initialized = 1;
  dgs.ambiguity_test.amb_check.num_matching_ndxs = 1;

  double b[3];
  u8 num_used;
  u8 num_sdiffs = 6;


  s8 valid = _dgnss_low_latency_IAR_baseline(&dgs, num_sdiffs, sdiffs,
                                              ref_ecef, &num_used, b);
  fail_unless(valid == -1);

  dgs.ambiguity_test.amb_check.num_matching_ndxs = 0;

  _dgnss_low_latency_float_baseline(&dgs, num_sdiffs, sdiffs,
                                   ref_ecef, &num_used, b);
  fail_unless(valid == -1);

  dgs.ambiguity_test.amb_check.initialized = 0;
  dgs.ambiguity_test.amb_check.num_matching_ndxs = 4;

  _dgnss_low_latency_float_baseline(&dgs, num_sdiffs, sdiffs,
                                   ref_ecef, &num_used, b);
  fail_unless(valid == -1);
}
END_TEST

START_TEST(test_dgnss_low_latency_IAR_baseline_uninitialized) {
  dgs.ambiguity_test.amb_check.initialized = 0;
  dgs.ambiguity_test.amb_check.num_matching_ndxs = 5;

  double b[3];
  u8 num_used;
  u8 num_sdiffs = 6;


  s8 valid = _dgnss_low_latency_IAR_baseline(&dgs, num_sdiffs, sdiffs,
                                              ref_ecef, &num_used, b);
  fail_unless(valid == -1);
}
END_TEST

START_TEST(test_dgnss_low_latency_baseline_uninitialized) {
  dgs.ambiguity_test.amb_check.initialized = 0;
  dgs.ambiguity_test.amb_check.num_matching_ndxs = 5;

  double b[3];
  u8 num_used;
  u8 num_sdiffs = 6;


  s8 valid = _dgnss_low_latency_IAR_baseline(&dgs, num_sdiffs, sdiffs,
                                              ref_ecef, &num_used, b);
  fail_unless(valid == -1);
}
END_TEST

Suite* dgnss_management_test_suite(void)
{
  Suite *s = suite_create("DGNSS Management");

  TCase *tc_core = tcase_create("Core");
  tcase_add_checked_fixture (tc_core, check_dgnss_management_setup,
                                      check_dgnss_management_teardown);
  tcase_add_test(tc_core, test_dgnss_low_latency_float_baseline_ref_first);
  tcase_add_test(tc_core, test_dgnss_low_latency_float_baseline_ref_middle);
  tcase_add_test(tc_core, test_dgnss_low_latency_float_baseline_ref_end);
  tcase_add_test(tc_core, test_dgnss_low_latency_float_baseline_fixed_point);
  tcase_add_test(tc_core, test_dgnss_low_latency_float_baseline_few_sats);

  tcase_add_test(tc_core, test_dgnss_low_latency_IAR_baseline_ref_first);
  tcase_add_test(tc_core, test_dgnss_low_latency_IAR_baseline_ref_middle);
  tcase_add_test(tc_core, test_dgnss_low_latency_IAR_baseline_ref_end);
  tcase_add_test(tc_core, test_dgnss_low_latency_IAR_baseline_fixed_point);
  tcase_add_test(tc_core, test_dgnss_low_latency_IAR_baseline_few_sats);
  tcase_add_test(tc_core, test_dgnss_low_latency_IAR_baseline_uninitialized);
  tcase_add_test(tc_core, test_dgnss_low_latency_baseline_uninitialized);
  suite_add_tcase(s, tc_core);

  return s;
}
