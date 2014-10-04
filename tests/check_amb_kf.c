
#include <check.h>
#include <stdio.h>
#include "amb_kf.h"
#include "single_diff.h"
#include "check_utils.h"

START_TEST(test_lsq) {
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

  least_squares_solve_b(&kf, sdiffs, &meas[0], ref_ecef, b);
  fail_unless(within_epsilon(b[0], 0));
  fail_unless(within_epsilon(b[1], 0));
  fail_unless(within_epsilon(b[2], 0));

  meas[1] = 1;
  least_squares_solve_b(&kf, sdiffs, &meas[0], ref_ecef, b);
  fail_unless(within_epsilon(b[0], -0.324757)); /* check that it matches computation made elsewhere */
  fail_unless(within_epsilon(b[1], -0.134519));
  fail_unless(within_epsilon(b[2], -0.324757));

  meas[1] = 2;
  least_squares_solve_b(&kf, sdiffs, &meas[0], ref_ecef, b);
  fail_unless(within_epsilon(b[0], -0.324757 * 2)); /* check that it's linear */
  fail_unless(within_epsilon(b[1], -0.134519 * 2));
  fail_unless(within_epsilon(b[2], -0.324757 * 2));
}
END_TEST

Suite* amb_kf_test_suite(void)
{
  Suite *s = suite_create("Ambiguity Kalman Filter");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_lsq);
  suite_add_tcase(s, tc_core);

  return s;
}

