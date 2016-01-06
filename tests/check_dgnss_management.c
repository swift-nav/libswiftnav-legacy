
#include <check.h>
#include <stdio.h>
#include <libswiftnav/linear_algebra.h>
#include <libswiftnav/dgnss_management.h>
#include <libswiftnav/ambiguity_test.h>
#include <libswiftnav/amb_kf.h>
#include <libswiftnav/printing_utils.h>

#include "check_utils.h"

extern sats_management_t sats_management;
extern nkf_t nkf;
extern ambiguity_test_t ambiguity_test;

START_TEST(test_dgnss_update_ambiguity_state_1)
{
  sats_management.num_sats = 5;
  sats_management.sids[0].sat = 1;
  sats_management.sids[1].sat = 2;
  sats_management.sids[2].sat = 3;
  sats_management.sids[3].sat = 4;
  sats_management.sids[4].sat = 5;
  nkf.state_dim = 4;
  nkf.state_mean[0] = 1;
  nkf.state_mean[1] = 2;
  nkf.state_mean[2] = 3;
  nkf.state_mean[3] = 4;


  ambiguity_test.amb_check.initialized = 1;
  ambiguity_test.amb_check.num_matching_ndxs = 4;
  ambiguity_test.amb_check.matching_ndxs[0] = 0;
  ambiguity_test.amb_check.matching_ndxs[1] = 2;
  ambiguity_test.amb_check.matching_ndxs[2] = 3;
  ambiguity_test.amb_check.matching_ndxs[3] = 5;
  ambiguity_test.sats.num_sats = 7;
  ambiguity_test.sats.sids[0].sat = 1;
  ambiguity_test.sats.sids[1].sat = 2;
  ambiguity_test.sats.sids[2].sat = 3;
  ambiguity_test.sats.sids[3].sat = 4;
  ambiguity_test.sats.sids[4].sat = 5;
  ambiguity_test.sats.sids[5].sat = 6;
  ambiguity_test.sats.sids[6].sat = 7;
  ambiguity_test.amb_check.ambs[0] = 20;
  ambiguity_test.amb_check.ambs[1] = 21;
  ambiguity_test.amb_check.ambs[2] = 22;
  ambiguity_test.amb_check.ambs[3] = 23;

  ambiguity_state_t s = {
    .float_ambs = {
      .n = 4,
      .sids = {{.sat = 1}, {.sat = 2}, {.sat = 3}, {.sat = 4}, {.sat = 5}},
      .ambs = {1, 2, 3, 4}
    },
    .fixed_ambs = {
      .n = 4,
      .sids = {{.sat = 1}, {.sat = 2}, {.sat = 4}, {.sat = 5}, {.sat = 7}},
      .ambs = {20, 21, 22, 23}
    }
  };

  ambiguity_state_t s_out;
  memset(&s_out, 0, sizeof(s_out));

  dgnss_update_ambiguity_state(&s_out);

  fail_unless(memcmp(&s, &s_out, sizeof(s)) == 0);
}
END_TEST

START_TEST(test_dgnss_update_ambiguity_state_2)
{
  sats_management.num_sats = 5;
  sats_management.sids[0].sat = 1;
  sats_management.sids[1].sat = 2;
  sats_management.sids[2].sat = 3;
  sats_management.sids[3].sat = 4;
  sats_management.sids[4].sat = 5;
  nkf.state_dim = 4;
  nkf.state_mean[0] = 1;
  nkf.state_mean[1] = 2;
  nkf.state_mean[2] = 3;
  nkf.state_mean[3] = 4;


  ambiguity_test.amb_check.initialized = 1;
  ambiguity_test.amb_check.num_matching_ndxs = 4;
  ambiguity_test.amb_check.matching_ndxs[0] = 0;
  ambiguity_test.amb_check.matching_ndxs[1] = 2;
  ambiguity_test.amb_check.matching_ndxs[2] = 3;
  ambiguity_test.amb_check.matching_ndxs[3] = 5;
  ambiguity_test.sats.num_sats = 7;
  ambiguity_test.sats.sids[0].sat = 1;
  ambiguity_test.sats.sids[1].sat = 2;
  ambiguity_test.sats.sids[2].sat = 3;
  ambiguity_test.sats.sids[3].sat = 4;
  ambiguity_test.sats.sids[4].sat = 5;
  ambiguity_test.sats.sids[5].sat = 6;
  ambiguity_test.sats.sids[6].sat = 7;
  ambiguity_test.amb_check.ambs[0] = 20;
  ambiguity_test.amb_check.ambs[1] = 21;
  ambiguity_test.amb_check.ambs[2] = 22;
  ambiguity_test.amb_check.ambs[3] = 23;

  ambiguity_state_t s_out;

  /* No fixed solution. */

  /* Uninitialized. */
  ambiguity_test.amb_check.initialized = 0;
  dgnss_update_ambiguity_state(&s_out);
  fail_unless(s_out.fixed_ambs.n == 0);

  /* Too few sats. */
  ambiguity_test.amb_check.initialized = 1;
  ambiguity_test.amb_check.num_matching_ndxs = 0;
  dgnss_update_ambiguity_state(&s_out);
  fail_unless(s_out.fixed_ambs.n == 0);

  ambiguity_test.amb_check.initialized = 1;
  ambiguity_test.amb_check.num_matching_ndxs = 4;

  /* No float solution. */

  /* Too few sats. */
  sats_management.num_sats = 0;
  nkf.state_dim = 0;
  dgnss_update_ambiguity_state(&s_out);
  fail_unless(s_out.float_ambs.n == 0);

  sats_management.num_sats = 1;
  nkf.state_dim = 0;
  dgnss_update_ambiguity_state(&s_out);
  fail_unless(s_out.float_ambs.n == 0);

  /* Ensure we check num_sats first as state_dim may not be valid if num_sats
   * is too low. */
  sats_management.num_sats = 1;
  nkf.state_dim = 22;
  dgnss_update_ambiguity_state(&s_out);
  fail_unless(s_out.float_ambs.n == 0);
}
END_TEST

static sdiff_t sdiffs[5];
static u8 num_sdiffs = sizeof(sdiffs) / sizeof(sdiffs[0]);
static double ref_ecef[3];

/* Initialize sdiffs used in baseline() tests. */
static void check_dgnss_baseline_setup()
{
  memset(ref_ecef, 0, sizeof(ref_ecef));

  sdiffs[0].sid.sat = 1;
  sdiffs[0].sat_pos[0] = 1;
  sdiffs[0].sat_pos[1] = 1;
  sdiffs[0].sat_pos[2] = 0;
  sdiffs[0].carrier_phase = 1;
  sdiffs[0].snr = 1.0;

  sdiffs[1].sid.sat = 2;
  sdiffs[1].sat_pos[0] = 1;
  sdiffs[1].sat_pos[1] = 0;
  sdiffs[1].sat_pos[2] = 0;
  sdiffs[1].carrier_phase = 2;
  sdiffs[1].snr = 0.0;

  sdiffs[2].sid.sat = 3;
  sdiffs[2].sat_pos[0] = 0;
  sdiffs[2].sat_pos[1] = 1;
  sdiffs[2].sat_pos[2] = 0;
  sdiffs[2].carrier_phase = 3;
  sdiffs[2].snr = 0.0;

  sdiffs[3].sid.sat = 4;
  sdiffs[3].sat_pos[0] = 0;
  sdiffs[3].sat_pos[1] = 1;
  sdiffs[3].sat_pos[2] = 1;
  sdiffs[3].carrier_phase = 4;
  sdiffs[3].snr = 0.0;

  sdiffs[4].sid.sat = 5;
  sdiffs[4].sat_pos[0] = 0;
  sdiffs[4].sat_pos[1] = 0;
  sdiffs[4].sat_pos[2] = 1;
  sdiffs[4].carrier_phase = 5;
  sdiffs[4].snr = 0.0;
}

/* No teardown required. */
static void check_dgnss_baseline_teardown(void) {}

START_TEST(test_dgnss_baseline_1)
{
  ambiguity_state_t s = {
    .float_ambs = {
      .n = 4,
      .sids = {{.sat = 1}, {.sat = 2}, {.sat = 3}, {.sat = 4}, {.sat = 5}},
      .ambs = {0, 0, 0, 0}
    },
    .fixed_ambs = {
      .n = 4,
      .sids = {{.sat = 2}, {.sat = 1}, {.sat = 3}, {.sat = 4}, {.sat = 5}},
      .ambs = {0, 1, 0, 0}
    }
  };

  double b[3];
  u8 num_used;

  /* Float only */
  s.fixed_ambs.n = 0;
  s8 valid = dgnss_baseline(num_sdiffs, sdiffs, ref_ecef, &s, &num_used, b,
    false, DEFAULT_RAIM_THRESHOLD);
  fail_unless(valid == 2);
  fail_unless(num_used == 5);
  fail_unless(within_epsilon(b[0], -0.742242));
  fail_unless(within_epsilon(b[1], -0.492905));
  fail_unless(within_epsilon(b[2], -0.0533294));

  /* Fixed and float */
  s.fixed_ambs.n = 4;
  valid = dgnss_baseline(num_sdiffs, sdiffs, ref_ecef, &s, &num_used, b,
    false, DEFAULT_RAIM_THRESHOLD);
  fail_unless(valid == 1);
  fail_unless(num_used == 5);
  fail_unless(within_epsilon(b[0], -0.417486));
  fail_unless(within_epsilon(b[1], -0.358386));
  fail_unless(within_epsilon(b[2],  0.271427));

  /* No solution possible */
  s.fixed_ambs.n = 0;
  s.float_ambs.n = 0;
  valid = dgnss_baseline(num_sdiffs, sdiffs, ref_ecef, &s, &num_used, b,
    false, DEFAULT_RAIM_THRESHOLD);
  fail_unless(valid == -1);
  fail_unless(num_used == 5);
  fail_unless(within_epsilon(b[0], -0.417486));
  fail_unless(within_epsilon(b[1], -0.358386));
  fail_unless(within_epsilon(b[2],  0.271427));
}
END_TEST

Suite* dgnss_management_test_suite(void)
{
  Suite *s = suite_create("DGNSS Management");

  TCase *tc_amb_state = tcase_create("Ambiguity State");
  tcase_add_test(tc_amb_state, test_dgnss_update_ambiguity_state_1);
  tcase_add_test(tc_amb_state, test_dgnss_update_ambiguity_state_2);
  suite_add_tcase(s, tc_amb_state);

  TCase *tc_baseline = tcase_create("Baseline");
  tcase_add_checked_fixture (tc_baseline, check_dgnss_baseline_setup,
                                          check_dgnss_baseline_teardown);
  tcase_add_test(tc_baseline, test_dgnss_baseline_1);
  suite_add_tcase(s, tc_baseline);

  return s;
}
