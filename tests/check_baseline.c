#include <check.h>
#include <math.h>
#include <stdio.h>

#include <linear_algebra.h>
#include <constants.h>
#include <baseline.h>

#include "check_utils.h"
#include "printing_utils.h"

START_TEST(test_predict_carrier_obs)
{
  double N[] = {22.2, 23.3, 34.4, 123.4};
  u8 num_dds = sizeof(N)/sizeof(N[0]);

  double DE[] = {1, 0, 0,
                 0, 1, 0,
                 0, 0, 1,
                 1, 1, 1};
  double b[3] = {1, 1, 1};
  double dd_obs_expected[] = {
    1.0/0.19023800915688557 + 22.2,
    1.0/0.19023800915688557 + 23.3,
    1.0/0.19023800915688557 + 34.4,
    3.0/0.19023800915688557 + 123.4
  };
  double dd_obs[num_dds];

  predict_carrier_obs(num_dds, N, DE, b, dd_obs);

  for (u8 i=0; i<num_dds; i++) {
    fail_unless(within_epsilon(dd_obs[i], dd_obs_expected[i]),
                "Observation mismatch: %lf vs %lf",
                dd_obs[i], dd_obs_expected[i]);
  }
}
END_TEST

START_TEST(test_predict_carrier_obs2)
{
  double N[1] = {222.2};
  double DE[3] = {3, 4, 5};
  double b[3] = {7, 8, 9};
  double dd_obs_expected[1] = {
    (3*7 + 4*8 + 5*9)/0.19023800915688557 + 222.2,
  };
  double dd_obs[1];

  predict_carrier_obs(1, N, DE, b, dd_obs);

  fail_unless(within_epsilon(dd_obs[0], dd_obs_expected[0]),
              "Observation mismatch: %lf vs %lf",
              dd_obs[0], dd_obs_expected[0]);
}
END_TEST

START_TEST(test_amb_from_baseline)
{
  double N_true[] = {22, 23, 34, -123};
  u8 num_dds = sizeof(N_true)/sizeof(N_true[0]);
  s32 N[num_dds];

  double DE[] = {1, 0, 0,
                 0, 1, 0,
                 0, 0, 1,
                 1, 1, 1};
  double b[3] = {1, 1, 1};
  double dd_obs[num_dds];

  predict_carrier_obs(num_dds, N_true, DE, b, dd_obs);
  amb_from_baseline(num_dds, DE, dd_obs, b, N);

  for (u8 i=0; i<num_dds; i++) {
    fail_unless(within_epsilon(N_true[i], N[i]),
                "Ambiguity mismatch: %ld vs %lf",
                N[i], N_true[i]);
  }
}
END_TEST

START_TEST(test_lesq_solution)
{
  /* Over constrained. */
  double N[] = {22, 23, 34, -123};
  u8 num_dds = sizeof(N)/sizeof(N[0]);

  double DE[] = {1, 0, 0,
                 0, 1, 0,
                 0, 0, 1,
                 1, 1, 1};
  double b_true[3] = {1.234, 1.456, 1.789};
  double dd_obs[num_dds];

  predict_carrier_obs(num_dds, N, DE, b_true, dd_obs);

  double b[3];
  double resid[num_dds];
  s8 ret = lesq_solution_float(num_dds, dd_obs, N, DE, b, resid);

  fail_unless(ret == 0, "solution returned error %d", ret);

  for (u8 i=0; i<3; i++) {
    fail_unless(within_epsilon(b[i], b_true[i]),
                "Baseline mismatch: %lf vs %lf",
                b[i], b_true[i]);
  }

  /* Try with resid = NULL */
  ret = lesq_solution_float(num_dds, dd_obs, N, DE, b, 0);

  fail_unless(ret == 0, "solution returned error %d", ret);

  for (u8 i=0; i<3; i++) {
    fail_unless(within_epsilon(b[i], b_true[i]),
                "Baseline mismatch: %lf vs %lf",
                b[i], b_true[i]);
  }
}
END_TEST

START_TEST(test_lesq_solution2)
{
  /* Exactly constrained. */
  double N[] = {22, 23, 34};
  u8 num_dds = sizeof(N)/sizeof(N[0]);

  double DE[] = {1, 0, 0,
                 0, 1, 0,
                 1, 1, 1};
  double b_true[3] = {1.234, 1.456, 1.789};
  double dd_obs[num_dds];

  predict_carrier_obs(num_dds, N, DE, b_true, dd_obs);

  double b[3];
  double resid[num_dds];
  s8 ret = lesq_solution_float(num_dds, dd_obs, N, DE, b, resid);

  fail_unless(ret == 0, "solution returned error %d", ret);

  for (u8 i=0; i<3; i++) {
    fail_unless(within_epsilon(b[i], b_true[i]),
                "Baseline mismatch: %lf vs %lf",
                b[i], b_true[i]);
  }
}
END_TEST

START_TEST(test_lesq_solution3)
{
  /* Under constrained, should fail with correct return code. */
  double N[2];
  double DE[2*3];
  double b[3];
  double resid[2];
  double dd_obs[2];

  for (u8 num_dds = 0; num_dds < 3; num_dds++) {
    s8 ret = lesq_solution_float(num_dds, dd_obs, N, DE, b, resid);
    fail_unless(ret == -1, "solution under-constrained, "
                "should have returned error -1, got %d, dds = %d",
                ret, num_dds);
  }
}
END_TEST

START_TEST(test_lesq_solution4)
{
  /* Over constrained, integer valued ambiguity */
  double N[] = {22, 23, 34, -123};
  u8 num_dds = sizeof(N)/sizeof(N[0]);

  double DE[] = {1, 0, 0,
                 0, 1, 0,
                 0, 0, 1,
                 1, 1, 1};
  double b_true[3] = {1.234, 1.456, 1.789};
  double dd_obs[num_dds];

  predict_carrier_obs(num_dds, N, DE, b_true, dd_obs);

  double b[3];
  double resid[num_dds];
  double N_int[num_dds];

  for (u8 i=0; i<num_dds; i++) {
    N_int[i] = N[i];
  }
  s8 ret = lesq_solve_raim(num_dds, dd_obs, N_int, DE, b,
    false, DEFAULT_RAIM_THRESHOLD, 0, resid, 0);

  fail_unless(ret >= 0, "solution returned error %d", ret);

  for (u8 i=0; i<3; i++) {
    fail_unless(within_epsilon(b[i], b_true[i]),
                "Baseline mismatch: %lf vs %lf",
                b[i], b_true[i]);
  }

  /* Try with null resid */
  ret = lesq_solve_raim(num_dds, dd_obs, N_int, DE, b,
    false, DEFAULT_RAIM_THRESHOLD, 0, 0, 0);

  fail_unless(ret >= 0, "solution returned error %d", ret);

  for (u8 i=0; i<3; i++) {
    fail_unless(within_epsilon(b[i], b_true[i]),
                "Baseline mismatch: %lf vs %lf",
                b[i], b_true[i]);
  }
}
END_TEST

START_TEST(test_lesq_solution5)
{
  /* Over constrained with non-zero residuals */
  double N[] = {0, 0, 0, 0};
  u8 num_dds = sizeof(N)/sizeof(N[0]);

  double DE[] = {1, 0, 0,
                 0, 1, 0,
                 0, 0, 1,
                 1, 0, 0};
  double b_true[3] = {0, 0, 0};
  double dd_obs[] = {0, 0, 0, 1};

  double b[3];
  double resid[num_dds];
  s8 ret = lesq_solution_float(num_dds, dd_obs, N, DE, b, resid);

  double b_expected[3] = {0.5*GPS_L1_LAMBDA_NO_VAC, 0, 0};
  double resid_expected[] = {-0.5, 0, 0, 0.5};

  fail_unless(ret == 0, "solution returned error %d", ret);

  for (u8 i=0; i<3; i++) {
    fail_unless(within_epsilon(b[i], b_expected[i]),
                "Baseline mismatch: %lf vs %lf",
                b[i], b_true[i]);
  }
  for (u8 i=0; i<num_dds; i++) {
    fail_unless(within_epsilon(resid[i], resid_expected[i]),
                "Residual mismatch: %lf vs %lf",
                resid[i], resid_expected[i]);
  }
}
END_TEST

static sdiff_t sdiffs[5];
static u8 num_sdiffs = sizeof(sdiffs) / sizeof(sdiffs[0]);
static double ref_ecef[3];

/* Initialize sdiffs used in baseline() tests. */
static void check_baseline_setup()
{
  memset(ref_ecef, 0, sizeof(ref_ecef));

  sdiffs[0].sid.sat = 1;
  sdiffs[0].sat_pos[0] = 1;
  sdiffs[0].sat_pos[1] = 1;
  sdiffs[0].sat_pos[2] = 0;
  sdiffs[0].carrier_phase = 1;
  sdiffs[0].snr = 0.0;

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
static void check_baseline_teardown(void) {}

START_TEST(test_baseline_ref_first)
{
  /* Check that it works with the first sdiff as the reference sat.
   * This should verify that the loop can start correctly. */
  ambiguities_t ambs = {
    .n = 4,
    .prns = {{.sat = 1}, {.sat = 2}, {.sat = 3}, {.sat = 4}, {.sat = 5}},
    .ambs = {0, 0, 0, 0}
  };

  double b[3];
  u8 num_used;

  s8 valid = baseline(num_sdiffs, sdiffs, ref_ecef, &ambs, &num_used, b,
    false, DEFAULT_RAIM_THRESHOLD);

  fail_unless(valid == 0);
  fail_unless(num_used == 5);
  fail_unless(within_epsilon(b[0], -0.742242));
  fail_unless(within_epsilon(b[1], -0.492905));
  fail_unless(within_epsilon(b[2], -0.0533294));
}
END_TEST

START_TEST(test_baseline_ref_middle)
{
  /* Check that it works with a middle sdiff as the reference sat.
   * This should verify that the induction works. */
  ambiguities_t ambs = {
    .n = 4,
    .prns = {{.sat = 1}, {.sat = 2}, {.sat = 3}, {.sat = 4}, {.sat = 5}},
    .ambs = {0, 0, 0, 0}
  };

  double old_snr_2 = sdiffs[1].snr;
  sdiffs[1].snr = 22;

  double b[3];
  u8 num_used;

  s8 valid = baseline(num_sdiffs, sdiffs, ref_ecef, &ambs, &num_used, b,
    false, DEFAULT_RAIM_THRESHOLD);

  fail_unless(valid == 0);
  fail_unless(num_used == 5);
  fail_unless(within_epsilon(b[0], -0.622609));
  fail_unless(within_epsilon(b[1], -0.432371));
  fail_unless(within_epsilon(b[2], -0.00461595));

  sdiffs[1].snr = old_snr_2;
}
END_TEST

START_TEST(test_baseline_ref_end)
{
  /* Check that it works with the last sdiff as the reference sat.
   * This should verify that the loop can terminate correctly.*/
  ambiguities_t ambs = {
    .n = 4,
    .prns = {{.sat = 1}, {.sat = 2}, {.sat = 3}, {.sat = 4}, {.sat = 5}},
    .ambs = {0, 0, 0, 0}
  };

  u8 index = ambs.n;
  double old_snr_2 = sdiffs[index].snr;
  sdiffs[index].snr = 22;

  double b[3];
  u8 num_used;

  s8 valid = baseline(num_sdiffs, sdiffs, ref_ecef, &ambs, &num_used, b,
    false, DEFAULT_RAIM_THRESHOLD);

  fail_unless(valid == 0);
  fail_unless(num_used == 5);
  fail_unless(within_epsilon(b[0], -0.589178));
  fail_unless(within_epsilon(b[1], -0.35166));
  fail_unless(within_epsilon(b[2], 0.0288157));

  sdiffs[index].snr = old_snr_2;
}
END_TEST

START_TEST(test_baseline_fixed_point)
{
  /* Check that measurements generated from a baseline result in an estimate
   * matching the baseline. */
  ambiguities_t ambs = {
    .n = 4,
    .prns = {{.sat = 5}, {.sat = 1}, {.sat = 2}, {.sat = 3}, {.sat = 4}},
    .ambs = {0, 0, 0, 0}
  };

  double b_orig[3] = {1, 1, 1};

  ref_ecef[0] = 0; /* Done so that we can just do the vector operations on  */
  ref_ecef[1] = 0; /*  the sat_pos vectors themselves, instead of computing */
  ref_ecef[2] = 0; /*  the line of sight vectors for each sdiff.            */

  for (u8 i=0; i<5; i++) {
    sdiffs[i].carrier_phase = vector_dot(3, b_orig, sdiffs[i].sat_pos) /
                              vector_norm(3, sdiffs[i].sat_pos) /
                              GPS_L1_LAMBDA_NO_VAC;
  }

  double b[3];
  u8 num_used;

  s8 valid = baseline(num_sdiffs, sdiffs, ref_ecef, &ambs, &num_used, b,
    false, DEFAULT_RAIM_THRESHOLD);

  fail_unless(valid == 0);
  fail_unless(num_used == 5);
  fail_unless(within_epsilon(b[0], b_orig[0]));
  fail_unless(within_epsilon(b[1], b_orig[1]));
  fail_unless(within_epsilon(b[2], b_orig[2]));
}
END_TEST

START_TEST(test_baseline_few_sats)
{
  ambiguities_t ambs = {
    .n = 0,
    .prns = {{.sat = 5}, {.sat = 1}, {.sat = 2}, {.sat = 3}, {.sat = 4}},
    .ambs = {0, 0, 0, 0}
  };

  double b[3];
  u8 num_used;

  s8 valid = baseline(num_sdiffs, sdiffs, ref_ecef, &ambs, &num_used, b,
    false, DEFAULT_RAIM_THRESHOLD);

  fail_unless(valid == -1);
}
END_TEST

/* Test raim repair */
START_TEST(test_lesq_repair8)
{
  /* Over constrained with bad DE row. */
  double N[] = {0,0,0,0,0,0,0,0};
  u8 num_dds = sizeof(N)/sizeof(N[0]);

  double DE[] = {1, 0, 0,
                 0, 1, 0,
                 0, 0, 1,
                 1, 1, 1,
                 1, 1, 1,
                 1, 1, 1,
                 1, 1, 1,
                 22, 222, 2222};
  double dd_obs[] = {1, 1, 1, 3, 3,3,3,3};

  double b[3];
  u8 bad_index;
  s8 ret = lesq_solve_raim(num_dds, dd_obs, N, DE, b,
    false, DEFAULT_RAIM_THRESHOLD, 0, 0, &bad_index);

  fail_unless(ret == 1,
    "Expecting 1 for repaired solution, got: %i.\n", ret);
  fail_unless(bad_index == 7,
    "Expecting repaired solution (dropping index 4 of DE), got: %i.\n",
    bad_index);
}
END_TEST

/* Test raim repair */
START_TEST(test_lesq_repair1)
{
  /* Over constrained with bad DE row. */
  double N[] = {0, 0, 0, 0, 0};
  u8 num_dds = sizeof(N)/sizeof(N[0]);

  double DE[] = {1, 0, 0,
                 0, 1, 0,
                 0, 0, 1,
                 1, 1, 1,
                 22, 222, 2222};
  double dd_obs[] = {1, 1, 1, 3, 1};

  double b[3];
  u8 bad_index;
  s8 ret = lesq_solve_raim(num_dds, dd_obs, N, DE, b,
    false, DEFAULT_RAIM_THRESHOLD, 0, 0, &bad_index);

  fail_unless(ret == 1,
    "Expecting 1 for repaired solution, got: %i.\n", ret);
  fail_unless(bad_index == 4,
    "Expecting repaired solution (dropping index 4 of DE), got: %i.\n", bad_index);
}
END_TEST

/* Test raim disabling flag */
START_TEST(test_lesq_repair_disabled)
{
  /* Over constrained with bad DE row. */
  double N[] = {0, 0, 0, 0, 0};
  u8 num_dds = sizeof(N)/sizeof(N[0]);

  double DE[] = {1, 0, 0,
                 0, 1, 0,
                 0, 0, 1,
                 1, 1, 1,
                 22, 222, 2222};
  double dd_obs[] = {1, 1, 1, 3, 1};

  double b[3];
  u8 bad_index;
  /* DISABLE raim */
  s8 ret = lesq_solve_raim(num_dds, dd_obs, N, DE, b,
    true, DEFAULT_RAIM_THRESHOLD, 0, 0, &bad_index);

  fail_unless(ret == 2,
    "Expecting 1 for repaired solution, got: %i.\n", ret);
}
END_TEST

START_TEST(test_lesq_repair2)
{
  /* Bad DE row, not enough rows to repair. */
  double N[] = {0, 0, 0, 0};
  u8 num_dds = sizeof(N)/sizeof(N[0]);

  double DE[] = {1, 0, 0,
                 0, 1, 0,
                 0, 0, 1,
                 22, 222, 2222};
  double dd_obs[] = {1, 1, 1, 1};

  double b[3];
  s8 ret = lesq_solve_raim(num_dds, dd_obs, N, DE, b,
    false, DEFAULT_RAIM_THRESHOLD, 0, 0, 0);

  fail_unless(ret == -4,
    "Expecting -4 for not enough dds to repair, got: %i.\n", ret);
}
END_TEST

Suite* baseline_test_suite(void)
{
  Suite *s = suite_create("Baseline Calculations");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_predict_carrier_obs);
  tcase_add_test(tc_core, test_predict_carrier_obs2);
  tcase_add_test(tc_core, test_amb_from_baseline);
  tcase_add_test(tc_core, test_lesq_solution);
  tcase_add_test(tc_core, test_lesq_solution2);
  tcase_add_test(tc_core, test_lesq_solution3);
  tcase_add_test(tc_core, test_lesq_solution4);
  tcase_add_test(tc_core, test_lesq_solution5);
  tcase_add_test(tc_core, test_lesq_repair1);
  tcase_add_test(tc_core, test_lesq_repair2);
  tcase_add_test(tc_core, test_lesq_repair8);
  tcase_add_test(tc_core, test_lesq_repair_disabled);

  suite_add_tcase(s, tc_core);

  TCase *tc_baseline = tcase_create("Baseline");
  tcase_add_checked_fixture (tc_baseline, check_baseline_setup,
                                          check_baseline_teardown);
  tcase_add_test(tc_baseline, test_baseline_ref_first);
  tcase_add_test(tc_baseline, test_baseline_ref_middle);
  tcase_add_test(tc_baseline, test_baseline_ref_end);
  tcase_add_test(tc_baseline, test_baseline_fixed_point);
  tcase_add_test(tc_baseline, test_baseline_few_sats);
  suite_add_tcase(s, tc_baseline);

  return s;
}

