
#include <check.h>
#include <math.h>
#include <stdio.h>

#include <linear_algebra.h>
#include <constants.h>
#include <baseline.h>

#define TOL 1e-10

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
    fail_unless(fabs(dd_obs[i] - dd_obs_expected[i]) < TOL,
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

  fail_unless(fabs(dd_obs[0] - dd_obs_expected[0]) < TOL,
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
    fail_unless(fabs(N_true[i] - N[i]) < TOL,
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
    fail_unless(fabs(b[i] - b_true[i]) < TOL,
                "Baseline mismatch: %lf vs %lf",
                b[i], b_true[i]);
  }

  /* Try with resid = NULL */
  ret = lesq_solution_float(num_dds, dd_obs, N, DE, b, 0);

  fail_unless(ret == 0, "solution returned error %d", ret);

  for (u8 i=0; i<3; i++) {
    fail_unless(fabs(b[i] - b_true[i]) < TOL,
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
    fail_unless(fabs(b[i] - b_true[i]) < TOL,
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
  s32 N_int[num_dds];
  for (u8 i=0; i<num_dds; i++) {
    N_int[i] = N[i];
  }
  s8 ret = lesq_solution_int(num_dds, dd_obs, N_int, DE, b, resid);

  fail_unless(ret == 0, "solution returned error %d", ret);

  for (u8 i=0; i<3; i++) {
    fail_unless(fabs(b[i] - b_true[i]) < TOL,
                "Baseline mismatch: %lf vs %lf",
                b[i], b_true[i]);
  }

  /* Try with resid = NULL */
  ret = lesq_solution_int(num_dds, dd_obs, N_int, DE, b, 0);

  fail_unless(ret == 0, "solution returned error %d", ret);

  for (u8 i=0; i<3; i++) {
    fail_unless(fabs(b[i] - b_true[i]) < TOL,
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
    fail_unless(fabs(b[i] - b_expected[i]) < TOL,
                "Baseline mismatch: %lf vs %lf",
                b[i], b_true[i]);
  }
  for (u8 i=0; i<num_dds; i++) {
    fail_unless(fabs(resid[i] - resid_expected[i]) < TOL,
                "Residual mismatch: %lf vs %lf",
                resid[i], resid_expected[i]);
  }
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
  suite_add_tcase(s, tc_core);

  return s;
}

