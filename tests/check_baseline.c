
#include <check.h>
#include <math.h>

#include <linear_algebra.h>
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
                "Observation mismatch: %ld vs %lf",
                N[i], N_true[i]);
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
  suite_add_tcase(s, tc_core);

  return s;
}

