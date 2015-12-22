#include <check.h>
#include <stdio.h>
#include <math.h>

#include <libswiftnav/time.h>
#include <libswiftnav/constants.h>

#include "check_utils.h"

START_TEST(test_gpsdifftime)
{
  struct gpsdifftime_testcase {
    gps_time_t a, b;
    double dt;
  } testcases[] = {
    {.a = {567890.0, 1234}, .b = {567890.0, 1234}, .dt = 0},
    {.a = {567890.0, 1234}, .b = {0.0, 1234}, .dt = 567890},
    {.a = {567890.0, WN_UNKNOWN}, .b = {0.0, 1234}, .dt = -36910},
    {.a = {222222.0, 2222}, .b = {2222.0, WN_UNKNOWN}, .dt = 220000},
    {.a = {444444.0, WN_UNKNOWN}, .b = {2222.0, WN_UNKNOWN}, .dt = -162578},
    {.a = {604578.0, 1000}, .b = {222.222, 1001}, .dt = -444.222},
    {.a = {604578.0, 1001}, .b = {222.222, 1000}, .dt = 1209155.778},
  };
  const double tow_tol = 1e-10;
  for (size_t i = 0;
       i < sizeof(testcases) / sizeof(struct gpsdifftime_testcase); i++) {
    double dt = gpsdifftime(testcases[i].a, testcases[i].b);
    fail_unless(fabs(dt - testcases[i].dt) < tow_tol,
                "gpsdifftime test case %d failed, dt = %.12f", i, dt);
  }
}
END_TEST

START_TEST(test_normalize_gps_time)
{
  gps_time_t testcases[] = {
    {0, 1234},
    {3 * DAY_SECS, 1234},
    {WEEK_SECS + DAY_SECS, 1234},
    {0 - DAY_SECS, 1234}
  };
  const double tow_tol = 1e-10;
  for (size_t i = 0;
       i < sizeof(testcases) / sizeof(gps_time_t); i++) {
    double t_original = testcases[i].wn * WEEK_SECS + testcases[i].tow;
    gps_time_t normalized = normalize_gps_time(testcases[i]);
    double t_normalized = normalized.wn * WEEK_SECS + normalized.tow;
    fail_unless(fabs(t_original - t_normalized) < tow_tol,
                "normalize_gps_time test case %d failed, t_original = %.12f, "
                "t_normalized = %.12f", i, t_original, t_normalized);
  }
}
END_TEST

Suite* time_test_suite(void)
{
  Suite *s = suite_create("Time handling");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_gpsdifftime);
  tcase_add_test(tc_core, test_normalize_gps_time);
  suite_add_tcase(s, tc_core);

  return s;
}
