#include <check.h>
#include <stdio.h>
#include <math.h>
#include "check_utils.h"

#include "gpstime.h"

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

Suite* gpstime_test_suite(void)
{
  Suite *s = suite_create("GPS time handling");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_gpsdifftime);
  suite_add_tcase(s, tc_core);

  return s;
}

