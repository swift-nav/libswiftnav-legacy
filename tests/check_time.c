#include <check.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

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
    double dt = gpsdifftime(&testcases[i].a, &testcases[i].b);
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

START_TEST(test_gps2utc_time)
{
  utc_params_t p;
  memset(&p, 0, sizeof(p));

  time_t t_gps = 0;
  while (t_gps < (1 * DAY_SECS + 1)) {
    time_t t_unix = t_gps + 315964800;
    struct tm *date = gmtime(&t_unix);

    gps_time_t t = {.wn = t_gps / WEEK_SECS, .tow = t_gps % WEEK_SECS};

    utc_time_t u;
    gps2utc(&p, &t, &u);

    fail_unless((date->tm_year + 1900) == u.year,
                "gmtime test, expected year %d, got %d, time = %s",
                date->tm_year + 1900, u.year, asctime(date));
    fail_unless((date->tm_mon + 1) == u.month,
                "gmtime test, expected month %d, got %d, time = %s",
                date->tm_mon + 1, u.month, asctime(date));
    fail_unless((date->tm_yday + 1) == u.year_day,
                "gmtime test, expected year_day %d, got %d, time = %s",
                date->tm_yday + 1, u.year_day, asctime(date));
    fail_unless(date->tm_mday == u.month_day,
                "gmtime test, expected month_day %d, got %d, time = %s",
                date->tm_mday, u.month_day, asctime(date));
    fail_unless(date->tm_wday == (u.week_day % 7),
                "gmtime test, expected week_day %d, got %d, time = %s",
                date->tm_wday, u.week_day % 7, asctime(date));
    fail_unless(date->tm_hour == u.hour,
                "gmtime test, expected hour %d, got %d, time = %s",
                date->tm_hour, u.hour, asctime(date));
    fail_unless(date->tm_min == u.minute,
                "gmtime test, expected minute %d, got %d, time = %s",
                date->tm_min, u.minute, asctime(date));
    fail_unless(date->tm_sec == u.second_int,
                "gmtime test, expected second_int %d, got %d, time = %s",
                date->tm_sec, u.second_int, asctime(date));
    fail_unless(0.0 == u.second_frac,
                "gmtime test, expected second_frac %f, got %f, time = %s",
                0.0, u.second_frac, asctime(date));

    t_gps++;
  }
}
END_TEST

START_TEST(test_gps_time_match_weeks)
{
  struct gps_time_match_weeks_testcase {
    gps_time_t t, ref, ret;
  } testcases[] = {
    {.t = {0.0, WN_UNKNOWN}, .ref = {0.0, 1234}, .ret = {0.0, 1234}},
    {.t = {WEEK_SECS-1, WN_UNKNOWN}, .ref = {0.0, 1234}, .ret = {WEEK_SECS-1, 1233}},
    {.t = {0.0, WN_UNKNOWN}, .ref = {WEEK_SECS-1, 1234}, .ret = {0.0, 1235}},
    {.t = {WEEK_SECS-1, WN_UNKNOWN}, .ref = {WEEK_SECS-1, 1234}, .ret = {WEEK_SECS-1, 1234}},
    {.t = {2*DAY_SECS, WN_UNKNOWN}, .ref = {5*DAY_SECS, 1234}, .ret = {2*DAY_SECS, 1234}},
    {.t = {5*DAY_SECS, WN_UNKNOWN}, .ref = {2*DAY_SECS, 1234}, .ret = {5*DAY_SECS, 1234}},
    {.t = {0.0, WN_UNKNOWN}, .ref = {WEEK_SECS/2, 1234}, .ret = {0.0, 1234}},
    {.t = {WEEK_SECS/2, WN_UNKNOWN}, .ref = {0.0, 1234}, .ret = {WEEK_SECS/2, 1234}},
    {.t = {WEEK_SECS/2+1, WN_UNKNOWN}, .ref = {0.0, 1234}, .ret = {WEEK_SECS/2+1, 1233}},
    {.t = {0.0, WN_UNKNOWN}, .ref = {WEEK_SECS/2+1, 1234}, .ret = {0.0, 1235}},
  };
  for (size_t i = 0;
       i < sizeof(testcases) / sizeof(struct gps_time_match_weeks_testcase); i++) {
    gps_time_match_weeks(&testcases[i].t, &testcases[i].ref);
    fail_unless(testcases[i].t.wn == testcases[i].ret.wn,
                "gps_time_match_weeks test case %d failed, t.wn = %d, ret.wn = %d",
                i, testcases[i].t.wn, testcases[i].ret.wn);
    fail_unless(testcases[i].t.tow == testcases[i].ret.tow,
                "gps_time_match_weeks test case %d failed, t.tow = %.12f, ret.tow = %.12f",
                i, testcases[i].t.tow, testcases[i].ret.tow);
  }
}
END_TEST

START_TEST(test_gps2utc_date)
{
  utc_params_t p;
  memset(&p, 0, sizeof(p));

  time_t t_gps = 0;
  while (t_gps < 100L * 365 * DAY_SECS) {
    time_t t_unix = t_gps + 315964800;
    struct tm *date = gmtime(&t_unix);

    gps_time_t t = {.wn = t_gps / WEEK_SECS, .tow = t_gps % WEEK_SECS};

    utc_time_t u;
    gps2utc(&p, &t, &u);

    fail_unless((date->tm_year + 1900) == u.year,
                "gmtime test, expected year %d, got %d, time = %s",
                date->tm_year + 1900, u.year, asctime(date));
    fail_unless((date->tm_mon + 1) == u.month,
                "gmtime test, expected month %d, got %d, time = %s",
                date->tm_mon + 1, u.month, asctime(date));
    fail_unless((date->tm_yday + 1) == u.year_day,
                "gmtime test, expected year_day %d, got %d, time = %s",
                date->tm_yday + 1, u.year_day, asctime(date));
    fail_unless(date->tm_mday == u.month_day,
                "gmtime test, expected month_day %d, got %d, time = %s",
                date->tm_mday, u.month_day, asctime(date));
    fail_unless(date->tm_wday == (u.week_day % 7),
                "gmtime test, expected week_day %d, got %d, time = %s",
                date->tm_wday, u.week_day % 7, asctime(date));
    fail_unless(date->tm_hour == u.hour,
                "gmtime test, expected hour %d, got %d, time = %s",
                date->tm_hour, u.hour, asctime(date));
    fail_unless(date->tm_min == u.minute,
                "gmtime test, expected minute %d, got %d, time = %s",
                date->tm_min, u.minute, asctime(date));
    fail_unless(date->tm_sec == u.second_int,
                "gmtime test, expected second_int %d, got %d, time = %s",
                date->tm_sec, u.second_int, asctime(date));
    fail_unless(0.0 == u.second_frac,
                "gmtime test, expected second_frac %f, got %f, time = %s",
                0.0, u.second_frac, asctime(date));

    t_gps += DAY_SECS;
  }
}
END_TEST

Suite* time_test_suite(void)
{
  Suite *s = suite_create("Time handling");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_gpsdifftime);
  tcase_add_test(tc_core, test_normalize_gps_time);
  tcase_add_test(tc_core, test_gps_time_match_weeks);
  tcase_add_test(tc_core, test_gps2utc_time);
  tcase_add_test(tc_core, test_gps2utc_date);
  suite_add_tcase(s, tc_core);

  return s;
}
