#include <check.h>
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

void gmtime_test(const char *name, const utc_params_t *p,
                 time_t start, time_t end, time_t step,
                 time_t offset)
{
  time_t t_gps = start;
  while (t_gps < end) {
    time_t t_unix = t_gps + GPS_EPOCH - offset;
    struct tm *date = gmtime(&t_unix);

    gps_time_t t = {.wn = t_gps / WEEK_SECS,
                    .tow = t_gps % WEEK_SECS};

    utc_time_t u;
    gps2utc(p, &t, &u);

    fail_unless((date->tm_year + 1900) == u.year,
                "%s, expected year %d, got %d, time = %s",
                 name, date->tm_year + 1900, u.year, asctime(date));
    fail_unless((date->tm_mon + 1) == u.month,
                "%s, expected month %d, got %d, time = %s",
                 name, date->tm_mon + 1, u.month, asctime(date));
    fail_unless((date->tm_yday + 1) == u.year_day,
                "%s, expected year_day %d, got %d, time = %s",
                name, date->tm_yday + 1, u.year_day, asctime(date));
    fail_unless(date->tm_mday == u.month_day,
                "%s, expected month_day %d, got %d, time = %s",
                name, date->tm_mday, u.month_day, asctime(date));
    fail_unless(date->tm_wday == (u.week_day % 7),
                "%s, expected week_day %d, got %d, time = %s",
                name, date->tm_wday, u.week_day % 7, asctime(date));
    fail_unless(date->tm_hour == u.hour,
                "%s, expected hour %d, got %d, time = %s",
                name, date->tm_hour, u.hour, asctime(date));
    fail_unless(date->tm_min == u.minute,
                "%s, expected minute %d, got %d, time = %s",
                name, date->tm_min, u.minute, asctime(date));
    fail_unless(date->tm_sec == u.second_int,
                "%s, expected second_int %d, got %d, time = %s",
                name, date->tm_sec, u.second_int, asctime(date));
    fail_unless(0.0 == u.second_frac,
                "%s, expected second_frac %f, got %f, time = %s",
                name, 0.0, u.second_frac, asctime(date));

    t_gps += step;
  }
}

START_TEST(test_gps2utc_time)
{
  utc_params_t p;
  memset(&p, 0, sizeof(p));

  gmtime_test("test_gps2utc_time", &p, 0, 1 * DAY_SECS + 1, 1, 0);
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

START_TEST(test_gps2utc_time_default)
{
  gmtime_test("test_gps2utc_time_default", NULL, DAY_SECS,
              2 * DAY_SECS + 1, 1, GPS_MINUS_UTC_SECS);
}
END_TEST
// TODO add test for negative LS
// TODO make general test function that takes utc parms pointer, start, end and jump per loop

START_TEST(test_gps2utc_time_negative)
{
  utc_params_t p;
  memset(&p, 0, sizeof(p));
  p.dt_ls = -GPS_MINUS_UTC_SECS;
  p.dt_lsf = -GPS_MINUS_UTC_SECS;

  gmtime_test("test_gps2utc_time_negative", &p, DAY_SECS,
              2 * DAY_SECS + 1, 1, -GPS_MINUS_UTC_SECS);
}
END_TEST

START_TEST(test_gps2utc_date)
{
  utc_params_t p;
  memset(&p, 0, sizeof(p));

  gmtime_test("test_gps2utc_date", &p, 0,
              100L * 365 * DAY_SECS, DAY_SECS, 0);
}
END_TEST

START_TEST(test_gps2utc_date_default)
{
  gmtime_test("test_gps2utc_date_default", NULL, DAY_SECS,
              100L * 365 * DAY_SECS, DAY_SECS, GPS_MINUS_UTC_SECS);
}
END_TEST

START_TEST(test_gps2utc_date_negative)
{
  utc_params_t p;
  memset(&p, 0, sizeof(p));
  p.dt_ls = -GPS_MINUS_UTC_SECS;
  p.dt_lsf = -GPS_MINUS_UTC_SECS;

  gmtime_test("test_gps2utc_date_negative", &p, DAY_SECS,
              100L * 365 * DAY_SECS, DAY_SECS, -GPS_MINUS_UTC_SECS);
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
  tcase_add_test(tc_core, test_gps2utc_time_default);
  tcase_add_test(tc_core, test_gps2utc_time_negative);
  tcase_add_test(tc_core, test_gps2utc_date);
  tcase_add_test(tc_core, test_gps2utc_date_negative);
  tcase_add_test(tc_core, test_gps2utc_date_default);
  suite_add_tcase(s, tc_core);

  return s;
}
