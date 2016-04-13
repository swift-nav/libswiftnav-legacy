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
    normalize_gps_time(&testcases[i]);
    double t_normalized = testcases[i].wn * WEEK_SECS + testcases[i].tow;
    fail_unless(fabs(t_original - t_normalized) < tow_tol,
                "normalize_gps_time test case %d failed, t_original = %.12f, "
                "t_normalized = %.12f", i, t_original, t_normalized);
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

START_TEST(test_gps_adjust_week_cycle)
{
  struct gps_adjust_week_cycle_testcase {
    u16 wn_raw, ret;
  } testcases[] = {
    {.wn_raw = 0, .ret = 2048},
    {.wn_raw = 1023, .ret =2047},
    {.wn_raw = GPS_WEEK_REFERENCE % 1024, .ret = GPS_WEEK_REFERENCE},
    {.wn_raw = GPS_WEEK_REFERENCE % 1024 + 1, .ret = GPS_WEEK_REFERENCE + 1},
    {.wn_raw = GPS_WEEK_REFERENCE % 1024 - 1, .ret = GPS_WEEK_REFERENCE + 1023},
    {.wn_raw = GPS_WEEK_REFERENCE, .ret = GPS_WEEK_REFERENCE},
    {.wn_raw = GPS_WEEK_REFERENCE + 1, .ret = GPS_WEEK_REFERENCE + 1},
  };
  const u16 wn_ref = GPS_WEEK_REFERENCE;
  for (size_t i = 0;
       i < sizeof(testcases) / sizeof(struct gps_adjust_week_cycle_testcase); i++) {
    u16 wn = gps_adjust_week_cycle(testcases[i].wn_raw, wn_ref);
    fail_unless(wn == testcases[i].ret,
                "gps_adjust_week_cycle test case %d failed, wn = %d, ret = %d",
                i, wn, testcases[i].ret);
  }
}
END_TEST

START_TEST(test_is_leap_year)
{
  struct is_leap_year_testcase {
    u16 year;
    bool ret;
  } testcases[] = {
    {.year = 1900, .ret = false},
    {.year = 1901, .ret = false},
    {.year = 1904, .ret = true},
    {.year = 1980, .ret = true},
    {.year = 1981, .ret = false},
    {.year = 1982, .ret = false},
    {.year = 1983, .ret = false},
    {.year = 1984, .ret = true},
    {.year = 1985, .ret = false},
    {.year = 1986, .ret = false},
    {.year = 1987, .ret = false},
    {.year = 1988, .ret = true},
    {.year = 1989, .ret = false},
    {.year = 1990, .ret = false},
    {.year = 1991, .ret = false},
    {.year = 1992, .ret = true},
    {.year = 1993, .ret = false},
    {.year = 1994, .ret = false},
    {.year = 1995, .ret = false},
    {.year = 1996, .ret = true},
    {.year = 1997, .ret = false},
    {.year = 1998, .ret = false},
    {.year = 1999, .ret = false},
    {.year = 2000, .ret = true},
    {.year = 2001, .ret = false},
    {.year = 2002, .ret = false},
    {.year = 2003, .ret = false},
    {.year = 2004, .ret = true},
    {.year = 2005, .ret = false},
    {.year = 2006, .ret = false},
    {.year = 2007, .ret = false},
    {.year = 2008, .ret = true},
    {.year = 2009, .ret = false},
    {.year = 2010, .ret = false},
    {.year = 2011, .ret = false},
    {.year = 2012, .ret = true},
    {.year = 2013, .ret = false},
    {.year = 2014, .ret = false},
    {.year = 2015, .ret = false},
    {.year = 2016, .ret = true},
    {.year = 2017, .ret = false},
    {.year = 2018, .ret = false},
    {.year = 2019, .ret = false},
    {.year = 2020, .ret = true},
  };

  for (size_t i = 0;
       i < sizeof(testcases) / sizeof(struct is_leap_year_testcase);
       i++) {
    fail_unless(is_leap_year(testcases[i].year) == testcases[i].ret,
                "is_leap_year test case %d failed, year = %d",
                testcases[i].year);
  }
}
END_TEST

START_TEST(test_glo_time2gps_time)
{
  struct glo_time2gps_time_testcase {
    u16 nt;
    u8 n4;
    s8 h;
    s8 m;
    s8 s;
    gps_time_t ret;
  } testcases[] = {
    /* GLO time 29th Dec 2000 01:00:00 */
    {.nt = 364, .n4 = 2, .h = 1, .m = 0, .s = 0,
      .ret = {.wn = 1094, .tow = 424817}},
    /* GLO time 30th Dec 2000 01:00:00 */
    {.nt = 365, .n4 = 2, .h = 1, .m = 0, .s = 0,
      .ret = {.wn = 1094, .tow = 511217}},
    /* GLO time 31st Dec 2000 02:00:00 */
    {.nt = 366, .n4 = 2, .h = 2, .m = 0, .s = 0,
      .ret = {.wn = 1094, .tow = 601217}},
    /* GLO time 1st Jan  2001 02:00:00 */
    {.nt = 367, .n4 = 2, .h = 2, .m = 0, .s = 0,
      .ret = {.wn = 1095, .tow = 82817}},
    /* GLO time 2nd Jan  2001 02:00:00 */
    {.nt = 368, .n4 = 2, .h = 2, .m = 0, .s = 0,
      .ret = {.wn = 1095, .tow = 169217}},
    /* GLO time 31st Dec 2009 12:12:12 */
    {.nt = 731, .n4 = 4, .h = 12, .m = 12, .s = 12,
      .ret = {.wn = 1564, .tow = 378749}},
    /* GLO time 31st Dec 2010 12:12:12 */
    {.nt = 1096, .n4 = 4, .h = 12, .m = 12, .s = 12,
      .ret = {.wn = 1616, .tow = 465149}},
    /* GLO time 31st Dec 2011 12:12:12 */
    {.nt = 1461, .n4 = 4, .h = 12, .m = 12, .s = 12,
      .ret = {.wn = 1668, .tow = 551549}},
  };
  for (size_t i = 0;
       i < sizeof(testcases) / sizeof(struct glo_time2gps_time_testcase);
       i++) {
    gps_time_t ret = glo_time2gps_time(testcases[i].nt, testcases[i].n4,
                                       testcases[i].h, testcases[i].m,
                                       testcases[i].s);
    fail_unless(ret.wn == testcases[i].ret.wn &&
                ret.tow == testcases[i].ret.tow,
                "glo_time2gps_time test case %d failed, %d, %f", i, ret.wn, ret.tow);
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
  tcase_add_test(tc_core, test_gps_adjust_week_cycle);
  tcase_add_test(tc_core, test_is_leap_year);
  tcase_add_test(tc_core, test_glo_time2gps_time);
  suite_add_tcase(s, tc_core);

  return s;
}
