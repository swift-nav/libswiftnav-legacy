#include <check.h>
#include <libswiftnav/track.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define BW   1000
#define CN0_0 40
#define CUTOFF_FREQ 0.1
#define LOOP_FREQ 1000


START_TEST(test_lp1_init1)
{
  lp1_filter_t filter;
  lp1_filter_init(&filter, 0.f, 0.5f, 1e3f);

  fail_unless(filter.xn == 0.f);
  fail_unless(filter.yn == 0.f);
  fail_unless(filter.a != 0.f);
  fail_unless(filter.b != 0.f);
}
END_TEST

START_TEST(test_lp1_init2)
{
  lp1_filter_t filter;
  lp1_filter_init(&filter, 1.f, 0.5f, 1e3f);

  fail_unless(filter.xn == 1.f * filter.b);
  fail_unless(filter.yn == 1.f);
  fail_unless(filter.a != 0.f);
  fail_unless(filter.b != 0.f);
}
END_TEST

START_TEST(test_lp1_update1)
{
  lp1_filter_t filter;
  lp1_filter_init(&filter, 1.f, 0.5f, 1e3f);
  float r;

  r = lp1_filter_update(&filter, 1.);
  fail_unless(r == 1.f);
}
END_TEST

START_TEST(test_lp1_update2)
{
  lp1_filter_t filter;
  lp1_filter_init(&filter, 1.f, 0.5f, 1e3f);
  float r;

  r = lp1_filter_update(&filter, 0.f);
  fail_unless(r < 1.f);
}
END_TEST

START_TEST(test_lp1_update3)
{
  lp1_filter_t filter;
  lp1_filter_init(&filter, 1.f, 0.5f, 1e3f);
  float r;

  r = lp1_filter_update(&filter, 1.1f);
  fail_unless(r > 1.f);
}
END_TEST

START_TEST(test_bw2_init1)
{
  bw2_filter_t filter;
  bw2_filter_init(&filter, 0.f, 0.5f, 1e3f);

  fail_unless(filter.xn == 0.f);
  fail_unless(filter.xn_prev == 0.f);
  fail_unless(filter.yn == 0.f);
  fail_unless(filter.yn_prev == 0.f);
  fail_unless(filter.a2 != 0.f);
  fail_unless(filter.a3 != 0.f);
  fail_unless(filter.b != 0.f);
}
END_TEST

START_TEST(test_bw2_init2)
{
  bw2_filter_t filter;
  bw2_filter_init(&filter, 1.f, 0.5f, 1e3f);

  fail_unless(filter.xn == 1.f);
  fail_unless(filter.xn_prev == 1.f);
  fail_unless(filter.yn == 1.f);
  fail_unless(filter.yn_prev == 1.f);
  fail_unless(filter.a2 != 0.f);
  fail_unless(filter.a3 != 0.f);
  fail_unless(filter.b != 0.f);
}
END_TEST

START_TEST(test_bw2_update1)
{
  bw2_filter_t filter;
  float r;

  bw2_filter_init(&filter, 1.f, 0.5f, 1e3f);
  r = bw2_filter_update(&filter, 1.f);

  fail_unless(fabsf(r - 1.f) < 1e-6);
}
END_TEST

START_TEST(test_bw2_update2)
{
  bw2_filter_t filter;
  float r;

  bw2_filter_init(&filter, 1.f, 0.5f, 1e3f);
  r = bw2_filter_update(&filter, .9);

  fail_unless(r < 1.f);
}
END_TEST

START_TEST(test_bw2_update3)
{
  bw2_filter_t filter;
  float r;

  bw2_filter_init(&filter, 1.f, 0.5f, 1e3f);
  r = bw2_filter_update(&filter, 1.1f);

  fail_unless(r > 1.f);
}
END_TEST


Suite* cn0_filter_suite(void)
{
  Suite *s = suite_create("CN0 Filters");
  TCase *tc_core = tcase_create("Core");

  tcase_add_test(tc_core, test_lp1_init1);
  tcase_add_test(tc_core, test_lp1_init2);
  tcase_add_test(tc_core, test_lp1_update1);
  tcase_add_test(tc_core, test_lp1_update2);
  tcase_add_test(tc_core, test_lp1_update3);
  tcase_add_test(tc_core, test_bw2_init1);
  tcase_add_test(tc_core, test_bw2_init2);
  tcase_add_test(tc_core, test_bw2_update1);
  tcase_add_test(tc_core, test_bw2_update2);
  tcase_add_test(tc_core, test_bw2_update3);
  suite_add_tcase(s, tc_core);

  return s;
}
