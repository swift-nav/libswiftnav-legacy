
#include <check.h>
#include <string.h>
#include <stdio.h>

#include <constants.h>
#include <linear_algebra.h>
#include <observation.h>
#include <filter_utils.h>

START_TEST(test_assign_de_mtx_1)
{
  sdiff_t sdiffs[] = {
    {.sat_pos = {23, 22, 22}},
    {.sat_pos = {22, 23, 22}},
    {.sat_pos = {22, 22, 23}}
  };
  u8 num_sats = sizeof(sdiffs) / sizeof(sdiffs[0]);
  double ref_ecef[3] = {22, 22, 22};
  double DE[(num_sats-1) * 3];

  double DE_expected[] = {
    -1, 1, 0,
    -1, 0, 1
  };

  fail_unless(assign_de_mtx(num_sats, sdiffs, ref_ecef, DE) == 0);

  fail_unless(memcmp(DE, DE_expected, sizeof(DE)) == 0, "DE matrix incorrect");
}
END_TEST

START_TEST(test_assign_de_mtx_2)
{
  sdiff_t sdiffs[] = {
    {.sat_pos = {23, 22, 22}},
    {.sat_pos = {22, 23, 22}},
  };
  u8 num_sats = sizeof(sdiffs) / sizeof(sdiffs[0]);
  double ref_ecef[3] = {22, 22, 22};
  double DE[(num_sats-1) * 3];

  double DE_expected[] = {
    -1, 1, 0,
  };

  fail_unless(assign_de_mtx(num_sats, sdiffs, ref_ecef, DE) == 0);

  fail_unless(memcmp(DE, DE_expected, sizeof(DE)) == 0, "DE matrix incorrect");
}
END_TEST

START_TEST(test_assign_de_mtx_3)
{
  sdiff_t sdiffs[] = {};
  double ref_ecef[] = {};
  double DE[] = {};

  fail_unless(assign_de_mtx(1, sdiffs, ref_ecef, DE) == -1);
  fail_unless(assign_de_mtx(0, sdiffs, ref_ecef, DE) == -1);
}
END_TEST

START_TEST(test_simple_amb_measurement)
{
  fail_unless(simple_amb_measurement(0, 0) == 0);
  fail_unless(simple_amb_measurement(1234, 0) == 1234);
  fail_unless(simple_amb_measurement(0, GPS_L1_LAMBDA_NO_VAC) == 1);
  fail_unless(simple_amb_measurement(22, 33) == 22 + (33 / GPS_L1_LAMBDA_NO_VAC));
}
END_TEST

Suite* filter_utils_suite(void)
{
  Suite *s = suite_create("Filter Utils");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_assign_de_mtx_1);
  tcase_add_test(tc_core, test_assign_de_mtx_2);
  tcase_add_test(tc_core, test_assign_de_mtx_3);
  tcase_add_test(tc_core, test_simple_amb_measurement);
  suite_add_tcase(s, tc_core);

  return s;
}

