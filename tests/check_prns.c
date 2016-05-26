#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "check_utils.h"

#include <libswiftnav/constants.h>
#include <libswiftnav/prns.h>

/* Include prns.c here to have a chance to turn off asserts.
   Otherwise the code lines with asserts cannot be covered by
   tests and it will lower test coverage statistics. */
#define NDEBUG
#include <prns.c>

START_TEST(test_prns_sid_to_init_g1)
{
  u32 g1;
  gnss_signal_t sid;

  sid = construct_sid(CODE_GPS_L2CM, 25);
  g1 = sid_to_init_g1(sid);
  fail_unless(0700274134 == g1);

  sid = construct_sid(CODE_GPS_L1CA, 25);
  g1 = sid_to_init_g1(sid);
  fail_unless(0x3ff == g1);

  /* check unsupported branch for code coverage stats */
  sid = construct_sid(CODE_GLO_L2CA, 25);
  g1 = sid_to_init_g1(sid);
}
END_TEST

START_TEST(test_prns_ca_code)
{
  const u8* ptr;
  gnss_signal_t sid;

  sid = construct_sid(CODE_GPS_L2CM, 25);
  ptr = ca_code(sid);
  fail_unless(NULL == ptr);

  sid = construct_sid(CODE_GPS_L1CA, 1);
  ptr = ca_code(sid);
  fail_unless(NULL != ptr);
}
END_TEST

START_TEST(test_prns_get_chip)
{
  u8 code = 0xFF;

  s8 chip;

  chip = get_chip(&code, 0);
  fail_unless(-1 == chip);
}
END_TEST

Suite* prns_test_suite(void)
{
  Suite *s = suite_create("Prns");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_prns_sid_to_init_g1);
  tcase_add_test(tc_core, test_prns_ca_code);
  tcase_add_test(tc_core, test_prns_get_chip);
  suite_add_tcase(s, tc_core);

  return s;
}
