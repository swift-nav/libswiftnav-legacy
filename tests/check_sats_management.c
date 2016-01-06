#include <check.h>
#include <stdio.h>

#include <libswiftnav/linear_algebra.h>
#include <libswiftnav/amb_kf.h>

#include "check_utils.h"

/* Need static method set_reference_sat */
#include "sats_management.c"

START_TEST(test_rebase_1)
{
  gnss_signal_t prns[] = {{.sat = 2},{.sat = 1},{.sat = 3},{.sat = 4}};
  u8 num_sats = sizeof(prns)/sizeof(gnss_signal_t);
  gnss_signal_t new_ref = {.sat = 3};

  sats_management_t sats_management;
  sats_management.num_sats = 4;
  sats_management.sids[0].sat = 2;
  sats_management.sids[1].sat = 1;
  sats_management.sids[2].sat = 3;
  sats_management.sids[3].sat = 4;

  set_reference_sat_of_sids(new_ref, num_sats, prns);
  /* Just check the sats_management prns update */
  set_reference_sat(new_ref, &sats_management, 0, 0, 0);
  /* TODO make a better test */
}
END_TEST

Suite* sats_management_test_suite(void)
{
  Suite *s = suite_create("Sats Management");

  TCase *tc_rebase = tcase_create("rebase");
  tcase_add_test(tc_rebase, test_rebase_1);
  suite_add_tcase(s, tc_rebase);

  return s;
}
