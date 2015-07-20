#include <check.h>
#include <stdio.h>
#include "linear_algebra.h"
#include "check_utils.h"
#include "sats_management.c"
#include "amb_kf.h"


START_TEST(test_rebase_1)
{
  u8 prns[4] = {2,1,3,4};
  u8 num_sats = sizeof(prns);
  u8 new_ref = 3;

  sats_management_t sats_management;
  sats_management.num_sats = 4;
  sats_management.prns[0] = 2;
  sats_management.prns[1] = 1;
  sats_management.prns[2] = 3;
  sats_management.prns[3] = 4;

  set_reference_sat_of_prns(new_ref, num_sats, prns);
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
