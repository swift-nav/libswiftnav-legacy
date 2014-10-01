
#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ambiguity_test.h>

#include "check_utils.h"


// assure that when the sdiffs match amb_test's sats, amb_test's sats are unchanged
START_TEST(test_update_sats_same_sats)
{
  ambiguity_test_t amb_test = {.sats = {.num_sats = 4,
                                        .prns = {3,1,2,4}}};
  sdiff_t sdiffs[4] = {{.prn = 1},
                       {.prn = 2},
                       {.prn = 3},
                       {.prn = 4}};
  u8 num_sdiffs = 4;

  ambiguity_update_sats(&amb_test, num_sdiffs, sdiffs, NULL, NULL, NULL, NULL);

  fail_unless(amb_test.sats.prns[0] == 3);
  fail_unless(amb_test.sats.prns[1] == 1);
  fail_unless(amb_test.sats.prns[2] == 2);
  fail_unless(amb_test.sats.prns[3] == 4);
}
END_TEST

//assure that when we've lost the reference, we choose a new one and rebase everything
START_TEST(test_update_sats_rebase)
{
  ambiguity_test_t amb_test;
  create_ambiguity_test(&amb_test);

  amb_test.sats.num_sats = 4;
  amb_test.sats.prns[0] = 3;
  amb_test.sats.prns[1] = 1;
  amb_test.sats.prns[2] = 2;
  amb_test.sats.prns[3] = 4;

  sdiff_t sdiffs[3] = {{.prn = 1, .snr = 0},
                       {.prn = 2, .snr = 0}, 
                      // {.prn = 3, .snr = 0}, 
                       {.prn = 4, .snr = 1}};
  u8 num_sdiffs = 3;
  
  hypothesis_t *hyp = (hypothesis_t *)memory_pool_add(amb_test.pool);
  hyp->N[0] = 0;
  hyp->N[1] = 1;
  hyp->N[2] = 2;

  sats_management_t float_sats = {.num_sats = 3};

  ambiguity_update_sats(&amb_test, num_sdiffs, sdiffs, &float_sats, NULL, NULL, NULL);
  fail_unless(amb_test.sats.num_sats == 3);
  fail_unless(amb_test.sats.prns[0] == 4);
  fail_unless(amb_test.sats.prns[1] == 1);
  fail_unless(amb_test.sats.prns[2] == 2);
  fail_unless(hyp->N[0] == -2);
  fail_unless(hyp->N[1] == -1);
}
END_TEST

//void ambiguity_update_reference(ambiguity_test_t *amb_test, u8 num_sdiffs, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first);
START_TEST(test_ambiguity_update_reference)
{
  srandom(1);

  ambiguity_test_t amb_test = {.sats = {.num_sats = 4, 
                                        .prns = {3,1,2,4}}};
  create_ambiguity_test(&amb_test);

  sdiff_t sdiffs[4] = {{.prn = 1, .snr = 0},
                       {.prn = 2, .snr = 0}, 
                       // {.prn = 3, .snr = 0}, 
                       {.prn = 4, .snr = 1}};
  u8 num_sdiffs = 4;

  for (u32 i=0; i<3; i++) {
    hypothesis_t *hyp = (hypothesis_t *)memory_pool_add(amb_test.pool);
    fail_unless(hyp != 0, "Null pointer returned by memory_pool_add");
    for (u8 j=0; j<amb_test.sats.num_sats-1; j++) {
      hyp->N[j] = sizerand(5);
    }
    hyp->ll = frand(0, 1);
  }

  u8 num_dds = MAX(0,amb_test.sats.num_sats - 1);
  printf("Before rebase:\n");
  memory_pool_map(amb_test.pool, &num_dds, &print_hyp);

  sdiff_t sdiffs_with_ref_first[4];
  ambiguity_update_reference(&amb_test, num_sdiffs, sdiffs, sdiffs_with_ref_first);

  printf("After rebase:\n");
  memory_pool_map(amb_test.pool, &num_dds, &print_hyp);
}
END_TEST

//s8 sats_match(ambiguity_test_t *amb_test, u8 num_sdiffs, sdiff_t *sdiffs);
START_TEST(test_sats_match)
{
  ambiguity_test_t amb_test = {.sats = {.num_sats = 3,
                                        .prns = {3,1,2}}};
  sdiff_t sdiffs[4] = {{.prn = 1},
                       {.prn = 2},
                       {.prn = 3},
                       {.prn = 4}};
  u8 num_sdiffs = 4;
  fail_unless(!sats_match(&amb_test, num_sdiffs, sdiffs));

  num_sdiffs = 3;
  fail_unless(sats_match(&amb_test, num_sdiffs, sdiffs));

  sdiffs[0].prn = 22;
  fail_unless(!sats_match(&amb_test, num_sdiffs, sdiffs));

}
END_TEST


Suite* ambiguity_test_suite(void)
{
  Suite *s = suite_create("Ambiguity Test");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_sats_match);
  tcase_add_test(tc_core, test_ambiguity_update_reference);
  tcase_add_test(tc_core, test_update_sats_same_sats);
  tcase_add_test(tc_core, test_update_sats_rebase);
  suite_add_tcase(s, tc_core);

  return s;
}

