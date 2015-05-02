
#include <check.h>
#include <stdio.h>
#include "single_diff.h"

navigation_measurement_t nm1 = {.prn = 1};
navigation_measurement_t nm2 = {.prn = 2};
navigation_measurement_t nm3 = {.prn = 3};
navigation_measurement_t nm4 = {.prn = 4};

navigation_measurement_t nms_no_match1[2];
navigation_measurement_t nms_no_match2[2];

START_TEST(test_no_match)
{
    sdiff_t sds_out[6];

    /* Test for when they are interleaved */
    memcpy(&nms_no_match1[0], &nm1, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match2[0], &nm2, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match1[1], &nm3, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match2[1], &nm4, sizeof(navigation_measurement_t));

    u8 num_match = single_diff(2, nms_no_match1,
                               2, nms_no_match2,
                               sds_out);

    fail_unless(num_match == 0);

    /* Test when one set follows the other */
    memcpy(&nms_no_match1[0], &nm1, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match1[1], &nm2, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match2[0], &nm3, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match2[1], &nm4, sizeof(navigation_measurement_t));

    num_match = single_diff(2, nms_no_match1,
                            2, nms_no_match2,
                            sds_out);

    fail_unless(num_match == 0);

    /* Test it the other way */
    memcpy(&nms_no_match2[0], &nm1, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match2[1], &nm2, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match1[0], &nm3, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match1[1], &nm4, sizeof(navigation_measurement_t));

    num_match = single_diff(2, nms_no_match1,
                            2, nms_no_match2,
                            sds_out);

    fail_unless(num_match == 0);
}
END_TEST

START_TEST(test_beginning_matches)
{
    sdiff_t sds_out[3];

    /* Test for when they both have two */
    memcpy(&nms_no_match1[0], &nm1, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match1[1], &nm2, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match2[0], &nm1, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match2[1], &nm3, sizeof(navigation_measurement_t));

    u8 num_match = single_diff(2, nms_no_match1,
                               2, nms_no_match2,
                               sds_out);

    fail_unless(num_match == 1);
    fail_unless(sds_out[0].prn == 1);

    /* Test for both with two the other way */
    memcpy(&nms_no_match1[0], &nm1, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match1[1], &nm3, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match2[0], &nm1, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match2[1], &nm2, sizeof(navigation_measurement_t));

    num_match = single_diff(2, nms_no_match1,
                            2, nms_no_match2,
                            sds_out);

    fail_unless(num_match == 1);
    fail_unless(sds_out[0].prn == 1);

    /* Test when one has only one */
    memcpy(&nms_no_match1[0], &nm1, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match2[0], &nm1, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match2[1], &nm2, sizeof(navigation_measurement_t));

    num_match = single_diff(1, nms_no_match1,
                            2, nms_no_match2,
                            sds_out);

    fail_unless(num_match == 1);
    fail_unless(sds_out[0].prn == 1);

    /* Test when the other has only one */
    memcpy(&nms_no_match1[0], &nm1, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match1[1], &nm2, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match2[0], &nm1, sizeof(navigation_measurement_t));

    num_match = single_diff(2, nms_no_match1,
                            1, nms_no_match2,
                            sds_out);

    fail_unless(num_match == 1);
    fail_unless(sds_out[0].prn == 1);
}
END_TEST

START_TEST(test_end_matches)
{
    sdiff_t sds_out[3];

    /* Test for when they both have two */
    memcpy(&nms_no_match1[0], &nm1, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match1[1], &nm3, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match2[0], &nm2, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match2[1], &nm3, sizeof(navigation_measurement_t));

    u8 num_match = single_diff(2, nms_no_match1,
                               2, nms_no_match2,
                               sds_out);

    fail_unless(num_match == 1);
    fail_unless(sds_out[0].prn == 3);

    /* Test for both with two the other way */
    memcpy(&nms_no_match1[0], &nm2, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match1[1], &nm3, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match2[0], &nm1, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match2[1], &nm3, sizeof(navigation_measurement_t));

    num_match = single_diff(2, nms_no_match1,
                            2, nms_no_match2,
                            sds_out);

    fail_unless(num_match == 1);
    fail_unless(sds_out[0].prn == 3);

    /* Test when one has only one */
    memcpy(&nms_no_match1[0], &nm2, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match2[0], &nm1, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match2[1], &nm2, sizeof(navigation_measurement_t));

    num_match = single_diff(1, nms_no_match1,
                            2, nms_no_match2,
                            sds_out);

    fail_unless(num_match == 1);
    fail_unless(sds_out[0].prn == 2);

    /* Test when the other has only one */
    memcpy(&nms_no_match1[0], &nm1, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match1[1], &nm2, sizeof(navigation_measurement_t));
    memcpy(&nms_no_match2[0], &nm2, sizeof(navigation_measurement_t));

    num_match = single_diff(2, nms_no_match1,
                            1, nms_no_match2,
                            sds_out);

    fail_unless(num_match == 1);
    fail_unless(sds_out[0].prn == 2);
}
END_TEST

START_TEST(test_is_prn_set)
{

#define TEST_IS_SET(set, result) \
  fail_unless(is_prn_set(sizeof(set)/sizeof(set[0]), set) == result, \
              "is_prn_set(" #set ") != " #result);

  /* Normal set. */
  u8 prns1[] = {0, 1, 2, 33, 44, 200};
  TEST_IS_SET(prns1, true);

  /* Empty set. */
  fail_unless(is_prn_set(0, prns1) == true);

  /* Single element set. */
  u8 prns2[] = {22};
  TEST_IS_SET(prns2, true);

  /* Repeated elements. */

  u8 prns3[] = {22, 22};
  TEST_IS_SET(prns3, false);

  u8 prns4[] = {0, 1, 2, 3, 3};
  TEST_IS_SET(prns4, false);

  u8 prns5[] = {1, 1, 2, 3, 4};
  TEST_IS_SET(prns5, false);

  /* Incorrectly sorted. */

  u8 prns6[] = {22, 1, 2, 3, 4};
  TEST_IS_SET(prns6, false);

  u8 prns7[] = {0, 1, 2, 3, 1};
  TEST_IS_SET(prns7, false);

  u8 prns8[] = {0, 1, 22, 3, 4};
  TEST_IS_SET(prns8, false);
}
END_TEST

Suite* sdiff_test_suite(void)
{
  Suite *s = suite_create("Single Differences");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_no_match);
  tcase_add_test(tc_core, test_beginning_matches);
  tcase_add_test(tc_core, test_end_matches);
  tcase_add_test(tc_core, test_is_prn_set);
  suite_add_tcase(s, tc_core);

  return s;
}

