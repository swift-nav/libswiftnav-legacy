#include <stdio.h>
#include <check.h>
#include "check_utils.h"

#include <time.h>
#include <stdlib.h>

#include <set.h>
#include <iterator.h>

bool arr_eq(size_t len, prn *arr1, prn *arr2)
{
  for (size_t i = 0; i < len; i++) {
    if (arr1[i] != arr2[i]) {
      return false;
    }
  }
  return true;
}

START_TEST(test_set)
{
  prn prns1[] = {};
  set_t set1;
  mk_prn_set(&set1, 0, prns1);
  fail_unless(is_set(&set1), "Empty set should be a valid set");

  prn prns2[] = {0,1,3};
  set_t set2;
  mk_prn_set(&set2, 3, prns2);
  fail_unless(is_set(&set2), "0,1,3 should be a valid set");

  prn prns3[] = {1,0};
  set_t set3;
  mk_prn_set(&set3, 2, prns3);
  fail_unless(!is_set(&set3), "is_set(Disordered set) should return false.");

}
END_TEST

START_TEST(test_freeze)
{
  prn prns2[] = {0,1,3};
  set_t set2;
  mk_prn_set(&set2, 3, prns2);

  prn frozen2[3];
  iterator_t it2;
  set_state_t s2;
  mk_set_itr(&it2, &s2, &set2);
  freeze_arr(sizeof(prn), 22, frozen2, &it2);
  fail_unless(arr_eq(3, frozen2, prns2), "Frozen array does not match input to iterator.");

  fail_unless(freeze_arr(sizeof(prn), 2, frozen2, &it2) == -1, "Iterator exceeds max length?");
}
END_TEST
START_TEST(test_filter)
{
  prn prns1[] = {0,22,3};
  prn prns2[] = {0,2,3};
  set_t set;
  mk_prn_set(&set, 3, prns1);
  iterator_t it;
  set_state_t ss;
  mk_set_itr(&it, &ss, &set);

  iterator_t fitr;
  filter_state_t fs;

  prn tt = 22;

  mk_filter_itr(&fitr, &fs, &eq_ref, &tt, &it);
  for(reset(&fitr); more(&fitr); next(&fitr)) {
    printf("itr val1: %i\n", *((prn *)current(&fitr)));
  }
  mk_filter_itr(&fitr, &fs, &not_eq_ref, &tt, &it);
  for(reset(&fitr); more(&fitr); next(&fitr)) {
    printf("itr val2: %i\n", *((prn *)current(&fitr)));
  }

  iterator_t it2;
  set_state_t ss2;
  mk_prn_set(&set, 3, prns2);
  mk_set_itr(&it, &ss, &set);
  mk_set_itr(&it2, &ss2, &set);
  mk_filter_itr(&fitr, &fs, &not_eq_ref, &tt, &it);
  fail_unless(ptr_itr_equality(&fitr, &it2));
}
END_TEST


START_TEST(test_set_ref)
{
  prn prns[] = {0,1,2};
  ptd_set_t ptd;
  mk_ptd_set(&ptd, 3, sizeof(prn), prns, 22, &prn_key);
  /* Test 1 */
  fail_unless(NULL == ptd.ref, "Ref must be contained in set.");

  /* Test 2 */
  prn ref = 1;
  set_ref(&ptd, ref);
  /* do pointer comparison */
  fail_unless(&prns[1] == ptd.ref);
  prn new_ref = *(prn *)ptd.ref;
  fail_unless(new_ref == ref, "Don't match: %i %i", new_ref, ref);
  /* Test 3 */
  iterator_t filter_itr;
  filter_state_t f_state;
  iterator_t set_itr;
  set_state_t s_state;
  mk_set_itr(&set_itr, &s_state, &ptd.set);
  mk_without_ref_itr(&filter_itr, &f_state, &set_itr, &ptd);
  each(&filter_itr, &print_prn, NULL); // 0, 2
}
END_TEST

// TODO make this a test, not a bunch of printf's
START_TEST(test_intersect)
{
  prn prns1[] = {0,1,4,7};
  prn prns2[] = {0,4,5,6,7};
  prn prns3[] = {};
  prn prns4[] = {3,10,11,12};
  prn prns5[] = {0,1,2,3,4,5,7};
  set_t set1, set2, set3, set4, set5;
  iterator_t itr1, itr2;
  set_state_t state1, state2;
  mk_prn_set(&set1, 4, prns1);
  mk_prn_set(&set2, 5, prns2);
  mk_prn_set(&set3, 0, prns3);
  mk_prn_set(&set4, 4, prns4);
  mk_prn_set(&set5, 7, prns5);

  /* Test 1 */
  mk_set_itr(&itr1, &state1, &set1);
  mk_set_itr(&itr2, &state2, &set2);
  iterator_t intersection;
  intersection_state_t s;
  prn curr;
  size_t prn_size = sizeof(prn);
  mk_intersection_itr(&intersection, &s, &curr, &itr1, &itr2, &prn_key, &prn_key, &fst, &prn_size);
  each(&intersection, &print_prn, NULL);
  /* Test 2 */
  mk_set_itr(&itr2, &state2, &set3);
  mk_intersection_itr(&intersection, &s, &curr, &itr1, &itr2, &prn_key, &prn_key, &snd, &prn_size);
  each(&intersection, &print_prn, NULL);
  /* Test 3 */
  mk_set_itr(&itr2, &state2, &set4);
  mk_intersection_itr(&intersection, &s, &curr, &itr1, &itr2, &prn_key, &prn_key, &snd, &prn_size);
  each(&intersection, &print_prn, NULL);

  /* Test 4 - Subset */
  mk_set_itr(&itr1, &state1, &set1);
  mk_set_itr(&itr2, &state2, &set5);
  fail_unless(is_subset(&itr1, &itr2, &prn_key, &prn_key));
  /* Test 5 */
  mk_set_itr(&itr2, &state2, &set4);
  fail_unless(!is_subset(&itr1, &itr2, &prn_key, &prn_key));
  mk_set_itr(&itr2, &state2, &set3);
  fail_unless(!is_subset(&itr1, &itr2, &prn_key, &prn_key));

}
END_TEST

key int_key(const void *num)
{
  return *(int *)num;
}
void fold_sum_int(const void *arg, void *sum, const void *elem)
{
  (void) arg;
  *(int *)sum = *(int *)sum + *(int *)elem;
}
START_TEST(test_fold)
{
  int nums[] = {0,1,2};
  set_t nums_set;
  iterator_t nums_itr;
  set_state_t nums_state;
  mk_set(&nums_set, 3, sizeof(int), nums, &int_key);
  mk_set_itr(&nums_itr, &nums_state, &nums_set);
  each(&nums_itr, &print_prn, NULL);
  int sum = 0;
  fold(&fold_sum_int, NULL, &sum, &nums_itr);
  printf("sum: %i\n", sum);
}
END_TEST

START_TEST(test_sats_man)
{
  u8 num_sats = 10;
  sdiff_t sdiffs1[num_sats], sdiffs2[num_sats];
  seed_rng();
  /* Test 1: choose_ref */
  for(prn prn = 0; prn < num_sats; prn++) {
    sdiffs1[prn].prn = prn;
    sdiffs1[prn].snr = frand(0,1);
    sdiffs2[prn].prn = prn;
    sdiffs2[prn].snr = frand(0,1);
  }
  //fail_unless(box_choose_reference_sat(num_sats, sdiffs1) ==
  //            old_choose_reference_sat(num_sats, sdiffs1));
  /* Test 2: intersect_sats */
  //prn other[3] = {0, 5, 15};
  //sdiff_t inter1[num_sats], inter2[num_sats];
  //u8 num_inter1 = intersect_sats( 3, num_sats, other, sdiffs1, inter1);
  //u8 num_inter2 = intersect_sats2(3, num_sats, other, sdiffs1, inter2);
  //fail_unless(num_inter1 == num_inter2, "Intersections should be the same size");
  //for (u8 i = 0; i < num_inter1; i++) {
  //  printf("%u, %u\n", inter1[i].prn, inter2[i].prn);
  //}
  //fail_unless(memcmp(inter1, inter2, num_inter1 * sizeof(sdiff_t)) == 0);
}
END_TEST
void (*each_body)(const void *arg, const void *elem);
#define each_block(itr) \
  iterator_t itr##_it; \
  set_state_t itr##_state;\
  
void do_each()
{
  prn ps[] = {22, 44};
  iterator_t it;
  set_state_t ss;
  set_t set;
  mk_set(&set, 2, sizeof(prn), ps, &prn_key);
  mk_set_itr(&it, &ss, &set);

  each_body = &print_prn;

  each(&it, each_body, NULL);
}
START_TEST(test_macros)
{
  do_each();
  free(NULL);
}
END_TEST

Suite* set_suite(void)
{
  Suite *s = suite_create("Set Utils");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_set);
  tcase_add_test(tc_core, test_freeze);
  tcase_add_test(tc_core, test_filter);
  tcase_add_test(tc_core, test_set_ref);
  tcase_add_test(tc_core, test_intersect);
  tcase_add_test(tc_core, test_fold);
  tcase_add_test(tc_core, test_sats_man);
  tcase_add_test(tc_core, test_macros);
  suite_add_tcase(s, tc_core);

  return s;
}

