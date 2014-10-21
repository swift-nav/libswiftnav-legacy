#include <stdio.h>
#include <check.h>

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
  freeze_itr(sizeof(prn), 22, frozen2, &it2);
  fail_unless(arr_eq(3, frozen2, prns2), "Frozen array does not match input to iterator.");

  fail_unless(freeze_itr(sizeof(prn), 2, frozen2, &it2) == -1, "Iterator exceeds max length?");
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
  set_t set;
  mk_prn_set(&set, 3, prns);
  /* Test 1 */
  fail_unless(-1 == set_ref_prn(&set, 3), "Ref must be contained in set.");

  /* Test 2 */
  prn ref = 1;
  fail_unless(0 == set_ref_prn(&set, ref));
  prn new_ref = *(prn *)set.ref;
  fail_unless(new_ref == ref, "Don't match: %i %i", new_ref, ref);
  /* Test 3 */
  iterator_t filter_itr;
  filter_state_t f_state;
  iterator_t set_itr;
  set_state_t s_state;
  mk_set_itr(&set_itr, &s_state, &set);
  mk_without_ref_itr(&filter_itr, &f_state, &set_itr, &set);
  each(&filter_itr, &print_prn, NULL); // 0, 2
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
  suite_add_tcase(s, tc_core);

  return s;
}

