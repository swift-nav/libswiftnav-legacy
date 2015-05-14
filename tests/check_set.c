
#include <check.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <set.h>

#include "check_utils.h"

#define LEN(x) (sizeof(x) / sizeof(x[0]))

int cmp_s32(const void * a, const void * b)
{
   return *(s32*)a - *(s32*)b;
}

void test_map_f(void *context, u32 n, const void *a_, const void *b_)
{
  s32 *a = (s32 *)a_;
  s32 *b = (s32 *)b_;
  s32 *c = (s32 *)context;

  fail_unless(context != NULL);
  fail_unless(*a == *b, "Intersection items not equal");
  c[n] = *a;
}

#define TEST_INTERSECTION(a, b, c) {                         \
  s32 c_result[LEN(c)];                                      \
                                                             \
  qsort(a, LEN(a), sizeof(a[0]), cmp_s32);                   \
  qsort(b, LEN(b), sizeof(b[0]), cmp_s32);                   \
  qsort(c, LEN(c), sizeof(c[0]), cmp_s32);                   \
                                                             \
  u32 ret = intersection_map(LEN(a), sizeof(a[0]), a,        \
                             LEN(b), sizeof(b[0]), b,        \
                             cmp_s32, c_result, test_map_f); \
                                                             \
  fail_unless(ret == LEN(c),                                 \
      "Intersection length does not match test data");       \
                                                             \
  fail_unless(memcmp(c, c_result, sizeof(c)) == 0,           \
      "Output of intersection does not match test data");    \
}

START_TEST(test_intersection_map_1)
{
  /* Empty first set */
  s32 a[] = {};
  s32 b[] = {1, 2, 3};
  s32 c[] = {};
  TEST_INTERSECTION(a, b, c)
}
END_TEST

START_TEST(test_intersection_map_2)
{
  /* Empty second set */
  s32 a[] = {1, 2, 3};
  s32 b[] = {};
  s32 c[] = {};
  TEST_INTERSECTION(a, b, c)
}
END_TEST

START_TEST(test_intersection_map_3)
{
  /* Beginning intersects */
  s32 a[] = {1, 2, 3, 4, 5, 6, 7};
  s32 b[] = {1, 2, 3};
  s32 c[] = {1, 2, 3};
  TEST_INTERSECTION(a, b, c)
}
END_TEST

START_TEST(test_intersection_map_4)
{
  /* End intersects */
  s32 a[] = {1, 2, 3, 4, 5, 6, 7};
  s32 b[] = {5, 6, 7};
  s32 c[] = {5, 6, 7};
  TEST_INTERSECTION(a, b, c)
}
END_TEST

START_TEST(test_intersection_map_5)
{
  /* Same set */
  s32 a[] = {1, 2, 3, 4, 5, 6, 7};
  s32 b[] = {1, 2, 3, 4, 5, 6, 7};
  s32 c[] = {1, 2, 3, 4, 5, 6, 7};
  TEST_INTERSECTION(a, b, c)
}
END_TEST

START_TEST(test_intersection_map_6)
{
  /* Disjoint */
  s32 a[] = {1, 2, 3, 4};
  s32 b[] = {5, 6, 7, 8};
  s32 c[] = {};
  TEST_INTERSECTION(a, b, c)
}
END_TEST

START_TEST(test_intersection_map_7)
{
  /* Middle overlaps */
  s32 a[] = {1, 2, 3, 4, 5, 6, 7};
  s32 b[] = {5, 6};
  s32 c[] = {5, 6};
  TEST_INTERSECTION(a, b, c)
}
END_TEST

START_TEST(test_intersection_map_8)
{
  /* Overlapping but not subset */
  s32 a[] = {1, 2, 3, 4, 5};
  s32 b[] = {4, 5, 6, 7, 8};
  s32 c[] = {4, 5};
  TEST_INTERSECTION(a, b, c)
}
END_TEST

START_TEST(test_intersection_map_9)
{
  /* Alternating disjoint */
  s32 a[] = {2, 4, 6, 8, 10};
  s32 b[] = {1, 3, 7, 9, 11};
  s32 c[] = {};
  TEST_INTERSECTION(a, b, c)
}
END_TEST

START_TEST(test_intersection_map_10)
{
  /* Alternating with overlap */
  s32 a[] = {2, 4, 6, 8, 9, 10};
  s32 b[] = {1, 3, 7, 8, 9, 11};
  s32 c[] = {8, 9};
  TEST_INTERSECTION(a, b, c)
}
END_TEST

Suite* set_suite(void)
{
  Suite *s = suite_create("Set");

  TCase *tc_intersection = tcase_create("Intersection");
  tcase_add_test(tc_intersection, test_intersection_map_1);
  tcase_add_test(tc_intersection, test_intersection_map_2);
  tcase_add_test(tc_intersection, test_intersection_map_3);
  tcase_add_test(tc_intersection, test_intersection_map_4);
  tcase_add_test(tc_intersection, test_intersection_map_5);
  tcase_add_test(tc_intersection, test_intersection_map_6);
  tcase_add_test(tc_intersection, test_intersection_map_7);
  tcase_add_test(tc_intersection, test_intersection_map_8);
  tcase_add_test(tc_intersection, test_intersection_map_9);
  tcase_add_test(tc_intersection, test_intersection_map_10);
  suite_add_tcase(s, tc_intersection);

  return s;
}

