
#include <check.h>
#include <stdio.h>

#include <memory_pool.h>

#include "check_utils.h"

memory_pool_t *test_pool_seq;
memory_pool_t *test_pool_random;
memory_pool_t *test_pool_empty;

void setup()
{
  /* Create a new pool and fill it with a sequence of ints. */
  test_pool_seq = memory_pool_new(50, sizeof(s32));

  s32 *x;
  for (u32 i=0; i<22; i++) {
    x = (s32 *)memory_pool_add(test_pool_seq);
    fail_unless(x != 0, "Null pointer returned by memory_pool_add");
    *x = i;
  }
  /* Create a new pool and fill it entirely with random numbers. */
  test_pool_random = memory_pool_new(20, sizeof(s32));

  for (u32 i=0; i<20; i++) {
    x = (s32 *)memory_pool_add(test_pool_random);
    fail_unless(x != 0, "Null pointer returned by memory_pool_add");
    *x = sizerand(16);
  }

  /* Create a new pool and leave it empty. */
  test_pool_empty = memory_pool_new(50, sizeof(s32));
}

void teardown()
{
  memory_pool_destroy(test_pool_seq);
  memory_pool_destroy(test_pool_random);
}

void print_elem(element_t *elem)
{
  printf("%d ", *((s32 *)elem));
}

s32 isum(s32 x, element_t *elem)
{
  s32 *y = (s32 *)elem;
  return x + *y;
}

double dsum(double x, element_t *elem)
{
  s32 *y = (s32 *)elem;
  return x + *y;
}

float fsum(float x, element_t *elem)
{
  s32 *y = (s32 *)elem;
  return x + *y;
}

START_TEST(test_simple_folds)
{
  s32 sum;

  sum = memory_pool_ifold(test_pool_seq, 0, &isum);
  fail_unless(sum == 231,
      "Fold failed for isum function, expected 231, got %d", sum);

  sum = memory_pool_ifold(test_pool_seq, 22, &isum);
  fail_unless(sum == 253,
      "Fold failed for isum function, expected 253, got %d", sum);

  sum = (s32)memory_pool_dfold(test_pool_seq, 22, &dsum);
  fail_unless(sum == 253,
      "Fold failed for dsum function, expected 253, got %d", sum);

  sum = (s32)memory_pool_ffold(test_pool_seq, 22, &fsum);
  fail_unless(sum == 253,
      "Fold failed for fsum function, expected 253, got %d", sum);
}
END_TEST

typedef struct {
  s32 min;
  s32 max;
} min_max_t;

void min_max_finder(void *x, element_t *elem)
{
  s32 *y = (s32 *)elem;
  min_max_t *mm = (min_max_t *)x;

  if (*y > mm->max)
    mm->max = *y;

  if (*y < mm->min)
    mm->min = *y;
}

START_TEST(test_general_fold)
{
  min_max_t mm = { .min = 1000, .max = -1000 };
  memory_pool_fold(test_pool_seq, &mm, &min_max_finder);

  fail_unless(mm.min == 0,
      "Fold failed for min_max_finder function, expected min 0, got %d", mm.min);

  fail_unless(mm.max == 21,
      "Fold failed for min_max_finder function, expected min 21, got %d", mm.max);
}
END_TEST

START_TEST(test_full)
{
  fail_unless(memory_pool_add(test_pool_random) == 0,
      "Adding to a full pool should return a NULL pointer.");
}
END_TEST

START_TEST(test_n_free)
{
  s32 n = memory_pool_n_free(test_pool_random);
  fail_unless(n == 0,
      "Error calculating free space in pool, expected 0, got %d", n);

  n = memory_pool_n_free(test_pool_seq);
  fail_unless(n == 28,
      "Error calculating free space in pool, expected 28, got %d", n);

  n = memory_pool_n_free(test_pool_empty);
  fail_unless(n == 50,
      "Error calculating free space in pool, expected 50, got %d", n);

  s32 i = 49;
  while(memory_pool_add(test_pool_empty)) {
    n = memory_pool_n_free(test_pool_empty);
    fail_unless(n == i,
        "Error calculating free space in pool, expected %d, got %d", i, n);
    i--;
  }
}
END_TEST

START_TEST(test_n_allocated)
{
  s32 n = memory_pool_n_allocated(test_pool_random);
  fail_unless(n == 20,
      "Error calculating free space in pool, expected 20, got %d", n);

  n = memory_pool_n_allocated(test_pool_seq);
  fail_unless(n == 22,
      "Error calculating free space in pool, expected 22, got %d", n);

  n = memory_pool_n_allocated(test_pool_empty);
  fail_unless(n == 0,
      "Error calculating free space in pool, expected 0, got %d", n);

  s32 i = 1;
  while(memory_pool_add(test_pool_empty)) {
    n = memory_pool_n_allocated(test_pool_empty);
    fail_unless(n == i,
        "Error calculating free space in pool, expected %d, got %d", i, n);
    i++;
  }
}
END_TEST

START_TEST(test_pool_to_array)
{
  s32 xs[22];
  s32 test_xs[22] = {
    21, 20, 19, 18, 17, 16, 15, 14, 13, 12,
    11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
  };
  memory_pool_to_array(test_pool_seq, xs);

  fail_unless(memcmp(xs, test_xs, sizeof(xs)) == 0,
      "Output of memory_pool_to_array does not match test data");
}
END_TEST

void times_two(element_t *elem)
{
  s32 *x = (s32 *)elem;
  *x = *x * 2;
}

START_TEST(test_map)
{
  s32 xs[22];
  s32 test_xs[22] = {
    42, 40, 38, 36, 34, 32, 30, 28, 26, 24,
    22, 20, 18, 16, 14, 12, 10, 8, 6, 4, 2, 0
  };

  memory_pool_map(test_pool_seq, &times_two);

  memory_pool_to_array(test_pool_seq, xs);
  fail_unless(memcmp(xs, test_xs, sizeof(xs)) == 0,
      "Output of map operation does not match test data");
}
END_TEST

s8 less_than_12(element_t *elem)
{
  return *((s32 *)elem) < 12;
}

START_TEST(test_filter_1)
{
  s32 xs[22];
  s32 test_xs_beginning[12] = {
    11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
  };

  memory_pool_filter(test_pool_seq, &less_than_12);
  fail_unless(memory_pool_n_allocated(test_pool_seq) == 12,
      "Filtered length does not match");

  memory_pool_to_array(test_pool_seq, xs);
  fail_unless(memcmp(xs, test_xs_beginning, sizeof(test_xs_beginning)) == 0,
      "Output of filter operation does not match test data");
}
END_TEST

s8 between_15_7(element_t *elem)
{
  return *((s32 *)elem) <= 15 && *((s32 *)elem) >= 7;
}

START_TEST(test_filter_2)
{
  s32 xs[22];
  s32 test_xs_middle[9] = {
    15, 14, 13, 12, 11, 10, 9, 8, 7
  };

  memory_pool_filter(test_pool_seq, &between_15_7);
  fail_unless(memory_pool_n_allocated(test_pool_seq) == 9,
      "Filtered length does not match");

  memory_pool_to_array(test_pool_seq, xs);
  fail_unless(memcmp(xs, test_xs_middle, sizeof(test_xs_middle)) == 0,
      "Output of filter operation does not match test data");
}
END_TEST

s8 greater_than_11(element_t *elem)
{
  return *((s32 *)elem) > 11;
}

START_TEST(test_filter_3)
{
  s32 xs[22];
  s32 test_xs_end[10] = {
    21, 20, 19, 18, 17, 16, 15, 14, 13, 12
  };

  memory_pool_filter(test_pool_seq, &greater_than_11);
  fail_unless(memory_pool_n_allocated(test_pool_seq) == 10,
      "Filtered length does not match");

  memory_pool_to_array(test_pool_seq, xs);
  fail_unless(memcmp(xs, test_xs_end, sizeof(test_xs_end)) == 0,
      "Output of filter operation does not match test data");
}
END_TEST

s8 even(element_t *elem)
{
  return *((s32 *)elem) % 2 == 0;
}

START_TEST(test_filter_4)
{
  s32 xs[22];
  s32 test_xs_evens[11] = {
    20, 18, 16, 14, 12, 10, 8, 6, 4, 2, 0
  };

  memory_pool_filter(test_pool_seq, &even);
  fail_unless(memory_pool_n_allocated(test_pool_seq) == 11,
      "Filtered length does not match");

  memory_pool_to_array(test_pool_seq, xs);
  fail_unless(memcmp(xs, test_xs_evens, sizeof(test_xs_evens)) == 0,
      "Output of filter operation does not match test data");
}
END_TEST

Suite* memory_pool_suite(void)
{
  Suite *s = suite_create("Memory Pools");

  TCase *tc_core = tcase_create("Core");
  tcase_add_checked_fixture (tc_core, setup, teardown);
  tcase_add_test(tc_core, test_simple_folds);
  tcase_add_test(tc_core, test_general_fold);
  tcase_add_test(tc_core, test_full);
  tcase_add_test(tc_core, test_n_free);
  tcase_add_test(tc_core, test_n_allocated);
  tcase_add_test(tc_core, test_pool_to_array);
  tcase_add_test(tc_core, test_map);
  tcase_add_test(tc_core, test_filter_1);
  tcase_add_test(tc_core, test_filter_2);
  tcase_add_test(tc_core, test_filter_3);
  tcase_add_test(tc_core, test_filter_4);
  suite_add_tcase(s, tc_core);

  return s;
}

