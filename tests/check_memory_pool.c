
#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <memory_pool.h>

#include "check_utils.h"

memory_pool_t *test_pool_seq;
memory_pool_t *test_pool_random;
memory_pool_t *test_pool_empty;

void setup()
{
  /* Seed the random number generator with a specific seed for our test. */
  srandom(1);

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
  /* Test pools for leaks. */
  fail_unless(memory_pool_n_free(test_pool_seq)
              + memory_pool_n_allocated(test_pool_seq)
              == 50,
              "Memory leak! test_pool_seq lost elements!");

  fail_unless(memory_pool_n_free(test_pool_random)
              + memory_pool_n_allocated(test_pool_random)
              == 20,
              "Memory leak! test_pool_random lost elements!");

  fail_unless(memory_pool_n_free(test_pool_empty)
              + memory_pool_n_allocated(test_pool_empty)
              == 50,
              "Memory leak! test_pool_empty lost elements!");

  memory_pool_destroy(test_pool_seq);
  memory_pool_destroy(test_pool_random);
  memory_pool_destroy(test_pool_empty);
}

void print_s32(element_t *elem)
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

START_TEST(test_empty)
{
  fail_unless(!memory_pool_empty(test_pool_random),
      "Error checking if memory pool empty, should have been non-empty");

  fail_unless(memory_pool_empty(test_pool_empty),
      "Error checking if memory pool empty, should have been empty");
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

void times_two(void *arg, element_t *elem)
{
  (void) arg;
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

  memory_pool_map(test_pool_seq, NULL, &times_two);

  memory_pool_to_array(test_pool_seq, xs);
  fail_unless(memcmp(xs, test_xs, sizeof(xs)) == 0,
      "Output of map operation does not match test data");
}
END_TEST

s8 less_than_12(void *arg, element_t *elem)
{
  (void) arg;
  return *((s32 *)elem) < 12;
}

START_TEST(test_filter_1)
{
  s32 xs[22];
  s32 test_xs_beginning[12] = {
    11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
  };

  memory_pool_filter(test_pool_seq, NULL, &less_than_12);
  fail_unless(memory_pool_n_allocated(test_pool_seq) == 12,
      "Filtered length does not match");

  memory_pool_to_array(test_pool_seq, xs);
  fail_unless(memcmp(xs, test_xs_beginning, sizeof(test_xs_beginning)) == 0,
      "Output of filter operation does not match test data");
}
END_TEST

s8 between_15_7(void *arg, element_t *elem)
{
  (void) arg;
  return *((s32 *)elem) <= 15 && *((s32 *)elem) >= 7;
}

START_TEST(test_filter_2)
{
  s32 xs[22];
  s32 test_xs_middle[9] = {
    15, 14, 13, 12, 11, 10, 9, 8, 7
  };

  memory_pool_filter(test_pool_seq, NULL, &between_15_7);
  fail_unless(memory_pool_n_allocated(test_pool_seq) == 9,
      "Filtered length does not match");

  memory_pool_to_array(test_pool_seq, xs);
  fail_unless(memcmp(xs, test_xs_middle, sizeof(test_xs_middle)) == 0,
      "Output of filter operation does not match test data");
}
END_TEST

s8 greater_than_11(void *arg, element_t *elem)
{
  (void) arg;
  return *((s32 *)elem) > 11;
}

START_TEST(test_filter_3)
{
  s32 xs[22];
  s32 test_xs_end[10] = {
    21, 20, 19, 18, 17, 16, 15, 14, 13, 12
  };

  memory_pool_filter(test_pool_seq, NULL, &greater_than_11);
  fail_unless(memory_pool_n_allocated(test_pool_seq) == 10,
      "Filtered length does not match");

  memory_pool_to_array(test_pool_seq, xs);
  fail_unless(memcmp(xs, test_xs_end, sizeof(test_xs_end)) == 0,
      "Output of filter operation does not match test data");
}
END_TEST

s8 even(void *arg, element_t *elem)
{
  (void) arg;
  return *((s32 *)elem) % 2 == 0;
}

START_TEST(test_filter_4)
{
  s32 xs[22];
  s32 test_xs_evens[11] = {
    20, 18, 16, 14, 12, 10, 8, 6, 4, 2, 0
  };

  memory_pool_filter(test_pool_seq, NULL, &even);
  fail_unless(memory_pool_n_allocated(test_pool_seq) == 11,
      "Filtered length does not match");

  memory_pool_to_array(test_pool_seq, xs);
  fail_unless(memcmp(xs, test_xs_evens, sizeof(test_xs_evens)) == 0,
      "Output of filter operation does not match test data");
}
END_TEST

s32 cmp_s32s(void *arg, element_t *a_, element_t *b_)
{
  (void)arg;
  s32 *a = (s32 *)a_;
  s32 *b = (s32 *)b_;

  return *a - *b;
}

START_TEST(test_clear)
{
  memory_pool_clear(test_pool_seq);
  fail_unless(memory_pool_n_allocated(test_pool_seq) == 0,
      "Pool still not empty after clear");

  memory_pool_clear(test_pool_random);
  fail_unless(memory_pool_n_allocated(test_pool_random) == 0,
      "Pool still not empty after clear");

  memory_pool_clear(test_pool_empty);
  fail_unless(memory_pool_n_allocated(test_pool_empty) == 0,
      "Pool still not empty after clear");
}
END_TEST

START_TEST(test_sort)
{
  s32 xs[22];
  s32 test_xs_sorted[22] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
    12, 13, 14, 15, 16, 17, 18, 19, 20, 21
  };

  memory_pool_sort(test_pool_seq, 0, &cmp_s32s);
  fail_unless(memory_pool_n_allocated(test_pool_seq) == 22,
      "Sorted length does not match");

  memory_pool_to_array(test_pool_seq, xs);
  fail_unless(memcmp(xs, test_xs_sorted, sizeof(test_xs_sorted)) == 0,
      "Output of sort operation does not match test data");

  s32 ys[20];
  s32 test_ys_sorted[20] = {
    3, 4, 5, 6, 6, 7, 8, 9, 9, 10, 11,
    11, 12, 13, 13, 13, 14, 15, 15, 16
  };

  memory_pool_sort(test_pool_random, 0, &cmp_s32s);
  fail_unless(memory_pool_n_allocated(test_pool_random) == 20,
      "Sorted length does not match");

  memory_pool_to_array(test_pool_random, ys);
  fail_unless(memcmp(ys, test_ys_sorted, sizeof(test_ys_sorted)) == 0,
      "Output of sort operation does not match test data");

  /* Test sorting an empty list. */
  memory_pool_sort(test_pool_empty, 0, &cmp_s32s);
  fail_unless(memory_pool_n_allocated(test_pool_empty) == 0,
      "Sorted length does not match");
  fail_unless(memory_pool_n_free(test_pool_empty) == 50,
      "Sorted length does not match");

  /* Test sorting a single element list. */
  s32 *x = (s32 *)memory_pool_add(test_pool_empty); *x = 22;

  memory_pool_sort(test_pool_empty, 0, &cmp_s32s);
  fail_unless(memory_pool_n_allocated(test_pool_empty) == 1,
      "Sorted length does not match");
  fail_unless(memory_pool_n_free(test_pool_empty) == 49,
      "Sorted length does not match");
}
END_TEST

s32 group_evens(void *arg, element_t *a_, element_t *b_)
{
  (void)arg;
  s32 *a = (s32 *)a_;
  s32 *b = (s32 *)b_;

  return (*a % 2) - (*b % 2);
}

void agg_sum_s32s(element_t *new_, void *x, u32 n, element_t *elem_)
{
  (void)x;

  s32 *new = (s32 *)new_;
  s32 *elem = (s32 *)elem_;

  if (n == 0)
    *new = 0;

  *new += *elem;
}

START_TEST(test_groupby_1)
{
  s32 xs[2];
  s32 test_xs_reduced[2] = {121, 110};

  memory_pool_group_by(test_pool_seq, 0, &group_evens, 0, 0, &agg_sum_s32s);

  fail_unless(memory_pool_n_allocated(test_pool_seq) == 2,
      "Reduced length does not match");

  memory_pool_to_array(test_pool_seq, xs);
  fail_unless(memcmp(xs, test_xs_reduced, sizeof(test_xs_reduced)) == 0,
      "Output of groupby operation does not match test data");

  /*memory_pool_map(test_pool_seq, &print_s32); printf("\n");*/
}
END_TEST

typedef struct __attribute__((packed)) {
  float p;
  u8 len;
  u8 N[15];
} hypothesis_t;

void print_hyp_elem(element_t *hyp_)
{
  hypothesis_t *hyp = (hypothesis_t *)hyp_;
  printf("[ ");
  for (u8 i=0; i<hyp->len; i++)
    printf("%d ", hyp->N[i]);
  printf("] %.3f\n", hyp->p);
}

s32 group_by_N_i(void *i_, element_t *a_, element_t *b_)
{
  u8 *i = (u8 *)i_;
  hypothesis_t *a = (hypothesis_t *)a_;
  hypothesis_t *b = (hypothesis_t *)b_;

  for (u8 n=0; n<a->len; n++) {
    if (n == *i)
      continue;
    if (a->N[n] < b->N[n])
      return -1;
    if (a->N[n] > b->N[n])
      return 1;
  }
  return 0;
}

void agg_sum_p(element_t *new_, void *x_, u32 n, element_t *elem_)
{
  u8 *x = (u8 *)x_;
  hypothesis_t *new = (hypothesis_t *)new_;
  hypothesis_t *elem = (hypothesis_t *)elem_;

  if (n == 0) {
    new->len = elem->len - 1;
    new->p = 0;
    u8 j=0;
    for (u8 i=0; i<elem->len; i++) {
      if (i != *x) {
        new->N[j] = elem->N[i];
        j++;
      }
    }
  }
  new->p += elem->p;
}

START_TEST(test_groupby_2)
{
  /* Create a new pool. */
  memory_pool_t *test_pool_hyps = memory_pool_new(50, sizeof(hypothesis_t));

  for (u32 i=0; i<10; i++) {
    hypothesis_t *hyp = (hypothesis_t *)memory_pool_add(test_pool_hyps);
    fail_unless(hyp != 0, "Null pointer returned by memory_pool_add");
    hyp->len = 4;
    for (u8 j=0; j<hyp->len; j++) {
      hyp->N[j] = sizerand(2);
    }
    hyp->p = frand(0, 1);
  }

  u8 col = 1;

  /*memory_pool_map(test_pool_hyps, &print_hyp_elem); printf("\n");*/
  /*memory_pool_sort(test_pool_hyps, &col, &group_by_N_i);*/
  /*memory_pool_map(test_pool_hyps, &print_hyp_elem); printf("\n");*/

  memory_pool_group_by(test_pool_hyps, &col, &group_by_N_i,
                       &col, 1, &agg_sum_p);

  /*memory_pool_map(test_pool_hyps, &print_hyp_elem); printf("\n");*/

  fail_unless(memory_pool_n_allocated(test_pool_hyps) == 7,
      "Reduced length does not match");

  hypothesis_t reduced_hyps[7];
  memory_pool_to_array(test_pool_hyps, reduced_hyps);

  for (u32 i=0; i<7; i++) {
    fail_unless(reduced_hyps[i].len == 3,
        "Output of groupby operation does not match test data");
    fail_unless(reduced_hyps[i].p >= 0 && reduced_hyps[i].p <= 2,
        "Output of groupby operation does not match test data, p not in range: %f",
        reduced_hyps[i].p);
  }

  fail_unless(fabs(reduced_hyps[5].p - 0.374936) < 0.00001,
      "Output of groupby operation does not match test data, p[5] = %f, expected 0.375",
      reduced_hyps[5].p);

  /* Test pool for leaks. */
  fail_unless(memory_pool_n_free(test_pool_hyps)
              + memory_pool_n_allocated(test_pool_hyps)
              == 50,
              "Memory leak! test_pool_hyps lost elements!");

  memory_pool_destroy(test_pool_hyps);
}
END_TEST

void prod_N(element_t *new_, void *x_, u32 n_xs, u32 n, element_t *elem_)
{
  (void)n;
  u8 *x = (u8 *)x_;
  hypothesis_t *new = (hypothesis_t *)new_;
  hypothesis_t *elem = (hypothesis_t *)elem_;

  new->len = elem->len + 1;
  new->p = elem->p / n_xs;
  new->N[new->len-1] = *x;
}

START_TEST(test_prod)
{
  /* Create a new pool. */
  memory_pool_t *test_pool_hyps = memory_pool_new(50, sizeof(hypothesis_t));

  for (u32 i=0; i<3; i++) {
    hypothesis_t *hyp = (hypothesis_t *)memory_pool_add(test_pool_hyps);
    fail_unless(hyp != 0, "Null pointer returned by memory_pool_add");
    hyp->len = 4;
    for (u8 j=0; j<hyp->len; j++) {
      hyp->N[j] = sizerand(5);
    }
    hyp->p = frand(0, 1);
  }

  u8 new_N_vals[3] = {7, 8, 9};

  /*memory_pool_map(test_pool_hyps, &print_hyp_elem); printf("\n");*/

  memory_pool_product(test_pool_hyps, new_N_vals, 3, sizeof(u8), &prod_N);

  /*memory_pool_map(test_pool_hyps, &print_hyp_elem); printf("\n");*/

  fail_unless(memory_pool_n_allocated(test_pool_hyps) == 9,
      "Reduced length does not match");

  u8 col = 4;
  memory_pool_group_by(test_pool_hyps, &col, &group_by_N_i,
                       &col, 1, &agg_sum_p);

  fail_unless(memory_pool_n_allocated(test_pool_hyps) == 3,
      "Reduced length does not match");

  /*memory_pool_map(test_pool_hyps, &print_hyp_elem); printf("\n");*/

  memory_pool_destroy(test_pool_hyps);
}
END_TEST

typedef struct {
  u8 val;
  u8 n_vals;
  u8 i;
} test_gen_state_t;

s8 test_next(void *x_, u32 n)
{
  (void) n;
  test_gen_state_t *x = (test_gen_state_t *)x_;
  x->val = x->i * x->i;
  x->i++;
  if (x->i >= x->n_vals)
    return 0;
  return 1;
}

void prod_N_gen(element_t *new_, void *x_, u32 n, element_t *elem_)
{
  (void)n;
  test_gen_state_t *x = (test_gen_state_t *)x_;
  hypothesis_t *new = (hypothesis_t *)new_;
  hypothesis_t *elem = (hypothesis_t *)elem_;

  new->len = elem->len + 1;
  new->N[new->len-1] = x->val;
}

START_TEST(test_prod_generator)
{
  /* Create a new pool. */
  memory_pool_t *test_pool_hyps = memory_pool_new(50, sizeof(hypothesis_t));

  for (u32 i=0; i<3; i++) {
    hypothesis_t *hyp = (hypothesis_t *)memory_pool_add(test_pool_hyps);
    fail_unless(hyp != 0, "Null pointer returned by memory_pool_add");
    hyp->len = 4;
    for (u8 j=0; j<hyp->len; j++) {
      hyp->N[j] = sizerand(5);
    }
    hyp->p = frand(0, 1);
  }

  test_gen_state_t test_gen_state = {
    .val = 0,
    .n_vals = 3,
    .i = 0
  };

  memory_pool_product_generator(test_pool_hyps, &test_gen_state, 100, sizeof(test_gen_state_t), &test_next, &prod_N_gen);

  fail_unless(memory_pool_n_allocated(test_pool_hyps) == 9,
      "Reduced length does not match22");

  u8 col = 4;
  memory_pool_group_by(test_pool_hyps, &col, &group_by_N_i,
                       &col, 1, &agg_sum_p);

  fail_unless(memory_pool_n_allocated(test_pool_hyps) == 3,
      "Reduced length does not match33");

  memory_pool_destroy(test_pool_hyps);
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
  tcase_add_test(tc_core, test_empty);
  tcase_add_test(tc_core, test_pool_to_array);
  tcase_add_test(tc_core, test_map);
  tcase_add_test(tc_core, test_filter_1);
  tcase_add_test(tc_core, test_filter_2);
  tcase_add_test(tc_core, test_filter_3);
  tcase_add_test(tc_core, test_filter_4);
  tcase_add_test(tc_core, test_clear);
  tcase_add_test(tc_core, test_sort);
  tcase_add_test(tc_core, test_groupby_1);
  tcase_add_test(tc_core, test_groupby_2);
  tcase_add_test(tc_core, test_prod);
  tcase_add_test(tc_core, test_prod_generator);
  suite_add_tcase(s, tc_core);

  return s;
}

