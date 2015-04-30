#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <linear_algebra.h>
#include <ambiguity_test.h>
#include <printing_utils.h>

#include "check_utils.h"



/* Assure that when the sdiffs match amb_test's sats, amb_test's sats are unchanged. */
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

/* Assure that when we've lost the reference, we choose a new one and rebase everything. */
START_TEST(test_update_sats_rebase)
{
  ambiguity_test_t amb_test;
  create_empty_ambiguity_test(&amb_test);

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
  printf("N0: %i\n", hyp->N[0]);
  printf("N1: %i\n", hyp->N[1]);
  printf("N2: %i\n", hyp->N[2]);
  fail_unless(hyp->N[0] == -2);
  fail_unless(hyp->N[1] == -1);
}
END_TEST

START_TEST(test_ambiguity_update_reference)
{
  srandom(1);

  ambiguity_test_t amb_test = {.sats = {.num_sats = 4, 
                                        .prns = {3,1,2,4}}};
  create_empty_ambiguity_test(&amb_test);

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

void resize_matrix(u8 r1, u8 c1, u8 r2, u8 c2, const double *m1, double *m2)
{
  (void) r1;
  for (u8 i = 0; i < r2; i++) {
    for (u8 j = 0; j < c2; j++) {
      m2[i*c2 + j] = m1[i*c1 + j];
    }
  }
}

// Stephen K. Park and Keith W. Miller (1988)
int a = 16807;
int m = 2147483647;
unsigned int seed = 21474364;
float MINSTD()
{
  seed = (a * seed) % m;
  return (float)seed / (float)m;
}

START_TEST(test_amb_sat_inclusion)
{
  /* TODO Load matrix */
  u8 state_dim = 7;
  /* Test dimension */
  u8 dim = 7;
  double cov_mat[state_dim * state_dim];
  double multiplier[state_dim * state_dim];
  matrix_eye(state_dim, cov_mat);
  double diag = 0.08;

  for (u8 i = 0; i < state_dim; i++) {
    cov_mat[i*state_dim+i] = diag;
  }
  for (u8 i = 0; i < state_dim; i++) {
    for (u8 j = 0; j < state_dim; j++) {
      multiplier[i*state_dim+j] = MINSTD();
    }
  }
  for (u8 i = 0; i < state_dim; i++) {
    for (u8 j = 0; j < i; j++) {
      cov_mat[i*state_dim+j] = 0;
      cov_mat[j*state_dim+i] = 0;
    }
  }

  /* Compute a lightly randomized covariance matrix:
   *     multiplier * diag * multiplier^t
   */
  double multiplierT[state_dim * state_dim];
  double a[state_dim * state_dim];
  double b[state_dim * state_dim];
  matrix_transpose(state_dim, state_dim, multiplier, multiplierT);
  matrix_multiply(state_dim, state_dim, state_dim, multiplier, cov_mat, a);
  matrix_multiply(state_dim, state_dim, state_dim, a, multiplierT, b);
  matrix_copy(state_dim, state_dim, b, cov_mat);

  /* Take some block, factor */
  u8 prns[dim];
  for (u8 i = 0; i < dim+1; i++) {
    prns[i] = i;
  }
  double block[dim * dim];
  resize_matrix(state_dim, state_dim, dim, dim, cov_mat, block);

  printf("test covariance matrix:\n");
  print_double_mtx(block, dim, dim);

  double u[dim * dim];
  double d[dim * dim];
  matrix_udu(dim, block, u, d);
  double mean[dim];
  for (u8 i = 0; i < dim; i++) {
    mean[i] = 0;
  }

  /* Init amb_test */
  ambiguity_test_t amb_test;
  create_ambiguity_test(&amb_test);
  sats_management_t float_sats = {
    .num_sats = dim+1,
  };
  memcpy(float_sats.prns, prns, (dim+1) * sizeof(u8));

  u16 pool_size;
  u8 flag;
  pool_size = memory_pool_n_allocated(amb_test.pool);
  printf("pool size before: %i\n", pool_size);
  /* Include. This one should succeed and add 5 sats. */
  flag = ambiguity_sat_inclusion(&amb_test, 0, &float_sats, mean, u, d);
  printf("inclusion return code: %i\n", flag);
  pool_size = memory_pool_n_allocated(amb_test.pool);
  printf("pool size after 1: %i\n", pool_size);
  fail_unless(flag == 1);
  fail_unless(pool_size == 625);
  /* Include again. This one should succeed and add 1 more sat. */
  flag = ambiguity_sat_inclusion(&amb_test, 0, &float_sats, mean, u, d);
  printf("inclusion return code: %i\n", flag);
  pool_size = memory_pool_n_allocated(amb_test.pool);
  printf("pool size after 2: %i\n", pool_size);
  fail_unless(flag == 1);
  fail_unless(pool_size == 945);
  /* Include again. This one should fail. */
  flag = ambiguity_sat_inclusion(&amb_test, 0, &float_sats, mean, u, d);
  printf("inclusion return code: %i\n", flag);
  pool_size = memory_pool_n_allocated(amb_test.pool);
  printf("pool size after 2: %i\n", pool_size);
  fail_unless(flag == 0);
  fail_unless(pool_size == 945);
}
END_TEST


Suite* ambiguity_test_suite(void)
{
  Suite *s = suite_create("Ambiguity Test");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_sats_match);
  tcase_add_test(tc_core, test_ambiguity_update_reference);
  tcase_add_test(tc_core, test_update_sats_same_sats);
  // TODO add back
  //tcase_add_test(tc_core, test_update_sats_rebase);
  (void) test_update_sats_rebase;
  tcase_add_test(tc_core, test_amb_sat_inclusion);
  suite_add_tcase(s, tc_core);

  return s;
}

