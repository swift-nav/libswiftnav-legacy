#include <check.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include <linear_algebra.h>
#include <ambiguity_test.h>
#include <printing_utils.h>
#include <lambda.h>

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

  ambiguity_update_sats(&amb_test, num_sdiffs, sdiffs, NULL, NULL, NULL, NULL, false);

  fail_unless(amb_test.sats.prns[0] == 3);
  fail_unless(amb_test.sats.prns[1] == 1);
  fail_unless(amb_test.sats.prns[2] == 2);
  fail_unless(amb_test.sats.prns[3] == 4);
}
END_TEST

/* Assure that when the measurement is bad, we drop (but don't add) sats */
START_TEST(test_bad_measurements)
{
  ambiguity_test_t amb_test;

  sdiff_t sdiffs[5] = {{.prn = 1},
                       {.prn = 2},
                       {.prn = 3},
                       {.prn = 5},
                       {.prn = 6}};
  u8 num_sdiffs = 5;


  sats_management_t float_sats = {.num_sats = 5,
                                  .prns = {3, 1, 2, 5, 6}};
  double U[16];
  matrix_eye(4, U);
  double D[4] = {1, 1, 1, 1};
  double est[5] = {1, 2, 5, 6};

  sats_management_t amb_sats_init = {.num_sats = 5,
                                     .prns = {3, 1, 2, 4, 5}};
  hypothesis_t hyp_init = {.N = {1,2,4,5}};

  create_empty_ambiguity_test(&amb_test);
  memcpy(&amb_test.sats, &amb_sats_init, sizeof(sats_management_t));
  hypothesis_t *hyp = (hypothesis_t *)memory_pool_add(amb_test.pool);
  memcpy(hyp, &hyp_init, sizeof(hypothesis_t));
  /* Test that with a good measurement, we get a projection and inclusion.
   * It should have dropped PRN 4 and include PRN 6. */
  ambiguity_update_sats(&amb_test, num_sdiffs, sdiffs, &float_sats, est, U, D, false);
  fail_unless(amb_test.sats.num_sats == 5);
  fail_unless(amb_test.sats.prns[0] == 3);
  fail_unless(amb_test.sats.prns[1] == 1);
  fail_unless(amb_test.sats.prns[2] == 2);
  fail_unless(amb_test.sats.prns[3] == 5);
  fail_unless(amb_test.sats.prns[4] == 6);
  /* Reset the amb_test to what it was before ambiguity_update_sats */
  create_empty_ambiguity_test(&amb_test);
  memcpy(&amb_test.sats, &amb_sats_init, sizeof(sats_management_t));
  hyp = (hypothesis_t *)memory_pool_add(amb_test.pool);
  memcpy(hyp, &hyp_init, sizeof(hypothesis_t));
  /* Test that with a bad measurement, we get (only) a projection.
   * It should have dropped PRN 4 and NOT include PRN 6. */
  ambiguity_update_sats(&amb_test, num_sdiffs, sdiffs, &float_sats, est, U, D, true);
  fail_unless(amb_test.sats.num_sats == 4);
  fail_unless(amb_test.sats.prns[0] == 3);
  fail_unless(amb_test.sats.prns[1] == 1);
  fail_unless(amb_test.sats.prns[2] == 2);
  fail_unless(amb_test.sats.prns[3] == 5);
  /* And it should have dropped PRN's value from the hypothesis,
   * which should still be there and still be the only one. */
  hyp = (hypothesis_t *) amb_test.pool->allocated_nodes_head->elem;
  fail_unless(ambiguity_test_n_hypotheses(&amb_test) == 1);
  fail_unless(hyp->N[0] == 1);
  fail_unless(hyp->N[1] == 2);
  fail_unless(hyp->N[2] == 5);
}
END_TEST

static void amb_test_setup(ambiguity_test_t *amb_test)
{
  create_empty_ambiguity_test(amb_test);
  memcpy(&amb_test->sats, &((sats_management_t) {.num_sats = 5,
                                                .prns = {3, 1, 2, 4, 5}}),
    sizeof(sats_management_t));
  hypothesis_t *hyp = (hypothesis_t *)memory_pool_add(amb_test->pool);
  memcpy(hyp, &((hypothesis_t) {.N = {1,2,4,5}}), sizeof(hypothesis_t));
  hyp = (hypothesis_t *)memory_pool_add(amb_test->pool);
  memcpy(hyp, &((hypothesis_t) {.N = {1,2,4,6}}), sizeof(hypothesis_t));
}

START_TEST(test_update_sats_reset)
{
  /* Tests that end with a reset */
  ambiguity_test_t amb_test;
  /* Make sure the ambiguity test fails when given only one sdiff */
  amb_test_setup(&amb_test);
  fail_unless(0 == ambiguity_update_sats(&amb_test, 1,
                         NULL, NULL, NULL, NULL, NULL, false));
  fail_unless(memory_pool_n_allocated(amb_test.pool) == 1);
  fail_unless(amb_test.sats.num_sats == 0);
  /* Make sure the ambiguity test fails when given zero sdiffs */
  amb_test_setup(&amb_test);
  fail_unless(0 == ambiguity_update_sats(&amb_test, 0,
                         NULL, NULL, NULL, NULL, NULL, false));
  fail_unless(memory_pool_n_allocated(amb_test.pool) == 1);
  fail_unless(amb_test.sats.num_sats == 0);
}
END_TEST

static bool sats_man_match(sats_management_t *sats1, sats_management_t *sats2)
{
  if (sats1->num_sats != sats2->num_sats) return false;
  for (u8 i=0; i < sats1->num_sats; i++) {
    if (sats1->prns[i] != sats2->prns[i]) return false;
  }
  return true;
}

static bool hyp_pool_match(ambiguity_test_t *amb_test1, ambiguity_test_t *amb_test2)
{
  if (memory_pool_n_allocated(amb_test1->pool) != memory_pool_n_allocated(amb_test2->pool))
      return false;
  /* TODO use some kind of "and $ zipWith (==) pool1 pool2" */
  return true;
}

static bool amb_tests_match(ambiguity_test_t *amb_test1, ambiguity_test_t *amb_test2)
{
  if (!sats_man_match(&amb_test1->sats, &amb_test2->sats)) return false;
  return hyp_pool_match(amb_test1, amb_test2);
}

START_TEST(test_update_sats_intermediate_reset)
{
  /* Tests that reset somewhere along the way.
     (Any of these should give the same output as an empty pool */
  ambiguity_test_t amb_test, amb_test2;
  double u[5*5];
  matrix_eye(5, u);
  double d[5] = {1e-5,1e-5,1e-5,1e-5,1e-5};
  double float_mean[5] = {0,-1,-2,-3,-4};
  // double float_mean[5] = {0,0,0,0,0};
  sdiff_t sdiffs[6];
  sdiffs[0].prn = 0; sdiffs[0].snr = 1;
  sdiffs[1].prn = 1; sdiffs[1].snr = 0;
  sdiffs[2].prn = 2; sdiffs[2].snr = 0;
  sdiffs[3].prn = 3; sdiffs[3].snr = 0;
  sdiffs[4].prn = 4; sdiffs[4].snr = 0;
  sdiffs[5].prn = 5; sdiffs[5].snr = 0;
  sats_management_t float_sats = {.num_sats = 6, .prns={3, 0,1,2,4,5}};
  /* num_sats == 0 should trigger a reset */
  create_empty_ambiguity_test(&amb_test);
  memcpy(&amb_test.sats, &((sats_management_t) {.num_sats = 0}),
    sizeof(sats_management_t));
  memory_pool_add(amb_test.pool);
  create_ambiguity_test(&amb_test2);
  ambiguity_update_sats(&amb_test, 6, sdiffs,
    &float_sats, float_mean, u, d, false);
  ambiguity_update_sats(&amb_test2, 6, sdiffs,
    &float_sats, float_mean, u, d, false);
  fail_unless(amb_tests_match(&amb_test, &amb_test2));
  /* num_sats == 1 should trigger a reset */
  create_empty_ambiguity_test(&amb_test);
  memcpy(&amb_test.sats, &((sats_management_t) {.num_sats = 1}),
    sizeof(sats_management_t));
  memory_pool_add(amb_test.pool);
  create_ambiguity_test(&amb_test2);
  ambiguity_update_sats(&amb_test, 6, sdiffs,
    &float_sats, float_mean, u, d, false);
  ambiguity_update_sats(&amb_test2, 6, sdiffs,
    &float_sats, float_mean, u, d, false);
  fail_unless(amb_tests_match(&amb_test, &amb_test2));
  /* Should also fail if the KF and the amb test diverge */
  amb_test_setup(&amb_test);
  create_ambiguity_test(&amb_test2);
  ambiguity_update_sats(&amb_test, 6, sdiffs,
    &float_sats, float_mean, u, d, false);
  ambiguity_update_sats(&amb_test2, 6, sdiffs,
    &float_sats, float_mean, u, d, false);
  fail_unless(amb_tests_match(&amb_test, &amb_test2));
  /* No DDs in intersection should trigger a reset */
  sdiffs[0].prn = 0;
  sdiffs[1].prn = 6;
  sdiffs[2].prn = 7;
  sdiffs[3].prn = 8;
  sdiffs[4].prn = 9;
  sdiffs[5].prn = 10;
  amb_test_setup(&amb_test);
  create_ambiguity_test(&amb_test2);
  ambiguity_update_sats(&amb_test, 6, sdiffs,
    &float_sats, float_mean, u, d, false);
  ambiguity_update_sats(&amb_test2, 6, sdiffs,
    &float_sats, float_mean, u, d, false);
  fail_unless(amb_tests_match(&amb_test, &amb_test2));
  for (double x = -100; x < 100; x++) {
    fail_unless(within_epsilon(log(1 + exp(x)), log1p(exp(x))));
  }
}
END_TEST

static s32 cmp_hyps_lex(void *arg, element_t *a, element_t *b)
{
  u8 num_sats = *((u8 *) arg);
  hypothesis_t *hyp1 = (hypothesis_t *)a;
  hypothesis_t *hyp2 = (hypothesis_t *)b;

  for (u8 i=0; i < MAX(1,num_sats)-1; i++) {
    if (hyp1->N[i] > hyp2->N[i]) return 1;
    if (hyp1->N[i] < hyp2->N[i]) return -1;
  }
  return 0;
}

static void set_prns(u8 n, sdiff_t *sdiffs, s32 *prns)
{
  for (u8 i=0; i < n; i++) {
    sdiffs[i].prn = prns[i];
  }
}

static void set_snrs(u8 n, sdiff_t *sdiffs, double *snrs)
{
  for (u8 i=0; i < n; i++) {
    sdiffs[i].snr = snrs[i];
  }
}

static void test_ambiguity_update_sats(ambiguity_test_t *amb_test,
  sats_management_t float_sats, double *d, double *u, double *float_mean,
  sdiff_t *sdiffs, bool is_bad_meas,
  sats_management_t expected_sats, s32 expected_num_hyps, s32 *expected_hyps)
{
  amb_test_setup(amb_test);
  ambiguity_update_sats(amb_test, float_sats.num_sats, sdiffs,
    &float_sats, float_mean, u, d, is_bad_meas);
  memory_pool_sort(amb_test->pool, &amb_test->sats.num_sats, &cmp_hyps_lex);
  hypothesis_t hyps[expected_num_hyps];
  memory_pool_to_array(amb_test->pool, hyps);
  fail_unless(memory_pool_n_allocated(amb_test->pool) == expected_num_hyps);
  fail_unless(amb_test->sats.num_sats == expected_sats.num_sats);
  fail_unless(!memcmp(amb_test->sats.prns, expected_sats.prns, expected_sats.num_sats * sizeof(u8)));
  u8 num_dds = CLAMP_DIFF(expected_sats.num_sats, 1);
  for (u8 i=0; i < expected_num_hyps; i++) {
    fail_unless(!memcmp(hyps[i].N, &expected_hyps[i*num_dds], num_dds*sizeof(s32)));
  }
}

START_TEST(test_ambiguity_update_sats_cases)
{
  ambiguity_test_t amb_test;

  double u[5*5];
  matrix_eye(4, u);
  double d[5] = {1e-5,1e-5,1e-5,1e-5,1e-5};
  double float_mean[6];
  sdiff_t sdiffs[6];

  /* Drop a sat and combine two hypotheses */
  set_prns(4, sdiffs, (s32[4]){1,2,3,4});
  test_ambiguity_update_sats(&amb_test,
    (sats_management_t) {.num_sats = 4, .prns={3, 1,2,4}}, d, u, float_mean,
    sdiffs, false, (sats_management_t) {.num_sats = 4, .prns = {3, 1,2,4}},
    1, (s32[3]){1,2,4});

  /* Drop a sat but keep both hypotheses */
  set_prns(4, sdiffs, (s32[4]) {1,2,3,5});
  test_ambiguity_update_sats(&amb_test,
    (sats_management_t) {.num_sats = 4, .prns={3, 1,2,5}}, d, u, float_mean,
    sdiffs, false, (sats_management_t) {.num_sats = 4, .prns = {3, 1,2,5}},
    2, (s32[2*3]){1,2,5,
                  1,2,6});

  /* Drop a sat and lose the reference */
  set_prns(4, sdiffs,    (s32[4]) {1,2,4,5});
  set_snrs(4, sdiffs, (double[4]) {0,0,1,0});
  test_ambiguity_update_sats(&amb_test,
    (sats_management_t) {.num_sats = 4, .prns={5, 1,2,4}}, d, u, float_mean,
    sdiffs, false, (sats_management_t) {.num_sats = 4, .prns = {4, 1,2,5}},
    2, (s32[2*3]){-3,-2,1,
                  -3,-2,2});

  /* Do nothing */
  set_prns(5, sdiffs, (s32[5]){1,2,3,4,5});
  test_ambiguity_update_sats(&amb_test,
    (sats_management_t) {.num_sats = 5, .prns={3, 1,2,4,5}}, d, u, float_mean,
    sdiffs, false, (sats_management_t) {.num_sats = 5, .prns={3, 1,2,4,5}},
    2, (s32[2*4]){1,2,4,5,
                  1,2,4,6});

  /* Add a sat */
  matrix_eye(5, u);
  memcpy(d, (double[5]) {1e-5,1e-5,1e-5,2,1e-5}, 5*sizeof(double));
  memcpy(float_mean, (double [5]) {1,2,3,4,5}, 5*sizeof(double));
  set_prns(6, sdiffs, (s32[6]){1,2,3,4,5,6});
  test_ambiguity_update_sats(&amb_test,
    (sats_management_t) {.num_sats = 6, .prns={3, 1,2,4,5,6}}, d, u, float_mean,
    sdiffs, false, (sats_management_t) {.num_sats = 6, .prns={3, 1,2,4,5,6}},
    6, (s32[6*5]){1,2,4,5,4,
                  1,2,4,5,5,
                  1,2,4,5,6,
                  1,2,4,6,4,
                  1,2,4,6,5,
                  1,2,4,6,6});

  /* Don't add a sat that we otherwise would, because it was a bad measurement. */
  test_ambiguity_update_sats(&amb_test,
    (sats_management_t) {.num_sats = 6, .prns={3, 1,2,4,5,6}}, d, u, float_mean,
    sdiffs, true, (sats_management_t) {.num_sats = 5, .prns={3, 1,2,4,5}},
    2, (s32[2*4]){1,2,4,5,
                  1,2,4,6});

  /* Start over */
}
END_TEST

// START_TEST(test_update_sats_inclusion_projection_identity)
// {
//   THIS DOES NOT HOLD ANYMORE :(, because the KF and amb test can diverge.
// }

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

  ambiguity_update_sats(&amb_test, num_sdiffs, sdiffs, &float_sats, NULL, NULL, NULL, false);
  fail_unless(amb_test.sats.num_sats == 3);
  fail_unless(amb_test.sats.prns[0] == 4);
  fail_unless(amb_test.sats.prns[1] == 1);
  fail_unless(amb_test.sats.prns[2] == 2);
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

  amb_test.sats.num_sats = 4;

  sdiff_t sdiffs[4] = {{.prn = 1, .snr = 0},
                       {.prn = 2, .snr = 0},
                       {.prn = 4, .snr = 1}};
  u8 num_sdiffs = 3;

  for (u32 i=0; i<3; i++) {
    hypothesis_t *hyp = (hypothesis_t *)memory_pool_add(amb_test.pool);
    fail_unless(hyp != 0, "Null pointer returned by memory_pool_add");
    for (u8 j=0; j<amb_test.sats.num_sats-1; j++) {
      hyp->N[j] = sizerand(5);
    }
    hyp->ll = frand(0, 1);
  }

  sdiff_t sdiffs_with_ref_first[4];
  ambiguity_update_reference(&amb_test, num_sdiffs, sdiffs, sdiffs_with_ref_first);

  /* TODO: Check results of rebase. */
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
  fail_unless(pool_size == 1);

  /* Include. This one should succeed and add 5 sats. */
  flag = ambiguity_sat_inclusion(&amb_test, 0, &float_sats, mean, u, d);
  pool_size = memory_pool_n_allocated(amb_test.pool);
  fail_unless(flag == 1);
  fail_unless(pool_size == 625);

  /* Include again. This one should succeed and add 1 more sat. */
  flag = ambiguity_sat_inclusion(&amb_test, 0, &float_sats, mean, u, d);
  pool_size = memory_pool_n_allocated(amb_test.pool);
  fail_unless(flag == 1);
  fail_unless(pool_size == 945);

  /* Include again. This one should fail. */
  flag = ambiguity_sat_inclusion(&amb_test, 0, &float_sats, mean, u, d);
  pool_size = memory_pool_n_allocated(amb_test.pool);
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
  tcase_add_test(tc_core, test_bad_measurements);
  // TODO add back
  //tcase_add_test(tc_core, test_update_sats_rebase);
  (void) test_update_sats_rebase;
  tcase_add_test(tc_core, test_amb_sat_inclusion);
  tcase_add_test(tc_core, test_update_sats_reset);
  tcase_add_test(tc_core, test_update_sats_intermediate_reset);
  tcase_add_test(tc_core, test_ambiguity_update_sats_cases);
  suite_add_tcase(s, tc_core);

  return s;
}

