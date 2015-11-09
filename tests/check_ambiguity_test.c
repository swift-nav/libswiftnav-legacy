#include <check.h>
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
                                        {{.sat = 3},{.sat = 1},{.sat = 2},{.sat = 4}}}};
  sdiff_t sdiffs[4] = {{.sid = {.sat = 1}},
                       {.sid = {.sat = 2}},
                       {.sid = {.sat = 3}},
                       {.sid = {.sat = 4}}};
  u8 num_sdiffs = 4;

  ambiguity_update_sats(&amb_test, num_sdiffs, sdiffs, NULL, NULL, NULL, NULL, false);

  fail_unless(amb_test.sats.prns[0].sat == 3);
  fail_unless(amb_test.sats.prns[1].sat == 1);
  fail_unless(amb_test.sats.prns[2].sat == 2);
  fail_unless(amb_test.sats.prns[3].sat == 4);
}
END_TEST

/* Assure that when the measurement is bad, we drop (but don't add) sats */
START_TEST(test_bad_measurements)
{
  ambiguity_test_t amb_test;

  sdiff_t sdiffs[5] = {{.sid = {.sat = 1}},
                       {.sid = {.sat = 2}},
                       {.sid = {.sat = 3}},
                       {.sid = {.sat = 5}},
                       {.sid = {.sat = 6}}};
  u8 num_sdiffs = 5;


  sats_management_t float_sats = {.num_sats = 5,
                                  .prns = {{.sat = 3}, {.sat = 1}, {.sat = 2}, {.sat = 5}, {.sat = 6}}};
  double U[16];
  matrix_eye(4, U);
  double D[4] = {1, 1, 1, 1};
  double est[5] = {1, 2, 5, 6};

  sats_management_t amb_sats_init = {.num_sats = 5,
                                     .prns = {{.sat = 3}, {.sat = 1}, {.sat = 2}, {.sat = 4}, {.sat = 5}}};
  hypothesis_t hyp_init = {.N = {1,2,4,5}};

  create_empty_ambiguity_test(&amb_test);
  memcpy(&amb_test.sats, &amb_sats_init, sizeof(sats_management_t));
  hypothesis_t *hyp = (hypothesis_t *)memory_pool_add(amb_test.pool);
  memcpy(hyp, &hyp_init, sizeof(hypothesis_t));
  /* Test that with a good measurement, we get a projection and inclusion.
   * It should have dropped PRN 4 and include PRN 6. */
  ambiguity_update_sats(&amb_test, num_sdiffs, sdiffs, &float_sats, est, U, D, false);
  fail_unless(amb_test.sats.num_sats == 5);
  fail_unless(amb_test.sats.prns[0].sat == 3);
  fail_unless(amb_test.sats.prns[1].sat == 1);
  fail_unless(amb_test.sats.prns[2].sat == 2);
  fail_unless(amb_test.sats.prns[3].sat == 5);
  fail_unless(amb_test.sats.prns[4].sat == 6);
  /* Reset the amb_test to what it was before ambiguity_update_sats */
  create_empty_ambiguity_test(&amb_test);
  memcpy(&amb_test.sats, &amb_sats_init, sizeof(sats_management_t));
  hyp = (hypothesis_t *)memory_pool_add(amb_test.pool);
  memcpy(hyp, &hyp_init, sizeof(hypothesis_t));
  /* Test that with a bad measurement, we get (only) a projection.
   * It should have dropped PRN 4 and NOT include PRN 6. */
  ambiguity_update_sats(&amb_test, num_sdiffs, sdiffs, &float_sats, est, U, D, true);
  fail_unless(amb_test.sats.num_sats == 4);
  fail_unless(amb_test.sats.prns[0].sat == 3);
  fail_unless(amb_test.sats.prns[1].sat == 1);
  fail_unless(amb_test.sats.prns[2].sat == 2);
  fail_unless(amb_test.sats.prns[3].sat == 5);
  /* And it should have dropped PRN's value from the hypothesis,
   * which should still be there and still be the only one. */
  hyp = (hypothesis_t *) amb_test.pool->allocated_nodes_head->elem;
  fail_unless(ambiguity_test_n_hypotheses(&amb_test) == 1);
  fail_unless(hyp->N[0] == 1);
  fail_unless(hyp->N[1] == 2);
  fail_unless(hyp->N[2] == 5);
}
END_TEST

/* Assure that when we've lost the reference, we choose a new one and rebase everything. */
START_TEST(test_update_sats_rebase)
{
  ambiguity_test_t amb_test;
  create_empty_ambiguity_test(&amb_test);

  amb_test.sats.num_sats = 4;
  amb_test.sats.prns[0].sat = 3;
  amb_test.sats.prns[1].sat = 1;
  amb_test.sats.prns[2].sat = 2;
  amb_test.sats.prns[3].sat = 4;

  sdiff_t sdiffs[3] = {{.sid = {.sat = 1}, .snr = 0},
                       {.sid = {.sat = 2}, .snr = 0},
                       {.sid = {.sat = 4}, .snr = 1}};
  u8 num_sdiffs = 3;

  hypothesis_t *hyp = (hypothesis_t *)memory_pool_add(amb_test.pool);
  hyp->N[0] = 0;
  hyp->N[1] = 1;
  hyp->N[2] = 2;

  sats_management_t float_sats = {.num_sats = 3};

  ambiguity_update_sats(&amb_test, num_sdiffs, sdiffs, &float_sats, NULL, NULL, NULL, false);
  fail_unless(amb_test.sats.num_sats == 3);
  fail_unless(amb_test.sats.prns[0].sat == 4);
  fail_unless(amb_test.sats.prns[1].sat == 1);
  fail_unless(amb_test.sats.prns[2].sat == 2);
  fail_unless(hyp->N[0] == -2);
  fail_unless(hyp->N[1] == -1);
}
END_TEST

START_TEST(test_ambiguity_update_reference)
{
  srandom(1);

  ambiguity_test_t amb_test = {.sats = {.num_sats = 4,
                                        .prns = {{.sat = 3},{.sat = 1},{.sat = 2},{.sat = 4}}}};
  create_empty_ambiguity_test(&amb_test);

  amb_test.sats.num_sats = 4;

  sdiff_t sdiffs[4] = {{.sid = {.sat = 1}, .snr = 0},
                       {.sid = {.sat = 2}, .snr = 0},
                       {.sid = {.sat = 4}, .snr = 1}};
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
                                        .prns = {{.sat = 3},{.sat = 1},{.sat = 2}}}};
  sdiff_t sdiffs[4] = {{.sid = {.sat = 1}},
                       {.sid = {.sat = 2}},
                       {.sid = {.sat = 3}},
                       {.sid = {.sat = 4}}};
  u8 num_sdiffs = 4;
  fail_unless(!sats_match(&amb_test, num_sdiffs, sdiffs));

  num_sdiffs = 3;
  fail_unless(sats_match(&amb_test, num_sdiffs, sdiffs));

  sdiffs[0].sid.sat = 22;
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
  signal_t prns[dim];
  memset(prns, 0, sizeof(prns));
  for (u8 i = 0; i < dim+1; i++) {
    prns[i].sat = i;
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
  memcpy(float_sats.prns, prns, (dim+1) * sizeof(signal_t));

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
  suite_add_tcase(s, tc_core);

  return s;
}

