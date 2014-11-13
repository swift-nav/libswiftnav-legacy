/*
 * Copyright (C) 2014 Swift Navigation Inc.
 * Contact: Ian Horn <ian@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#include <clapack.h>
#include <inttypes.h>
#include <cblas.h>
#include <stdio.h>
#include <string.h>
#include "ambiguity_test.h"
#include "common.h"
#include "constants.h"
#include "linear_algebra.h"
#include "single_diff.h"
#include "amb_kf.h"
#include "lambda.h"

#define RAW_PHASE_BIAS_VAR 0
#define DECORRELATED_PHASE_BIAS_VAR 0
#define NUM_SEARCH_STDS 5
#define LOG_PROB_RAT_THRESHOLD -90

#define DEBUG_AMBIGUITY_TEST 0


/** \defgroup ambiguity_test Integer Ambiguity Resolution
 * Integer ambiguity resolution using bayesian hypothesis testing.
 * \{ */

void create_ambiguity_test(ambiguity_test_t *amb_test)
{
  static u8 pool_buff[MAX_HYPOTHESES*(sizeof(hypothesis_t) + sizeof(void *))];
  static memory_pool_t pool;
  amb_test->pool = &pool;
  memory_pool_init(amb_test->pool, MAX_HYPOTHESES, sizeof(hypothesis_t), pool_buff);
  amb_test->sats.num_sats = 0;
  amb_test->amb_check.initialized = 0;
}

/** A memory pool filter for clearing out all elements of a pool.
 */
s8 filter_all(void *arg, element_t *elem)
{
  (void) arg; (void) elem;
  return false;
}


void reset_ambiguity_test(ambiguity_test_t *amb_test) //TODO is this even necessary? we may only need create_ambiguity_test
{
  if (DEBUG_AMBIGUITY_TEST) {
      printf("<RESET_AMBIGUITY_TEST>\n");
  }
  u8 x = 0;
  memory_pool_filter(amb_test->pool, (void *) &x, &filter_all);
  /* Initialize pool with single element with num_dds = 0, i.e.
   * zero length N vector, i.e. no satellites. When we take the
   * product of this single element with the set of new satellites
   * we will just get a set of elements corresponding to the new sats. */
  hypothesis_t *empty_element = (hypothesis_t *)memory_pool_add(amb_test->pool);
  /* Start with ll = 0, just for the sake of argument. */
  empty_element->ll = 0;
  amb_test->sats.num_sats = 0;
  amb_test->amb_check.initialized = 0;
  if (DEBUG_AMBIGUITY_TEST) {
      printf("</RESET_AMBIGUITY_TEST>\n");
  }
}


void destroy_ambiguity_test(ambiguity_test_t *amb_test)
{
  memory_pool_destroy(amb_test->pool);
}


/** Prints a matrix of doubles.
 *
 * Used primarily for printing matrices in a gdb session.
 * Assumes row-major (C-style) storage.
 *
 * \param m     The matrix to be printed.
 * \param _r    The number of rows to be printed.
 * \param _c    The number of columns in the matrix.
 */
void print_double_mtx(double *m, u32 _r, u32 _c) {
    for (u32 _i = 0; _i < (_r); _i++) {
      printf(" [% 12lf", (m)[_i*(_c) + 0]);
      for (u32 _j = 1; _j < (_c); _j++)
        printf(" % 12lf", (m)[_i*(_c) + _j]);
      printf("]\n");
    }
}


/** Prints the matrix of Pearson correlation coefficients of a covariance matrix.
 *
 * Takes in a double valued covariance matrix and prints the pearson correlation
 * coefficients of each pair of variables:
 *
 * \f[
 *    \rho_{ij} = \frac{cov_{ij}}{\sigma_i \sigma_j}
 * \f]
 *
 * \param m     The covariance matrix to be transformed.
 * \param dim   The dimension of the covariance matrix.
 */
void print_pearson_mtx(double *m, u32 dim) {
    for (u32 _i = 0; _i < dim; _i++) {
      printf(" [% 12lf", m[_i*dim + 0] / sqrt(m[_i*dim + _i]) / sqrt(m[0]));
      for (u32 _j = 1; _j < dim; _j++)
        printf(" % 12lf", m[_i*dim + _j] / sqrt(m[_i*dim + _i]) / sqrt(m[_j * dim + _j]));
      printf("]\n");
    }
}


/** Prints the difference between two s32 valued matrices.
 *
 * Given two s32 valued matrices of the same shape, prints their difference
 * (first matrix minus second).
 *
 * \param m     The number of rows in the matrices.
 * \param n     The number of columns in the matrices.
 * \param mat1  The matrix to be subtracted from.
 * \param mat2  The matrix to be subtracted.
 */
void print_s32_mtx_diff(u32 m, u32 n, s32 *mat1, s32 *mat2)
{
  for (u32 i=0; i < m; i++) {
    for (u32 j=0; j < n; j++) {
      printf("%"PRId32", ", mat1[i*n + j] - mat2[i*n + j]);
    }
    printf("\n");
  }
  printf("\n");
}

/** Prints a s32 valued matrix.
 *
 * \param m     The number of rows to be printed in the matrices.
 * \param n     The number of columns in the matrix.
 * \param mat1  The matrix to be printed.
 */
void print_s32_mtx(u32 m, u32 n, s32 *mat)
{
  for (u32 i=0; i < m; i++) {
    for (u32 j=0; j < n; j++) {
      printf("%"PRId32", ", mat[i*n + j]);
    }
    printf("\n");
  }
  printf("\n");
}

/** Prints the result of a matrix-vector product of s32's.
 * Given a matrix M and vector v of s32's, prints the result of M*v.
 *
 * \param m   The number of rows in the matrix M.
 * \param n   The number of columns in the matrix M and elements of v.
 * \param M   The matrix M to be multiplied.
 * \param v   The vector v to be multiplied.
 */
void print_s32_gemv(u32 m, u32 n, s32 *M, s32 *v)
{
  s32 mv[m];
  memset(mv, 0, m * sizeof(s32));
  printf("[");
  for (u32 i=0; i < m; i++) {
    for (u32 j=0; j < n; j++) {
      mv[i] += M[i*n + j] * v[j];
    }
    if (i+1 == m) {
      printf("%"PRId32" ]\n", mv[i]);
    }
    else {
      printf("%"PRId32", ", mv[i]);
    }
  }
}

/** Gets the hypothesis out of an ambiguity test struct, if there is only one.
 *
 * Given an ambiguity_test_t, if that test has only one hypothesis allocated,
 * copies it to the hyp_N and returns 0. Otherwise, returns -1.
 *
 * \param amb_test    A pointer to the ambiguity_test_t struct to query against.
 * \param hyp_N       A pointer to the array of s32's to represent the output hyp.
 * \return            0 iff there only one hypothesis in the pool, -1 otherwise.
 */
s8 get_single_hypothesis(ambiguity_test_t *amb_test, s32 *hyp_N)
{
  if (memory_pool_n_allocated(amb_test->pool) == 1) {
    hypothesis_t hyp;
    memory_pool_to_array(amb_test->pool, &hyp);
    memcpy(hyp_N, hyp.N, (amb_test->sats.num_sats-1) * sizeof(s32));
    return 0;
  }
  return -1;
}

/** A struct to be used as the x in a memory pool fold, checking set membership.
 * Used in fold_contains().
 */
typedef struct {
  u8 num_dds;             /**< The number of double differences in the hypotheses. */
  s32 N[MAX_CHANNELS-1];  /**< The ambiguity vector being searched for. */
  u8 found;               /**< Whether or not the ambiguity vector is found yet. */
} fold_contains_t;

/** A memory pool fold method to check if a hypothesis is in the pool.
 *
 * Should be used in memory_pool_fold().
 *
 * \param x    The accumulator. Contains the hypothesis and whether it's found.
 * \param elem The hypothesis being folded against (checked for match).
 */
void fold_contains(void *x, element_t *elem)
{
  fold_contains_t *acc = (fold_contains_t *) x;
  hypothesis_t *hyp = (hypothesis_t *) elem;
  if (acc->found == 1) {
    return;
  }
  for (u8 i=0; i<acc->num_dds; i++) {
    if (hyp->N[i] != acc->N[i]) {
      return;
    }
  }
  acc->found = 1;
}

/** Tests whether an ambiguity test has a particular hypothesis.
 *
 * \param amb_test    The test to check against.
 * \param ambs        The ambiguity hypothesis to look for.
 * \return            Whether or not the pool contains ambs.
 */
u8 ambiguity_test_pool_contains(ambiguity_test_t *amb_test, double *ambs)
{
  fold_contains_t acc;
  acc.num_dds = amb_test->sats.num_sats-1;
  for (u8 i=0; i<acc.num_dds; i++) {
    acc.N[i] = lround(ambs[i]);
  }
  acc.found = 0;
  memory_pool_fold(amb_test->pool, (void *) &acc, &fold_contains);
  return acc.found;
}


/** A struct to be used in a memory pool fold to find the most likely hypothesis.
 * Should be used in fold_mle().
 */
typedef struct {
  u8 started;            /**< Whether the test has actually been started yet. */
  double max_ll;         /**< The likelihood of the most likely hypothesis so far. */
  u8 num_dds;            /**< The number of double differences in the hypotheses being tested. */
  s32 N[MAX_CHANNELS-1]; /**< The most likely hypothesis so far. */
} fold_mle_t; // fold omelette

/** A memory pool fold method to find the max likelihood estimate of the ambiguities.
 *
 * Should be used in memory_pool_fold().
 *
 * \param x     The fold accumulator. Contains the MLE and supporting fields.
 * \param elem  The hypothesis to fold.
 */
void fold_mle(void *x, element_t *elem)
{
  hypothesis_t *hyp = (hypothesis_t *) elem;
  fold_mle_t *mle = (fold_mle_t *) x;
  if (mle->started == 0 || hyp->ll > mle->max_ll) {
    mle->started = 1;
    mle->max_ll = hyp->ll;
    memcpy(mle->N, hyp->N, mle->num_dds * sizeof(s32));
  }
}

/** Performs max likelihood estimation on an ambiguity test.
 *
 * Assuming an ambiguity test already has hypotheses, finds the MLE hypothesis.
 * WARNING: Does not handle the case where the test is not populated.
 *
 * \param amb_test  The ambiguity test to perform MLE in.
 * \param ambs      The output MLE hypothesis.
 */
void ambiguity_test_MLE_ambs(ambiguity_test_t *amb_test, s32 *ambs)
{
  fold_mle_t mle;
  mle.started = 0;
  mle.num_dds = MAX(1,amb_test->sats.num_sats)-1;
  memory_pool_fold(amb_test->pool, (void *) &mle, &fold_mle);
  memcpy(ambs, mle.N, mle.num_dds * sizeof(s32));
}

/** Starts a hypothesis test for integer ambiguity resolution (IAR).
 *
 * Initializes an IAR test.
 *
 * Adds an empty element the ambiguity test pool, then adds
 * hypotheses for various satellites if they fit.
 *
 * Assumes that the IAR hypothesis pool is empty.
 * Assumes that the dimension of the float state is the number of DDs.
 *
 * \param amb_test    The amb_test struct to be updated
 * \param state_dim   The length of the float state/number of DDs.
 * \param float_prns  The list of prns associated with the float DD estimate.
 * \param sdiffs      Not currently used.
 * \param float_mean  The float DD ambiguity estimate.
 * \param float_cov   The float DD ambiguity estimate's covariance.
 * \param DE_mtx      The differenced ECEF unit vectors pointing to the satellite.
 * \param obs_cov     The observation covariance matrix of the phase.
 */
void init_ambiguity_test(ambiguity_test_t *amb_test, u8 state_dim, u8 *float_prns, sdiff_t *sdiffs, double *float_mean,
                         double *float_cov, double *DE_mtx, double *obs_cov)
{
  (void) sdiffs;
  u8 num_dds = state_dim;
  double float_cov_N[num_dds * num_dds];

  /* Re-shape float_cov to contain just ambiguity states. */
  for (u8 i=0; i<num_dds; i++) {
    for (u8 j=0; j<num_dds; j++) {
      float_cov_N[i*num_dds + j] = float_cov[i*state_dim + j]; //TODO this is just a memcpy
    }
  }

  /* Initialize pool with single element with num_dds = 0, i.e.
   * zero length N vector, i.e. no satellites. When we take the
   * product of this single element with the set of new satellites
   * we will just get a set of elements corresponding to the new sats. */
  hypothesis_t *empty_element = (hypothesis_t *)memory_pool_add(amb_test->pool); // only in init
  /* Start with ll = 0, just for the sake of argument. */
  empty_element->ll = 0; // only in init
  amb_test->sats.num_sats = 0; // only in init
  s32 Z_inv[num_dds * num_dds];
  s32 lower_bounds[num_dds];
  s32 upper_bounds[num_dds];
  u8 num_dds_to_add;
  s8 add_any_sats =  determine_sats_addition(amb_test,
                                             float_cov_N, num_dds, &float_mean[6],
                                             lower_bounds, upper_bounds, &num_dds_to_add,
                                             Z_inv);

  if (add_any_sats == 1) {
    add_sats(amb_test, float_prns[0], num_dds_to_add, &float_prns[1], lower_bounds, upper_bounds, Z_inv);
    /* Update the rest of the amb_test state with the new sats. */
    init_residual_matrices(&amb_test->res_mtxs, num_dds, DE_mtx, obs_cov); // only in init
  }

}


/** Updates the IAR process with new measurements.
 *
 * Updates the satellites being tested, adding and removing hypotheses as needed.
 *
 * Assumes that the float estimates include all of the sats in the sdiffs.
 * Assumes that the dimension of the float estimates is one less than the number of sats.
 * Assumes that the sdiffs are arranged entirely in increasing order by PRN.
 *
 * \param ref_ecef    The ecef coordinate to pretend we are at to use relative to the sats.
 * \param phase_var   The variance of the carrier phase measurements.
 * \param code_var    The variance of the code pseudorange measurements.
 * \param amb_test    The ambiguity test to update.
 * \param state_dim   The dimension of the float state.
 * \param float_sats  The satellites being represented in the float state.
 * \param sdiffs      The single differenced measurements/sat positions of all sats tracked.
 * \param float_mean  The float estimate of the integer ambiguities.
 * \param float_cov_U The U in the UDU' decomposition of the covariance of the float estimate.
 * \param float_cov_D The D in the UDU' decomposition of the covariance of the float estimate.
 *
 *  INVALIDATES unanimous ambiguities
 */
void update_ambiguity_test(double ref_ecef[3], double phase_var, double code_var,
                           ambiguity_test_t *amb_test, u8 state_dim, sdiff_t *sdiffs,
                           u8 changed_sats)
{
  if (DEBUG_AMBIGUITY_TEST) {
    printf("<UPDATE_AMBIGUITY_TEST>\n");
  }
  u8 num_sdiffs = state_dim + 1;

  if (amb_test->sats.num_sats < 5) {
    if (DEBUG_AMBIGUITY_TEST) {
      printf("</UPDATE_AMBIGUITY_TEST>\n");
    }
    return;
  }

  sdiff_t ambiguity_sdiffs[amb_test->sats.num_sats];
  double ambiguity_dd_measurements[2*(amb_test->sats.num_sats-1)];
  s8 valid_sdiffs = make_ambiguity_dd_measurements_and_sdiffs(
      amb_test, num_sdiffs, sdiffs, ambiguity_dd_measurements, ambiguity_sdiffs);

  /* Error */
  if (valid_sdiffs != 0) {
    printf("update_ambiguity_test: Invalid sdiffs. return code: %i\n", valid_sdiffs);
    for (u8 k=0; k < num_sdiffs; k++) {
      printf("%u, ", sdiffs[k].prn);
    }
    printf("}\n");
    print_sats_management_short(&amb_test->sats);
    return;
  }

  if (1 == 1 || changed_sats == 1) { //TODO add logic about when to update DE
    double DE_mtx[(amb_test->sats.num_sats-1) * 3];
    assign_de_mtx(amb_test->sats.num_sats, ambiguity_sdiffs, ref_ecef, DE_mtx);
    double obs_cov[(amb_test->sats.num_sats-1) * (amb_test->sats.num_sats-1) * 4];
    memset(obs_cov, 0, (amb_test->sats.num_sats-1) * (amb_test->sats.num_sats-1) * 4 * sizeof(double));
    u8 num_dds = amb_test->sats.num_sats-1;
    for (u8 i=0; i<num_dds; i++) {
      for (u8 j=0; j<num_dds; j++) {
        u8 i_ = i+num_dds;
        u8 j_ = j+num_dds;
        if (i==j) {
          obs_cov[i*2*num_dds + j] = phase_var * 2;
          obs_cov[i_*2*num_dds + j_] = code_var * 2;
        }
        else {
          obs_cov[i*2*num_dds + j] = phase_var;
          obs_cov[i_*2*num_dds + j_] = code_var;
        }
      }
    }
    // MAT_PRINTF(DE_mtx, ((u32) amb_test->sats.num_sats-1), 3);
    // MAT_PRINTF(obs_cov, 2*num_dds, 2*num_dds);
    init_residual_matrices(&amb_test->res_mtxs, amb_test->sats.num_sats-1, DE_mtx, obs_cov);
  }

  test_ambiguities(amb_test, ambiguity_dd_measurements);

  if (DEBUG_AMBIGUITY_TEST) {
    printf("</UPDATE_AMBIGUITY_TEST>\n");
  }
}


/** Returns the number of hypotheses currently in the ambiguity test.
 * \param amb_test    The ambiguity test whose number of hypotheses we want.
 * \return            The number of hypotheses currently in the ambiguity test.
 */
u32 ambiguity_test_n_hypotheses(ambiguity_test_t *amb_test)
{
  return memory_pool_n_allocated(amb_test->pool);
}

/** A memory pool filter to do most everything related to the estimation algorithm.
 * To be used by update_and_get_max_ll() and filter_and_renormalize().
 */
typedef struct {
  u8 num_dds;                                 /**< Number of ambiguities. */
  double r_vec[2*MAX_CHANNELS-5];             /**< Transformed measurement to check hypotheses against. */
  double max_ll;                              /**< The greatest log likelihood in the pool so far. */
  residual_mtxs_t *res_mtxs;                  /**< Matrices necessary for testing hypotheses. */
  unanimous_amb_check_t *unanimous_amb_check; /**< A struct to check which int ambs are agreed upon among all hyps. */
} hyp_filter_t;

/** A mapAccum styled function to update the hypothesis log-likelihoods and find the greatest LL.
 * Simultaneously performs a map, doing a Bayesian update of the log likelihoods
 * of each hypothesis, while performing a fold on those updated log likelihoods
 * to find the likelihood of the MLE hypothesis.
 *
 * To be given to memory_pool_fold().
 *
 * \param x     Points to a hyp_filter_t containing the accumulator and everything needed for the map.
 * \param elem  The hypothesis to be mapAccum'd.
 */
void update_and_get_max_ll(void *x_, element_t *elem) {
  hyp_filter_t *x = (hyp_filter_t *) x_;
  hypothesis_t *hyp = (hypothesis_t *) elem;
  double hypothesis_N[x->num_dds];

  for (u8 i=0; i < x->num_dds; i++) {
    hypothesis_N[i] = hyp->N[i];
  }

  hyp->ll += get_quadratic_term(x->res_mtxs, x->num_dds, hypothesis_N, x->r_vec);
  x->max_ll = MAX(x->max_ll, hyp->ll);
}

/** Keeps track of which integer ambiguities are uninimously agreed upon in the pool.
 * \param num_dds   The number of DDs in each hypothesis. (Used to initialize amb_check).
 * \param hyp       The hypothesis to be checked against.
 * \param amb_check Keeps track of which ambs are still unanimous and their values.
 */
void check_unanimous_ambs(u8 num_dds, hypothesis_t *hyp,
                          unanimous_amb_check_t *amb_check)
{
  if (amb_check->initialized) {
    u8 j = 0; // index in newly constructed amb_check matches
    for (u8 i = 0; i < amb_check->num_matching_ndxs; i++) {
      if (amb_check->ambs[i] == hyp->N[amb_check->matching_ndxs[i]]) {
        if (i != j) { //  j <= i necessarily
          amb_check->matching_ndxs[j] = amb_check->matching_ndxs[i];
          amb_check->ambs[j] = amb_check->ambs[i];
        }
        j++;
      }
    }
    amb_check->num_matching_ndxs = j;
  } else {
    amb_check->initialized = 1;
    amb_check->num_matching_ndxs = num_dds;
    for (u8 i=0; i < num_dds; i++) {
      amb_check->matching_ndxs[i] = i;
    }
    memcpy(amb_check->ambs, hyp->N, num_dds * sizeof(s32));
  }
}

/** A combination fold/map/filter to renormalize the log likelihoods and remove unlikely hyps.
 *
 * The log likelihood of the hypotheses are filtered against a threshold.
 * Those hypotheses that make the cut are (1) normalized such that the MLE has
 * value 0, making them logs of the probability ratio against the MLE hyp,
 * and (2) checked to see if any of their integer ambiguities are unanimous
 * amongst all hypotheses in the test pool.
 *
 * The thresholding is done before the normalization for both numerical
 * stability, and so that hypotheses which are just REALLY BAD are removed,
 * even if they are the best we have. This is a kinda arbitrary choice of how
 * to do things. Maybe we should see if it has practical implications?
 *
 * This function is used by memory_pool_filter().
 *
 * \param arg     The accumulator/map/filter information
 * \param elem    The hypothesis to be mapped/filtered/folded against.
 * \return        Whether or not this hypothesis made the cut.
 */
s8 filter_and_renormalize(void *arg, element_t *elem) {
  hypothesis_t *hyp = (hypothesis_t *) elem;

  u8 keep_it = (hyp->ll > LOG_PROB_RAT_THRESHOLD);
  if (keep_it) {
    hyp->ll -= ((hyp_filter_t *) arg)->max_ll;
  }
  return keep_it;
}

static void _check_unanimous(void *arg, element_t *elem)
{
  hypothesis_t *hyp = (hypothesis_t *) elem;

  check_unanimous_ambs(((hyp_filter_t *) arg)->num_dds, hyp,
                       ((hyp_filter_t *) arg)->unanimous_amb_check);
}

void update_unanimous_ambiguities(ambiguity_test_t *amb_test)
{
  hyp_filter_t x;
  if (amb_test->sats.num_sats <= 1) {
    amb_test->amb_check.num_matching_ndxs = 0;
    return;
  }
  x.num_dds = amb_test->sats.num_sats-1;
  x.unanimous_amb_check = &amb_test->amb_check;
  x.unanimous_amb_check->initialized = 0;

  memory_pool_map(amb_test->pool, (void *) &x, &_check_unanimous);
}

/* Updates the IAR hypothesis pool log likelihood ratios and filters them.
 *  It assumes that the observations are structured to match the amb_test sats.
 *  INVALIDATES unanimous ambiguities
 */
void test_ambiguities(ambiguity_test_t *amb_test, double *dd_measurements)
{
  if (DEBUG_AMBIGUITY_TEST) {
    printf("<TEST_AMBIGUITIES>\n");
  }
  hyp_filter_t x;
  x.num_dds = amb_test->sats.num_sats-1;
  assign_r_vec(&amb_test->res_mtxs, x.num_dds, dd_measurements, x.r_vec);
  // VEC_PRINTF(x.r_vec, amb_test->res_mtxs.res_dim);
  x.max_ll = -1e20; //TODO get the first element, or use this as threshold to restart test
  x.res_mtxs = &amb_test->res_mtxs;
  x.unanimous_amb_check = &amb_test->amb_check;
  x.unanimous_amb_check->initialized = 0;

  memory_pool_fold(amb_test->pool, (void *) &x, &update_and_get_max_ll);
  /*memory_pool_map(amb_test->pool, &x.num_dds, &print_hyp);*/
  memory_pool_filter(amb_test->pool, (void *) &x, &filter_and_renormalize);
  if (memory_pool_empty(amb_test->pool)) {
      /* Initialize pool with single element with num_dds = 0, i.e.
     * zero length N vector, i.e. no satellites. When we take the
     * product of this single element with the set of new satellites
     * we will just get a set of elements corresponding to the new sats. */
    hypothesis_t *empty_element = (hypothesis_t *)memory_pool_add(amb_test->pool);
    /* Start with ll = 0, just for the sake of argument. */
    empty_element->ll = 0;
    amb_test->sats.num_sats = 0;
    amb_test->amb_check.initialized = 0;
  }
  if (DEBUG_AMBIGUITY_TEST) {
    memory_pool_map(amb_test->pool, &x.num_dds, &print_hyp);
    printf("num_unanimous_ndxs=%u\n</TEST_AMBIGUITIES>\n", x.unanimous_amb_check->num_matching_ndxs);
  }
}

/* This says whether we can use the ambiguity_test to resolve a position in 3-space.
 */
u8 ambiguity_iar_can_solve(ambiguity_test_t *amb_test)
{
  return amb_test->amb_check.initialized &&
         amb_test->amb_check.num_matching_ndxs >= 3;
}

bool is_prn_set(u8 len, u8 *prns)
{
  if (len == 0) {
    return true;
  }
  u8 current = prns[0];
  for (u8 i = 1; i < len; i++) {
    if (prns[i] <= current) {
      return false;
    }
    current = prns[i];
  }
  return true;
}

/* Constructs the double differenced measurements and sdiffs needed for IAR.
 * This requires that all IAR prns are in the sdiffs used (sdiffs is a superset
 * of IAR's PRNs), that the sdiffs are ordered by prn, and that the IAR
 * non_ref_prns are also ordered (both ascending).
 *
 * The ambiguity_dd_measurements output will be ordered like the non_ref_prns
 * with the whole vector of carrier phase elements coming before the
 * pseudoranges. The amb_sdiffs will have the ref sat first, then the rest in
 * ascending order like the non_ref_prns.
 *
 * \param ref_prn                    The current reference PRN for the IAR.
 * \param non_ref_prns               The rest of the current PRNs for the IAR.
 * \param num_dds                    The number of dds used in the IAR
 *                                   (length of non_ref_prns).
 * \param num_sdiffs                 The number of sdiffs being passed in.
 * \param sdiffs                     The sdiffs to pull measurements out of.
 * \param ambiguity_dd_measurements  The output vector of DD measurements
 *                                   to be used to update the IAR.
 * \param amb_sdiffs                 The sdiffs that correspond to the IAR PRNs.
 * \return 0 if the input sdiffs are superset of the IAR sats,
 *        -1 if they are not,
 *        -2 if non_ref_prns is not an ordered set.
 */
s8 make_dd_measurements_and_sdiffs(u8 ref_prn, u8 *non_ref_prns, u8 num_dds,
                                   u8 num_sdiffs, sdiff_t *sdiffs,
                                   double *ambiguity_dd_measurements, sdiff_t *amb_sdiffs)
{
  if (DEBUG_AMBIGUITY_TEST) {
    printf("<MAKE_DD_MEASUREMENTS_AND_SDIFFS>\n");
    printf("ref_prn = %u\nnon_ref_prns = {", ref_prn);
    for (u8 i=0; i < num_dds; i++) {
      printf("%d, ", non_ref_prns[i]);
    }
    printf("}\nnum_dds = %u\nnum_sdiffs = %u\nsdiffs[*].prn = {", num_dds, num_sdiffs);
    for (u8 i=0; i < num_sdiffs; i++) {
      printf("%d, ", sdiffs[i].prn);
    }
    printf("}\n");
  }

  if (!is_prn_set(num_dds, non_ref_prns)) {
    printf("There is disorder in the amb_test sats.\n");
    printf("amb_test sat prns = {%u, ", ref_prn);
    for (u8 k=0; k < num_dds; k++) {
      printf("%u, ", non_ref_prns[k]);
    }
    printf("}\n");
    return -2;
  }

  double ref_phase = 0;
  double ref_pseudorange = 0;
  u8 i=0;
  u8 j=0;
  u8 found_ref = 0;
  /* Go through the sdiffs, pulling out the measurements of the non-ref amb sats
   * and the reference sat. */
  while (i < num_dds) {
    if (non_ref_prns[i] == sdiffs[j].prn) {
      /* When we find a non-ref sat, we fill in the next measurement. */
      memcpy(&amb_sdiffs[i+1], &sdiffs[j], sizeof(sdiff_t));
      ambiguity_dd_measurements[i] = sdiffs[j].carrier_phase;
      ambiguity_dd_measurements[i+num_dds] = sdiffs[j].pseudorange;
      i++;
      j++;
    } else if (ref_prn == sdiffs[j].prn) {
      /* when we find the ref sat, we copy it over and raise the FOUND flag */
      memcpy(&amb_sdiffs[0], &sdiffs[j], sizeof(sdiff_t));
      ref_phase =  sdiffs[j].carrier_phase;
      ref_pseudorange = sdiffs[j].pseudorange;
      j++;
      found_ref = 1;
    }
    else if (non_ref_prns[i] > sdiffs[j].prn) {
      /* If both sets are ordered, and we increase j (and possibly i), and the
       * i prn is higher than the j one, it means that the i one might be in the 
       * j set for higher j, and that the current j prn isn't in the i set. */
      j++;
    } else {
      /* if both sets are ordered, and we increase j (and possibly i), and the
       * j prn is higher than the i one, it means that the j one might be in the 
       * i set for higher i, and that the current i prn isn't in the j set.
       * This means a sat in the IAR's sdiffs isn't in the sdiffs.
       * */
      return -1;
    }
  }
  /* This awkward case deals with the situation when sdiffs and sats have the
   * same satellites only the ref of amb_test.sats is the last PRN in sdiffs.
   * This case is never checked for j = num_dds as i only runs to num_dds-1. */
  /* TODO: This function could be refactored to be a lot clearer. */
  while (!found_ref && j < num_sdiffs ) {
    if (ref_prn == sdiffs[j].prn) {
      memcpy(&amb_sdiffs[0], &sdiffs[j], sizeof(sdiff_t));
      ref_phase =  sdiffs[j].carrier_phase;
      ref_pseudorange = sdiffs[j].pseudorange;
      found_ref = 1;
    }
    j++;
  }

  if (found_ref == 0) {
    return -1;
  }
  for (u8 i=0; i < num_dds; i++) {
    ambiguity_dd_measurements[i] -= ref_phase;
    ambiguity_dd_measurements[i+num_dds] -= ref_pseudorange;
  }
  if (DEBUG_AMBIGUITY_TEST) {
    printf("amb_sdiff_prns = {");
    for (i = 0; i < num_dds; i++) {
      printf("%u, ", amb_sdiffs[i].prn);
    }
    printf("}\ndd_measurements = {");
    for (i=0; i < 2 * num_dds; i++) {
      printf("%f, \t", ambiguity_dd_measurements[i]);
    }
    printf("}\n</MAKE_DD_MEASUREMENTS_AND_SDIFFS>\n");
  }
  return 0;
}


/** Make the DD measurements and sdiffs that correspond to the resolved DDs in amb_test.
 * Assuming that amb_test has not been modified since the last check to see
 * whether it can resolve a fixed position (that is, amb_test->amb_check is up
 * to date, this will make a set of sdiffs and DD measurements that correspond
 * to the resolved sats.
 *
 * Returns 0 if the input sdiff sats are a superset of the resolved IAR sats.
 * Returns -1 otherwise.
 *
 * \param amb_test                  The local amb_test struct. Must have a
 *                                  current amb_check.
 * \param num_sdiffs                The number of sdiffs passed in.
 * \param sdiffs                    The sdiffs from which we make amb_sdiffs.
 * \param ambiguity_dd_measurements DD measurement vector for the amb_test's
 *                                  unanimously resolved sats.
 * \param amb_sdiffs                sdiffs that match the amb_test's unanimously
 *                                  resolved sats constructed from the input
 *                                  sdiffs.
 * \return error code; see make_dd_measurements_and_sdiffs docstring.
 */
s8 make_ambiguity_resolved_dd_measurements_and_sdiffs(
            ambiguity_test_t *amb_test, u8 num_sdiffs, sdiff_t *sdiffs,
            double *ambiguity_dd_measurements, sdiff_t *amb_sdiffs)
{
  if (DEBUG_AMBIGUITY_TEST) {
    printf("<MAKE_AMBIGUITY_RESOLVED_DD_MEASUREMENTS_AND_SDIFFS>\n");
    printf("amb_test->sats.prns = {");
    for (u8 i=0; i < amb_test->sats.num_sats; i++) {
      printf("%u, ", amb_test->sats.prns[i]);
    }
    printf("}\n");
  }

  u8 ref_prn = amb_test->sats.prns[0];
  u8 num_dds = amb_test->amb_check.num_matching_ndxs;
  u8 non_ref_prns[num_dds];
  for (u8 i=0; i < num_dds; i++) {
    non_ref_prns[i] = amb_test->sats.prns[1 + amb_test->amb_check.matching_ndxs[i]];
    if (DEBUG_AMBIGUITY_TEST) {
      printf("non_ref_prns[%u] = %u, \t (ndx=%u) \t amb[%u] = %"PRId32"\n",
             i, non_ref_prns[i], amb_test->amb_check.matching_ndxs[i],
             i, amb_test->amb_check.ambs[i]);
    }
  }
  s8 valid_sdiffs =
    make_dd_measurements_and_sdiffs(ref_prn, non_ref_prns, num_dds, num_sdiffs,
                                    sdiffs, ambiguity_dd_measurements,
                                    amb_sdiffs);
  if (DEBUG_AMBIGUITY_TEST) {
    printf("</MAKE_AMBIGUITY_RESOLVED_DD_MEASUREMENTS_AND_SDIFFS>\n");
  }
  return valid_sdiffs;
}


/** Make the DD measurements and sdiffs that correspond to the DDs in amb_test.
 * This will make a set of sdiffs and DD measurements that correspond to the
 * amb_test's sats.
 *
 * \param amb_test                  The local amb_test struct. Must have the
 *                                  current amb_check.
 * \param num_sdiffs                The number of sdiffs passed in.
 * \param sdiffs                    The sdiffs from which we make amb_sdiffs.
 * \param ambiguity_dd_measurements DD measurement vector for the amb_test's
 *                                  sats.
 * \param amb_sdiffs                sdiffs that match the amb_test's sats
 *                                  constructed from the input sdiffs.
 * \return error code. see make_dd_measurements_and_sdiffs docstring.
 */
s8 make_ambiguity_dd_measurements_and_sdiffs(ambiguity_test_t *amb_test, u8 num_sdiffs, sdiff_t *sdiffs,
                                               double *ambiguity_dd_measurements, sdiff_t *amb_sdiffs)
{
  if (DEBUG_AMBIGUITY_TEST) {
    printf("<MAKE_AMBIGUITY_DD_MEASUREMENTS_AND_SDIFFS>\n");
  }
  u8 ref_prn = amb_test->sats.prns[0];
  u8 *non_ref_prns = &amb_test->sats.prns[1];
  u8 num_dds = MAX(1, amb_test->sats.num_sats)-1;
  s8 valid_sdiffs = make_dd_measurements_and_sdiffs(ref_prn, non_ref_prns,
                                  num_dds, num_sdiffs, sdiffs,
                                  ambiguity_dd_measurements, amb_sdiffs);
  if (DEBUG_AMBIGUITY_TEST) {
    printf("</MAKE_AMBIGUITY_DD_MEASUREMENTS_AND_SDIFFS>\n");
  }
  return valid_sdiffs;
}


s8 sats_match(ambiguity_test_t *amb_test, u8 num_sdiffs, sdiff_t *sdiffs)
{
  if (DEBUG_AMBIGUITY_TEST) {
    printf("<SATS_MATCH>\namb_test.sats.prns = {");
    for (u8 i=0; i< amb_test->sats.num_sats; i++) {
      printf("%u, ", amb_test->sats.prns[i]);
    }
    printf("}\nsdiffs[*].prn      = {");
    for (u8 i=0; i < num_sdiffs; i++) {
      printf("%u, ", sdiffs[i].prn);
    }
    printf("}\n");
  }
  if (amb_test->sats.num_sats != num_sdiffs) {
    if (DEBUG_AMBIGUITY_TEST) {
      printf("sats don't match (different length)\n</SATS_MATCH>\n");
    }
    return 0;
  }
  u8 *prns = amb_test->sats.prns;
  u8 amb_ref = amb_test->sats.prns[0];
  u8 j=0;
  for (u8 i = 1; i<amb_test->sats.num_sats; i++) { //TODO will not having a j condition cause le fault du seg?
    if (j >= num_sdiffs) {
      if (DEBUG_AMBIGUITY_TEST) {
        printf("sats don't match\n</SATS_MATCH>\n");
      }
      return 0;
    }
    if (prns[i] == sdiffs[j].prn) {
      j++;
    }
    else if (amb_ref == sdiffs[j].prn) {
      j++;
      i--;
    }
    else {
      if (DEBUG_AMBIGUITY_TEST) {
        printf("sats don't match\n</SATS_MATCH>\n");
      }
      return 0;
    }
  }
  if (DEBUG_AMBIGUITY_TEST) {
    printf("sats match\n</SATS_MATCH>\n");
  }
  return 1;
}

typedef struct {
  u8 num_sats;
  u8 old_prns[MAX_CHANNELS];
  u8 new_prns[MAX_CHANNELS];
} rebase_prns_t;

void rebase_hypothesis(void *arg, element_t *elem) //TODO make it so it doesn't have to do all these lookups every time
{
  rebase_prns_t *prns = (rebase_prns_t *) arg;
  u8 num_sats = prns->num_sats;
  u8 *old_prns = prns->old_prns;
  u8 *new_prns = prns->new_prns;

  hypothesis_t *hypothesis = (hypothesis_t *)elem;

  u8 old_ref = old_prns[0];
  u8 new_ref = new_prns[0];

  s32 new_N[num_sats-1];
  s32 index_of_new_ref_in_old = find_index_of_element_in_u8s(num_sats, new_ref, &old_prns[1]);
  s32 val_for_new_ref_in_old_basis = hypothesis->N[index_of_new_ref_in_old];
  for (u8 i=0; i<num_sats-1; i++) {
    u8 new_prn = new_prns[1+i];
    if (new_prn == old_ref) {
      new_N[i] = - val_for_new_ref_in_old_basis;
    }
    else {
      s32 index_of_this_sat_in_old_basis = find_index_of_element_in_u8s(num_sats, new_prn, &old_prns[1]);
      new_N[i] = hypothesis->N[index_of_this_sat_in_old_basis] - val_for_new_ref_in_old_basis;
    }
  }
  memcpy(hypothesis->N, new_N, (num_sats-1) * sizeof(s32));
}

u8 ambiguity_update_reference(ambiguity_test_t *amb_test, u8 num_sdiffs, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  if (DEBUG_AMBIGUITY_TEST) {
    printf("<AMBIGUITY_UPDATE_REFERENCE>\n");
  }
  u8 changed_ref = 0;
  u8 old_prns[amb_test->sats.num_sats];
  memcpy(old_prns, amb_test->sats.prns, amb_test->sats.num_sats * sizeof(u8));

  // print_sats_management_short(&amb_test->sats);
  s8 sats_management_code = rebase_sats_management(&amb_test->sats, num_sdiffs, sdiffs, sdiffs_with_ref_first);
  // print_sats_management_short(&amb_test->sats);
  if (sats_management_code != OLD_REF) {
    if (DEBUG_AMBIGUITY_TEST) {
      printf("updating iar reference sat\n");
    }
    changed_ref = 1;
    u8 new_prns[amb_test->sats.num_sats];
    memcpy(new_prns, amb_test->sats.prns, amb_test->sats.num_sats * sizeof(u8));

    rebase_prns_t prns = {.num_sats = amb_test->sats.num_sats};
    memcpy(prns.old_prns, old_prns, amb_test->sats.num_sats * sizeof(u8));
    memcpy(prns.new_prns, new_prns, amb_test->sats.num_sats * sizeof(u8));
    memory_pool_map(amb_test->pool, &prns, &rebase_hypothesis);
  }
  if (DEBUG_AMBIGUITY_TEST) {
    printf("</AMBIGUITY_UPDATE_REFERENCE>\n");
  }
  return changed_ref;
}

typedef struct {
  u8 num_ndxs;
  u8 intersection_ndxs[MAX_CHANNELS-1];
} intersection_ndxs_t;

s32 projection_comparator(void *arg, element_t *a, element_t *b)
{
  intersection_ndxs_t *intersection_struct = (intersection_ndxs_t *) arg;
  u8 num_ndxs = intersection_struct->num_ndxs;
  u8 *intersection_ndxs = intersection_struct->intersection_ndxs;

  hypothesis_t *hyp_a = (hypothesis_t *) a;
  hypothesis_t *hyp_b = (hypothesis_t *) b;

  for (u8 i=0; i<num_ndxs; i++) {
    u8 ndxi = intersection_ndxs[i];
    s32 ai = hyp_a->N[ndxi];
    s32 bi = hyp_b->N[ndxi];
    if (ai == bi) {
      continue;
    }
    if (ai < bi) {
      return -1;
    }
    if (ai > bi) {
      return 1;
    }
  }
  return 0;
}

void projection_aggregator(element_t *new_, void *x_, u32 n, element_t *elem_)
{
  intersection_ndxs_t *x = (intersection_ndxs_t *)x_;
  hypothesis_t *new = (hypothesis_t *)new_;
  hypothesis_t *elem = (hypothesis_t *)elem_;

  if (n == 0) {
    for (u8 i=0; i<x->num_ndxs; i++) {
      u8 ndxi = x->intersection_ndxs[i];
      new->N[i] = elem->N[ndxi];
    }
    new->ll = elem->ll;
  }
  else {
    new->ll += log(1 + exp(elem->ll - new->ll));
    // new->ll = MAX(new->ll, elem->ll);
  }

}

u8 ambiguity_sat_projection(ambiguity_test_t *amb_test, u8 num_dds_in_intersection, u8 *dd_intersection_ndxs)
{
  if (DEBUG_AMBIGUITY_TEST) {
    printf("<AMBIGUITY_SAT_PROJECTION>\n");
  }
  u8 num_dds_before_proj = MAX(1, amb_test->sats.num_sats) - 1;
  if (num_dds_before_proj == num_dds_in_intersection) {
    if (DEBUG_AMBIGUITY_TEST) {
      printf("no need for projection\n</AMBIGUITY_SAT_PROJECTION>\n");
    }
    return 0;
  }

  intersection_ndxs_t intersection = {.num_ndxs = num_dds_in_intersection};
  memcpy(intersection.intersection_ndxs, dd_intersection_ndxs, num_dds_in_intersection * sizeof(u8));


  printf("IAR: %"PRIu32" hypotheses before projection\n", memory_pool_n_allocated(amb_test->pool));
  /*memory_pool_map(amb_test->pool, &num_dds_before_proj, &print_hyp);*/
  memory_pool_group_by(amb_test->pool,
                       &intersection, &projection_comparator,
                       &intersection, sizeof(intersection),
                       &projection_aggregator);
  printf("IAR: updates to %"PRIu32"\n", memory_pool_n_allocated(amb_test->pool));
  /*memory_pool_map(amb_test->pool, &num_dds_in_intersection, &print_hyp);*/
  u8 work_prns[MAX_CHANNELS];
  memcpy(work_prns, amb_test->sats.prns, amb_test->sats.num_sats * sizeof(u8));
  for (u8 i=0; i<num_dds_in_intersection; i++) {
    amb_test->sats.prns[i+1] = work_prns[dd_intersection_ndxs[i]+1];
  }
  amb_test->sats.num_sats = num_dds_in_intersection+1;
  if (DEBUG_AMBIGUITY_TEST) {
    printf("</AMBIGUITY_SAT_PROJECTION>\n");
  }
  return 1;
}


u8 ambiguity_sat_inclusion(ambiguity_test_t *amb_test, u8 num_dds_in_intersection,
                             sats_management_t *float_sats, double *float_mean, double *float_cov_U, double *float_cov_D)
{
  if (DEBUG_AMBIGUITY_TEST) {
    printf("<AMBIGUITY_SAT_INCLUSION>\n");
  }
  if (float_sats->num_sats <= num_dds_in_intersection + 1 || float_sats->num_sats < 2) {
    if (DEBUG_AMBIGUITY_TEST) {
      printf("no need for inclusion\n</AMBIGUITY_SAT_INCLUSION>\n");
    }
    return 0;
  }

  u32 state_dim = float_sats->num_sats-1;
  double float_cov[state_dim * state_dim];
  matrix_reconstruct_udu(state_dim, float_cov_U, float_cov_D, float_cov);
  u8 float_prns[float_sats->num_sats];
  memcpy(float_prns, float_sats->prns, float_sats->num_sats * sizeof(u8));
  double N_mean[float_sats->num_sats-1];
  memcpy(N_mean, float_mean, (float_sats->num_sats-1) * sizeof(double));
  if (amb_test->sats.num_sats >= 2 && amb_test->sats.prns[0] != float_sats->prns[0]) {
    u8 old_prns[float_sats->num_sats];
    memcpy(old_prns, float_sats->prns, float_sats->num_sats * sizeof(u8));
    // memcpy(N_mean, &float_mean[6], (float_sats->num_sats-1) * sizeof(double));
    set_reference_sat_of_prns(amb_test->sats.prns[0], float_sats->num_sats, float_prns);
    rebase_mean_N(N_mean, float_sats->num_sats, old_prns, float_prns);
    rebase_covariance_sigma(float_cov, float_sats->num_sats, old_prns, float_prns);
  }
  double N_cov[(float_sats->num_sats-1) * (float_sats->num_sats-1)];
  memcpy(N_cov, float_cov, state_dim * state_dim * sizeof(double)); //TODO we can just use N_cov throughout
  //by now float_prns has the correct reference, as do N_cov and N_mean

  // MAT_PRINTF(float_cov, state_dim, state_dim);
  // MAT_PRINTF(N_cov, (u8)(float_sats->num_sats-1), (u8)(float_sats->num_sats-1));

  // printf("pearson mtx of float cov\n");
  // print_pearson_mtx(float_cov, state_dim);
  // printf("\n");

  // printf("pearson mtx of N cov\n");
  // print_pearson_mtx(N_cov, float_sats->num_sats-1);
  // printf("\n");

  //next we add the new sats
  //loop through the new sat sets in decreasing number (for now, in prn order), adding when it suits us

  /* first set up the largest prn lists of added sats
   * then get those prns covariances
   * then loop through:
   *    test if we'll add these
   *    if so:
   *      add them
   *      break loop
   *    else:
   *      decrease size of cov set
   */

  //first get all the prns we might add, and their covariances
  u8 i = 1;
  u8 j = 1;
  u8 num_addible_dds = 0;
  u8 ndxs_of_new_dds_in_float[MAX_CHANNELS-1];
  u8 new_dd_prns[MAX_CHANNELS-1];
  while (j < float_sats->num_sats) {
    if (i < amb_test->sats.num_sats && amb_test->sats.prns[i] == float_prns[j]) {
      i++;
      j++;
    } else { //else float_sats[j] is a new one
      ndxs_of_new_dds_in_float[num_addible_dds] = j-1;
      new_dd_prns[num_addible_dds] = float_prns[j];
      num_addible_dds++;
      j++;
    }
  }

  double addible_float_cov[num_addible_dds * num_addible_dds];
  double addible_float_mean[num_addible_dds];
  for (i=0; i < num_addible_dds; i++) {
    for (j=0; j < num_addible_dds; j++) {
      addible_float_cov[i*num_addible_dds + j] = N_cov[ndxs_of_new_dds_in_float[i]*(float_sats->num_sats-1) + ndxs_of_new_dds_in_float[j]];
    }
    addible_float_mean[i] = N_mean[ndxs_of_new_dds_in_float[i]];
  }
  // MAT_PRINTF(addible_float_cov, num_addible_dds, num_addible_dds);
  /*VEC_PRINTF(addible_float_mean, num_addible_dds);*/

  s32 Z_inv[num_addible_dds * num_addible_dds];
  s32 lower_bounds[num_addible_dds];
  s32 upper_bounds[num_addible_dds];
  u8 num_dds_to_add;
  s8 add_any_sats =  determine_sats_addition(amb_test,
                                             addible_float_cov, num_addible_dds, addible_float_mean,
                                             lower_bounds, upper_bounds, &num_dds_to_add,
                                             Z_inv);

  if (add_any_sats == 1) {
    add_sats(amb_test, float_prns[0], num_dds_to_add, new_dd_prns, lower_bounds, upper_bounds, Z_inv);
    if (DEBUG_AMBIGUITY_TEST) {
      printf("adding sats\n<AMBIGUITY_SAT_INCLUSION>\n");
    }
    return 1;
  } else {
    if (DEBUG_AMBIGUITY_TEST) {
      printf("not adding sats\n<AMBIGUITY_SAT_INCLUSION>\n");
    }
    return 0;
  }

}

u32 float_to_decor(ambiguity_test_t *amb_test,
                   double *addible_float_cov, u8 num_addible_dds,
                   double *addible_float_mean,
                   u8 num_dds_to_add,
                   s32 *lower_bounds, s32 *upper_bounds, double *Z)
{
  (void) amb_test;
  double added_float_cov[num_dds_to_add * num_dds_to_add];
  for (u8 i=0; i<num_dds_to_add; i++) {
    for (u8 j=0; j<num_dds_to_add; j++) {
      added_float_cov[i*num_dds_to_add + j] = addible_float_cov[i*num_addible_dds + j];
      #if RAW_PHASE_BIAS_VAR != 0
      if (i == j) {
        added_float_cov[i*num_dds_to_add + j] += RAW_PHASE_BIAS_VAR;
      }
      #endif
    }
  }

  // double Z[num_dds_to_add * num_dds_to_add];
  lambda_reduction(num_dds_to_add, added_float_cov, Z);

  double decor_float_cov_diag[num_dds_to_add];

  memset(decor_float_cov_diag, 0, num_dds_to_add * sizeof(double));

  for (u8 i=0; i < num_dds_to_add; i++) {
    for (u8 j=0; j < num_dds_to_add; j++) {
      for (u8 k=0; k < num_dds_to_add; k++) {
        decor_float_cov_diag[i] += Z[i*num_dds_to_add + j] * added_float_cov[j * num_dds_to_add + k] * Z[i*num_dds_to_add + k];
      }
    }
    #if DECORRELATED_PHASE_BIAS_VAR != 0
    decor_float_cov_diag[i] += DECORRELATED_PHASE_BIAS_VAR;
    #endif
  }

  double decor_float_mean[num_dds_to_add];
  memset(decor_float_mean, 0, num_dds_to_add * sizeof(double));
  for (u8 i=0; i < num_dds_to_add; i++) {
    for (u8 j=0; j < num_dds_to_add; j++) {
      decor_float_mean[i] += Z[i*num_dds_to_add + j] * addible_float_mean[j];
    }
    // printf("decor_float_mean[%u] = %f\n", i, decor_float_mean[i]);
  }

  u32 new_hyp_set_cardinality = 1;
  // s32 lower_bounds[num_dds_to_add];
  // s32 upper_bounds[num_dds_to_add];
  for (u8 i=0; i<num_dds_to_add; i++) {
    double search_distance = NUM_SEARCH_STDS * sqrt(decor_float_cov_diag[i]);
    // upper_bounds[i] = MAX(floor(float_mean[i] + search_distance), ceil(float_mean[i]));
    // lower_bounds[i] = MIN(ceil(float_mean[i] - search_distance), floor(float_mean[i]));
    upper_bounds[i] = lround(ceil(decor_float_mean[i] + search_distance));
    lower_bounds[i] = lround(floor(decor_float_mean[i] - search_distance));
    new_hyp_set_cardinality *= upper_bounds[i] - lower_bounds[i] + 1;
  }
  return new_hyp_set_cardinality;
}

s8 determine_sats_addition(ambiguity_test_t *amb_test,
                           double *float_N_cov, u8 num_float_dds, double *float_N_mean,
                           s32 *lower_bounds, s32 *upper_bounds, u8 *num_dds_to_add,
                           s32 *Z_inv)
{
  u8 num_current_dds = MAX(1, amb_test->sats.num_sats) - 1;
  u8 min_dds_to_add = MAX(1, 4 - num_current_dds); // num_current_dds + min_dds_to_add = 4 so that we have a nullspace projector

  u32 max_new_hyps_cardinality;
  s32 current_num_hyps = memory_pool_n_allocated(amb_test->pool);
  u32 max_num_hyps = memory_pool_n_elements(amb_test->pool);
  if (current_num_hyps <= 0) {
    max_new_hyps_cardinality = max_num_hyps;
  } else {
    max_new_hyps_cardinality = max_num_hyps / current_num_hyps;
  }

  // printf("pearson mtx of float N cov\n");
  // print_pearson_mtx(float_N_cov, num_float_dds);
  // printf("\n");

  *num_dds_to_add = num_float_dds;
  double Z[num_float_dds * num_float_dds];
  while (*num_dds_to_add >= min_dds_to_add) {
    u32 new_hyp_set_cardinality = float_to_decor(amb_test,
                                                 float_N_cov, num_float_dds,
                                                 float_N_mean,
                                                 *num_dds_to_add,
                                                 lower_bounds, upper_bounds, Z);
    if (new_hyp_set_cardinality <= max_new_hyps_cardinality) {
      double Z_inv_[*num_dds_to_add * *num_dds_to_add];
      matrix_inverse(*num_dds_to_add, Z, Z_inv_);
      /* TODO: Check return value of matrix_inverse to handle singular matrix. */
      for (u8 i=0; i < *num_dds_to_add; i++) {
        for (u8 j=0; j < *num_dds_to_add; j++) {
          Z_inv[i* *num_dds_to_add + j] = lround(Z_inv_[i* *num_dds_to_add + j]);
        }
      }
      return 1;
    }
    else {
      *num_dds_to_add -= 1;
    }
  }
  return -1;
}


// input/output: amb_test
// input:        num_sdiffs
// input:        sdiffs
// input:        float_sats
// input:        float_mean
// input:        float_cov_U
// input:        float_cov_D
// INVALIDATES unanimous ambiguities
u8 ambiguity_update_sats(ambiguity_test_t *amb_test, u8 num_sdiffs, sdiff_t *sdiffs,
                         sats_management_t *float_sats, double *float_mean, double *float_cov_U, double *float_cov_D)
{
  if (DEBUG_AMBIGUITY_TEST) {
    printf("<AMBIGUITY_UPDATE_SATS>\n");
  }
  if (num_sdiffs < 2) {
    create_ambiguity_test(amb_test);
    if (DEBUG_AMBIGUITY_TEST) {
      printf("< 2 sdiffs, starting over\n</AMBIGUITY_UPDATE_SATS>\n");
    }
    return 0; // I chose 0 because it doesn't lead to anything dynamic
  }
  //if the sats are the same, we're good
  u8 changed_sats = 0;
  if (!sats_match(amb_test, num_sdiffs, sdiffs)) {
    sdiff_t sdiffs_with_ref_first[num_sdiffs];
    if (amb_test->sats.num_sats >= 2) {
      if (ambiguity_update_reference(amb_test, num_sdiffs, sdiffs, sdiffs_with_ref_first)) {
       changed_sats=1;
      }
    } else {
      create_ambiguity_test(amb_test);//we don't have what we need
    }

    u8 intersection_ndxs[num_sdiffs];
    u8 num_dds_in_intersection = find_indices_of_intersection_sats(amb_test, num_sdiffs, sdiffs_with_ref_first, intersection_ndxs);

    if (amb_test->sats.num_sats > 1 && num_dds_in_intersection == 0) {
      create_ambiguity_test(amb_test); //TODO is create_ambiguity_test any better than reset_ambiguity_test here? does reset even need to exist
    }

    // u8 num_dds_in_intersection = ambiguity_order_sdiffs_with_intersection(amb_test, sdiffs, float_cov, intersection_ndxs);
    if (ambiguity_sat_projection(amb_test, num_dds_in_intersection, intersection_ndxs)) {
      changed_sats = 1;
    }
    if (ambiguity_sat_inclusion(amb_test, num_dds_in_intersection,
                                float_sats, float_mean, float_cov_U, float_cov_D)) {
      changed_sats = 1;
    }
  }
  if (DEBUG_AMBIGUITY_TEST) {
    printf("changed_sats = %u\n</AMBIGUITY_UPDATE_SATS>\n", changed_sats);
  }

  return changed_sats;
  /* TODO: Should we order by 'goodness'? Perhaps by volume of hyps? */
}

u8 find_indices_of_intersection_sats(ambiguity_test_t *amb_test, u8 num_sdiffs, sdiff_t *sdiffs_with_ref_first, u8 *intersection_ndxs)
{
  if (DEBUG_AMBIGUITY_TEST) {
    printf("<FIND_INDICES_OF_INTERSECTION_SATS>\n");
    printf("amb_test->sats.prns          = {");
    for (u8 i = 0; i < amb_test->sats.num_sats; i++) {
      printf("%u, ", amb_test->sats.prns[i]);
    }
    if (amb_test->sats.num_sats < 2) {
      printf("}\nsdiffs_with_ref_first not populated\n");
    }
    else {
      printf("}\nsdiffs_with_ref_first[*].prn = {");
      for (u8 i = 0; i < num_sdiffs; i++) {
        printf("%u, ", sdiffs_with_ref_first[i].prn);
      }
      printf("}\n");
    }
    printf("(i, j, k, amb_test->sats.prns[i], sdiffs_with_ref_first[j].prn)\n");
  }
  u8 i = 1;
  u8 j = 1;
  u8 k = 0;
  while (i < amb_test->sats.num_sats && j < num_sdiffs) {
    if (amb_test->sats.prns[i] == sdiffs_with_ref_first[j].prn) {
      if (DEBUG_AMBIGUITY_TEST) {
        printf("(%u, \t%u, \t%u, \t%u, \t%u)\t\t\tamb_test->sats.prns[i] == sdiffs_with_ref_first[j].prn; i,j,k++\n",
                i, j, k, amb_test->sats.prns[i], sdiffs_with_ref_first[j].prn);
      }
      intersection_ndxs[k] = i-1;
      i++;
      j++;
      k++;
    }
    else if (amb_test->sats.prns[i] < sdiffs_with_ref_first[j].prn) {
      if (DEBUG_AMBIGUITY_TEST) {
        printf("(%u, \t%u, \t%u, \t%u, \t%u)\t\t\tamb_test->sats.prns[i] <  sdiffs_with_ref_first[j].prn; i++\n",
                i, j, k, amb_test->sats.prns[i], sdiffs_with_ref_first[j].prn);
      }
      i++;
    }
    else {
      if (DEBUG_AMBIGUITY_TEST) {
        printf("(%u, \t%u, \t%u, \t%u, \t%u)\t\t\tamb_test->sats.prns[i] >  sdiffs_with_ref_first[j].prn; j++\n",
                i, j, k, amb_test->sats.prns[i], sdiffs_with_ref_first[j].prn);
      }
      j++;
    }
  }
  if (DEBUG_AMBIGUITY_TEST) {
    printf("</FIND_INDICES_OF_INTERSECTION_SATS>\n");
  }
  return k;
}

typedef struct {
  s32 upper_bounds[MAX_CHANNELS-1];
  s32 lower_bounds[MAX_CHANNELS-1];
  s32 counter[MAX_CHANNELS-1];
  u8 ndxs_of_old_in_new[MAX_CHANNELS-1];
  u8 ndxs_of_added_in_new[MAX_CHANNELS-1];
  u8 num_added_dds;
  u8 num_old_dds;
  s32 Z_inv[(MAX_CHANNELS-1) * (MAX_CHANNELS-1)];
} generate_hypothesis_state_t;

s8 generate_next_hypothesis(void *x_, u32 n)
{
  (void) n;
  generate_hypothesis_state_t *x = (generate_hypothesis_state_t *)x_;
  // printf("generate_next_hypothesis:\n");

  if (memcmp(x->upper_bounds, x->counter, x->num_added_dds * sizeof(s32)) == 0) {
    /* counter has reached upper_bound, terminate iteration. */
    return 0;
  }

  for (u8 i=0; i<x->num_added_dds; i++) {
    x->counter[i]++;
    if (x->counter[i] > x->upper_bounds[i]) {
      /* counter[i] has reached maximum, reset counter[i]
       * to lower[i] and 'carry' to next 'digit' */
      x->counter[i] = x->lower_bounds[i];
    } else {
      /* Incremented, so now we have the next counter value. */
      break;
    }
  }

  // printf("[");
  // for (u8 i=0; i<x->num_added_dds; i++) {
  //   printf("%d, ", x->counter[i]);
  // }
  // printf("]\n\n");

  return 1;
}

void hypothesis_prod(element_t *new_, void *x_, u32 n, element_t *elem_)
{
  (void) elem_;
  (void) n;
  generate_hypothesis_state_t *x = (generate_hypothesis_state_t *)x_;
  hypothesis_t *new = (hypothesis_t *) new_;

  u8 *ndxs_of_old_in_new = x->ndxs_of_old_in_new;
  u8 *ndxs_of_added_in_new = x->ndxs_of_added_in_new;

  s32 old_N[MAX_CHANNELS-1];
  memcpy(old_N, new->N, x->num_old_dds * sizeof(s32));

  // printf("counter: [");
  // for (u8 i=0; i < x->num_added_dds; i++) {
  //   printf("%d, ", x->counter[i]);
  // }
  // printf("]\n");
  // printf("old_N:ndx [");
  // for (u8 i=0; i < x->num_old_dds; i++) {
  //   printf("%d, ", ndxs_of_old_in_new[i]);
  // }
  // printf("]\n");
  // printf("old_N:    [");
  // for (u8 i=0; i < x->num_old_dds; i++) {
  //   printf("%d, ", old_N[i]);
  // }
  // printf("]\n");

  // printf("added_N:ndx [");
  // for (u8 i=0; i < x->num_added_dds; i++) {
  //   printf("%d, ", ndxs_of_added_in_new[i]);
  // }
  // printf("]\n");
  // printf("added_N:    [");
  // for (u8 i=0; i < x->num_added_dds; i++) {
  //   printf("%d, ", x->counter[i]);
  // }
  // printf("]\n");


  for (u8 i=0; i < x->num_old_dds; i++) {
    new->N[ndxs_of_old_in_new[i]] = old_N[i];
  }
  for (u8 i=0; i<x->num_added_dds; i++) {
    new->N[ndxs_of_added_in_new[i]] = 0;
    for (u8 j=0; j<x->num_added_dds; j++) {
      new->N[ndxs_of_added_in_new[i]] += x->Z_inv[i*x->num_added_dds + j] * x->counter[j];
      // printf("Z_inv[%u] = %d\n", i*x->num_added_dds + j, x->Z_inv[i*x->num_added_dds + j]);
      // printf("x->counter[[%u] = %d\n", j, x->counter[j]);
    }
  }
  // for (u8 i=0; i < x->num_added_dds; i++) {
  //   new->N[ndxs_of_added_in_new[i]] = x->counter[i];
  // }
  // printf("new_N:    [");
  // for (u8 i=0; i < x->num_added_dds + x->num_old_dds; i++) {
  //   printf("%d, ", new->N[i]);
  // }
  // printf("]\n\n");

  /* NOTE: new->ll remains the same as elem->ll as p := exp(ll) is invariant under a
   * constant multiplicative factor common to all hypotheses. TODO: reference^2 document (currently lives in page 3/5.6/2014 of ian's notebook) */
}

typedef struct {
  u8 num_added_dds;
  u8 num_old_dds;
  s32 Z_inv[(MAX_CHANNELS-1)*(MAX_CHANNELS-1)];
} recorrelation_params_t;

void recorrelate_added_sats(void *arg, element_t *elem_)
{
  hypothesis_t *elem = (hypothesis_t *) elem_;
  recorrelation_params_t *params = (recorrelation_params_t *)arg;

  /* elem->N <- [[I 0] [0 Z_inv]] . elem->N
   * where Z_inv is the inverse Lambda reduction matrix having
   * shape (num_added_dds, num_added_dds) and I has shape
   * (num_old_dds, num_old_dds) */

  s32 recorrelated_N[params->num_added_dds];
  memset(recorrelated_N, 0, params->num_added_dds * sizeof(s32));
  for (u8 i=0; i<params->num_added_dds; i++) {
    for (u8 j=0; j<params->num_added_dds; j++) {
      recorrelated_N[i] += params->Z_inv[i*params->num_added_dds + j] * elem->N[params->num_old_dds + j];
    }
  }
  memcpy(&elem->N[params->num_old_dds], recorrelated_N, params->num_added_dds * sizeof(s32));
}

void print_hyp(void *arg, element_t *elem)
{
  u8 num_dds = *( (u8 *) arg );

  hypothesis_t *hyp = (hypothesis_t *)elem;
  printf("[");
  for (u8 i=0; i< num_dds; i++) {
    printf("%"PRId32", ", hyp->N[i]);
  }
  printf("]: %f\n", hyp->ll);
}

void add_sats(ambiguity_test_t *amb_test,
              u8 ref_prn,
              u32 num_added_dds, u8 *added_prns,
              s32 *lower_bounds, s32 *upper_bounds,
              s32 *Z_inv)
{
  /* Make a generator that iterates over the new hypotheses. */
  generate_hypothesis_state_t x0;
  memcpy(x0.upper_bounds, upper_bounds, num_added_dds * sizeof(s32));
  memcpy(x0.lower_bounds, lower_bounds, num_added_dds * sizeof(s32));
  memcpy(x0.counter, lower_bounds, num_added_dds * sizeof(s32));
  // printf("upper = [");
  // for (u8 i=0; i<num_added_dds; i++) {
  //   printf("%d, ", x0.upper_bounds[i]);
  // }
  // printf("]\n");
  // printf("lower = [");
  // for (u8 i=0; i<num_added_dds; i++) {
  //   printf("%d, ", x0.lower_bounds[i]);
  // }
  // printf("]\n");

  x0.num_added_dds = num_added_dds;
  x0.num_old_dds = MAX(1,amb_test->sats.num_sats)-1;

  //then construct the mapping from the old prn indices into the new, and from the added prn indices into the new
  u8 i = 0;
  u8 j = 0;
  u8 k = 0;
  u8 old_prns[x0.num_old_dds];
  memcpy(old_prns, &amb_test->sats.prns[1], x0.num_old_dds * sizeof(u8));
  while (k < x0.num_old_dds + num_added_dds) { //TODO should this be one less, since its just DDs?
    if (j == x0.num_added_dds || (old_prns[i] < added_prns[j] && i != x0.num_old_dds)) {
      x0.ndxs_of_old_in_new[i] = k;
      amb_test->sats.prns[k+1] = old_prns[i];
      i++;
      k++;
    } else if (i == x0.num_old_dds || old_prns[i] > added_prns[j]) {
      x0.ndxs_of_added_in_new[j] = k;
      amb_test->sats.prns[k+1] = added_prns[j];
      j++;
      k++;
    } else {
      printf("This method is being used improperly. This shouldn't happen.\n");
      printf("old_prns = [");
      for (u8 ii=0; ii < x0.num_old_dds; ii++) {
        printf("%d, ",old_prns[ii]);
      }
      printf("]\n");
      printf("added_prns = [");
      for (u8 jj=0; jj < x0.num_old_dds; jj++) {
        printf("%d, ", added_prns[jj]);
      }
      printf("]\n");
      break;
    }
  }
  amb_test->sats.prns[0] = ref_prn;
  amb_test->sats.num_sats = k+1;

  if (x0.num_old_dds == 0 && memory_pool_n_allocated(amb_test->pool) == 0) {
    hypothesis_t *empty_element = (hypothesis_t *)memory_pool_add(amb_test->pool); // only in init
    /* Start with ll = 0, just for the sake of argument. */
    empty_element->ll = 0; // only in init
  }

  printf("IAR: %"PRIu32" hypotheses before inclusion\n", memory_pool_n_allocated(amb_test->pool));
  if (DEBUG_AMBIGUITY_TEST) {
    memory_pool_map(amb_test->pool, &x0.num_old_dds, &print_hyp);
  }
  memcpy(x0.Z_inv, Z_inv, num_added_dds * num_added_dds * sizeof(s32));
  /* Take the product of our current hypothesis state with the generator, recorrelating the new ones as we go. */
  memory_pool_product_generator(amb_test->pool, &x0, MAX_HYPOTHESES, sizeof(x0),
                                &generate_next_hypothesis, &hypothesis_prod);
  printf("IAR: updates to %"PRIu32"\n", memory_pool_n_allocated(amb_test->pool));
  if (DEBUG_AMBIGUITY_TEST) {
    memory_pool_map(amb_test->pool, &k, &print_hyp);
  }
}

void init_residual_matrices(residual_mtxs_t *res_mtxs, u8 num_dds, double *DE_mtx, double *obs_cov)
{
  res_mtxs->res_dim = num_dds + MAX(3, num_dds) - 3;
  res_mtxs->null_space_dim = MAX(3, num_dds) - 3;
  assign_phase_obs_null_basis(num_dds, DE_mtx, res_mtxs->null_projector);
  assign_residual_covariance_inverse(num_dds, obs_cov, res_mtxs->null_projector, res_mtxs->half_res_cov_inv);
  // MAT_PRINTF(res_mtxs->null_projector, res_mtxs->null_space_dim, num_dds);
  // MAT_PRINTF(res_mtxs->half_res_cov_inv, res_mtxs->res_dim, res_mtxs->res_dim);
}


// void QR_part1(integer m, integer n, double *A, double *tau)
// {
//   double w[1];
//   integer lwork = -1;
//   integer info;
//   integer jpvt[3];
//   memset(jpvt, 0, 3 * sizeof(integer));
//   dgeqp3_(&m, &n,
//           A, &m,
//           jpvt,
//           tau,
//           w, &lwork, &info);
//   lwork = round(w[0]);
//   double work[lwork];
//   dgeqp3_(&m, &n,
//           A, &m,
//           jpvt,
//           tau,
//           work, &lwork, &info); //set A = QR(A)
// }

// void QR_part2(integer m, integer n, double *A, double *tau)
// {
//   double w[1];
//   integer lwork = -1;
//   integer info;
//   dorgqr_(&m, &m, &n,
//           A, &m,
//           tau,
//           w, &lwork, &info);
//   lwork = round(w[0]);
//   double work[lwork];
//   dorgqr_(&m, &m, &n,
//           A, &m,
//           tau,
//           work, &lwork, &info);
// }

// void assign_phase_obs_null_basis(u8 num_dds, double *DE_mtx, double *q)
// {
//   // use either GEQRF or GEQP3. GEQP3 is the version with pivoting
//   // int dgeqrf_(__CLPK_integer *m, __CLPK_integer *n, __CLPK_doublereal *a, __CLPK_integer *
//   //       lda, __CLPK_doublereal *tau, __CLPK_doublereal *work, __CLPK_integer *lwork, __CLPK_integer *info)
//   // int dgeqp3_(__CLPK_integer *m, __CLPK_integer *n, __CLPK_doublereal *a, __CLPK_integer *
//   //       lda, __CLPK_integer *jpvt, __CLPK_doublereal *tau, __CLPK_doublereal *work, __CLPK_integer *lwork,
//   //        __CLPK_integer *info)

//   //DE is num_sats-1 by 3, need to transpose it to column major
//   double A[num_dds * num_dds];
//   for (u8 i=0; i<num_dds; i++) {
//     for (u8 j=0; j<3; j++) {
//       A[j*num_dds + i] = DE_mtx[i*3 + j]; //set A = Transpose(DE_mtx)
//     }
//   }
//   integer m = num_dds;
//   integer n = 3;
//   double tau[3];
//   QR_part1(m, n, A, tau);
//   QR_part2(m, n, A, tau);
//   memcpy(q, &A[3*num_dds], (num_dds-3) * num_dds * sizeof(double));
// }

void assign_residual_covariance_inverse(u8 num_dds, double *obs_cov, double *q, double *r_cov_inv) //TODO make this more efficient (e.g. via page 3/6.2-3/2014 of ian's notebook)
{
  integer dd_dim = 2*num_dds;
  integer res_dim = num_dds + MAX(3, num_dds) - 3;
  u32 nullspace_dim = MAX(3, num_dds) - 3;
  double q_tilde[res_dim * dd_dim];
  memset(q_tilde, 0, res_dim * dd_dim * sizeof(double));
  // MAT_PRINTF(obs_cov, dd_dim, dd_dim);

  // MAT_PRINTF(q, nullspace_dim, num_dds);
  for (u8 i=0; i<nullspace_dim; i++) {
    memcpy(&q_tilde[i*dd_dim], &q[i*num_dds], num_dds * sizeof(double));
  }
  // MAT_PRINTF(q_tilde, res_dim, dd_dim);
  for (u8 i=0; i<num_dds; i++) {
    q_tilde[(i+nullspace_dim)*dd_dim + i] = 1;
    q_tilde[(i+nullspace_dim)*dd_dim + i+num_dds] = -1 / GPS_L1_LAMBDA_NO_VAC;
  }
  // MAT_PRINTF(q_tilde, res_dim, dd_dim);

  //TODO make more efficient via the structure of q_tilde, and it's relation to the I + 1*1^T structure of the obs cov mtx
  double QC[res_dim * dd_dim];
  cblas_dsymm(CblasRowMajor, CblasRight, CblasUpper, //CBLAS_ORDER, CBLAS_SIDE, CBLAS_UPLO
              res_dim, dd_dim, // int M, int N
              1, obs_cov, dd_dim, // double alpha, double *A, int lda
              q_tilde, dd_dim, // double *B, int ldb
              0, QC, dd_dim); // double beta, double *C, int ldc
  // MAT_PRINTF(QC, res_dim, dd_dim);

  //TODO make more efficient via the structure of q_tilde, and it's relation to the I + 1*1^T structure of the obs cov mtx
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, // CBLAS_ORDER, CBLAS_TRANSPOSE transA, cBLAS_TRANSPOSE transB
              res_dim, res_dim, dd_dim, // int M, int N, int K
              2, QC, dd_dim, // double alpha, double *A, int lda
              q_tilde, dd_dim, //double *B, int ldb
              0, r_cov_inv, res_dim); //beta, double *C, int ldc
  // MAT_PRINTF(r_cov_inv, res_dim, res_dim);

  // dpotrf_(char *uplo, __CLPK_integer *n, __CLPK_doublereal *a, __CLPK_integer *
  //       lda, __CLPK_integer *info)
  //dpotri_(char *uplo, __CLPK_integer *n, __CLPK_doublereal *a, __CLPK_integer *
  //        lda, __CLPK_integer *info)
  char uplo = 'L'; //actually this makes it upper. the lapack stuff is all column major, but we're good since we're operiating on symmetric matrices
  integer info;
  dpotrf_(&uplo, &res_dim, r_cov_inv, &res_dim, &info);
  // printf("info: %i\n", (int) info);
  // MAT_PRINTF(r_cov_inv, res_dim, res_dim);
  dpotri_(&uplo, &res_dim, r_cov_inv, &res_dim, &info);
  for (u8 i=0; i < res_dim; i++) {
    for (u8 j=0; j < i; j++) {
      r_cov_inv[i*res_dim + j] = r_cov_inv[j*res_dim + i];
    }
  }
  // printf("info: %i\n", (int) info);
  // MAT_PRINTF(r_cov_inv, res_dim, res_dim);
}

void assign_r_vec(residual_mtxs_t *res_mtxs, u8 num_dds, double *dd_measurements, double *r_vec)
{
  cblas_dgemv(CblasRowMajor, CblasNoTrans,
              res_mtxs->null_space_dim, num_dds,
              1, res_mtxs->null_projector, num_dds,
              dd_measurements, 1,
              0, r_vec, 1);
  for (u8 i=0; i< num_dds; i++) {
    r_vec[i + res_mtxs->null_space_dim] = dd_measurements[i] - dd_measurements[i+num_dds] / GPS_L1_LAMBDA_NO_VAC;
  }
}

void assign_r_mean(residual_mtxs_t *res_mtxs, u8 num_dds, double *hypothesis, double *r_mean)
{
  cblas_dgemv(CblasRowMajor, CblasNoTrans,
              res_mtxs->null_space_dim, num_dds,
              1, res_mtxs->null_projector, num_dds,
              hypothesis, 1,
              0, r_mean, 1);
  memcpy(&r_mean[res_mtxs->null_space_dim], hypothesis, num_dds * sizeof(double));
}

double get_quadratic_term(residual_mtxs_t *res_mtxs, u8 num_dds, double *hypothesis, double *r_vec)
{
  // VEC_PRINTF(r_vec, res_mtxs->res_dim);
  double r[res_mtxs->res_dim];
  assign_r_mean(res_mtxs, num_dds, hypothesis, r);
  // VEC_PRINTF(r, res_mtxs->res_dim);
  for (u32 i=0; i<res_mtxs->res_dim; i++) {
    r[i] = r_vec[i] - r[i];
  }
  // VEC_PRINTF(r, res_mtxs->res_dim);
  double half_sig_dot_r[res_mtxs->res_dim];
  cblas_dsymv(CblasRowMajor, CblasUpper,
                 res_mtxs->res_dim,
                 1, res_mtxs->half_res_cov_inv, res_mtxs->res_dim,
                 r, 1,
                 0, half_sig_dot_r, 1);
  // VEC_PRINTF(half_sig_dot_r, res_mtxs->res_dim);
  double quad_term = 0;
  for (u32 i=0; i<res_mtxs->res_dim; i++) {
    quad_term -= half_sig_dot_r[i] * r[i];
  }
  return quad_term;
}

/** \} */
