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

#include <assert.h>
#include <clapack.h>
#include <inttypes.h>
#include <cblas.h>
#include <stdio.h>
#include <string.h>
#include "logging.h"
#include "ambiguity_test.h"
#include "common.h"
#include "constants.h"
#include "linear_algebra.h"
#include "single_diff.h"
#include "amb_kf.h"
#include "lambda.h"
#include "memory_pool.h"
#include "printing_utils.h"
#include "sats_management.h"

#define RAW_PHASE_BIAS_VAR 0
#define DECORRELATED_PHASE_BIAS_VAR 0
#define NUM_SEARCH_STDS 5
#define LOG_PROB_RAT_THRESHOLD -90
#define SINGLE_OBS_CHISQ_THRESHOLD 20

// TODO delete?
static void matrix_multiply_z_t(u32 n, u32 m, u32 p, const z_t *a,
                         const z_t *b, z_t *c)
{
  matrix_multiply_s64(n,m,p,a,b,c);
}

/** \defgroup ambiguity_test Integer Ambiguity Resolution
 * Integer ambiguity resolution using bayesian hypothesis testing.
 * \{ */
void create_empty_ambiguity_test(ambiguity_test_t *amb_test)
{
  static u8 pool_buff[MAX_HYPOTHESES*(sizeof(hypothesis_t) + sizeof(void *))];
  static memory_pool_t pool;
  amb_test->pool = &pool;
  memory_pool_init(amb_test->pool, MAX_HYPOTHESES, sizeof(hypothesis_t), pool_buff);

  amb_test->sats.num_sats = 0;
  amb_test->amb_check.initialized = 0;
}
void create_ambiguity_test(ambiguity_test_t *amb_test)
{
  create_empty_ambiguity_test(amb_test);

  /* Initialize pool with single element with num_dds = 0, i.e.
   * zero length N vector, i.e. no satellites. When we take the
   * product of this single element with the set of new satellites
   * we will just get a set of elements corresponding to the new sats. */
  hypothesis_t *empty_element = (hypothesis_t *)memory_pool_add(amb_test->pool);
  /* Start with ll = 0, just for the sake of argument. */
  empty_element->ll = 0;
}

void destroy_ambiguity_test(ambiguity_test_t *amb_test)
{
  memory_pool_destroy(amb_test->pool);
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
static void fold_contains(void *x, element_t *elem)
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

/** A struct to be used as the x in a memory pool fold, checking likelihood.
 * Used in fold_ll().
 */
typedef struct {
  u8 num_dds;             /**< The number of double differences in the hypotheses. */
  s32 N[MAX_CHANNELS-1];  /**< The ambiguity vector being searched for. */
  u8 found;               /**< Whether or not the ambiguity vector is found yet. */
  double ll;              /**< The pseudo log likelihood if it's been found. */
} fold_ll_t;

/** A memory pool fold method to check a hypothesis's pseudo log likelihood.
 *
 * Should be used in memory_pool_fold().
 *
 * \param x    The accumulator. Contains the hypothesis and whether it's found.
 * \param elem The hypothesis being folded against (checked for match).
 */
static void fold_ll(void *x, element_t *elem)
{
  fold_ll_t *acc = (fold_ll_t *) x;
  hypothesis_t *hyp = (hypothesis_t *) elem;
  if (acc->found == 1) {
    return;
  }
  for (u8 i=0; i<acc->num_dds; i++) {
    if (hyp->N[i] != acc->N[i]) {
      return;
    }
  }
  acc->ll = hyp->ll;
  acc->found = 1;
}

/** Finds the unnormalized log likelihood of an input ambiguity
 *
 * \param amb_test    The test to check against.
 * \param num_ambs    The length of the ambs vector.
 * \param ambs        The ambiguity hypothesis to look for.
 * \return            The pseudo log-likelihood of the hypothesis if it is
 *                    in the pool, and a positive number otherwise.
 */
double ambiguity_test_pool_ll(ambiguity_test_t *amb_test, u8 num_ambs, double *ambs)
{
  fold_ll_t acc;
  acc.num_dds = amb_test->sats.num_sats-1;
  acc.ll = 1;
  assert(acc.num_dds == num_ambs);
  for (u8 i=0; i<acc.num_dds; i++) {
    acc.N[i] = lround(ambs[i]);
  }
  memory_pool_fold(amb_test->pool, (void *) &acc, &fold_ll);
  return acc.ll;
}

/** A memory pool fold method to add probabilities.
 *
 * Should be used in memory_pool_fold().
 *
 * \param x    The accumulator. Contains the sum of the probabilities thus far.
 * \param elem The hypothesis being folded against (checked for match).
 */
static void fold_prob(void *x, element_t *elem)
{
  double *acc = (double *) x;
  hypothesis_t *hyp = (hypothesis_t *) elem;
  *acc += exp(hyp->ll);
}

/** Finds the probability of an input ambiguity
 *
 * \param amb_test    The test to check against.
 * \param num_ambs    The length of the ambs vector.
 * \param ambs        The ambiguity hypothesis to look for.
 * \return            The probability of the hypothesis if it is
 *                    in the pool, and a negative number otherwise.
 */
double ambiguity_test_pool_prob(ambiguity_test_t *amb_test, u8 num_ambs, double *ambs)
{
  fold_ll_t acc;
  acc.num_dds = amb_test->sats.num_sats-1;
  acc.ll = 1;
  assert(acc.num_dds == num_ambs);
  for (u8 i=0; i<acc.num_dds; i++) {
    acc.N[i] = lround(ambs[i]);
  }
  memory_pool_fold(amb_test->pool, (void *) &acc, &fold_ll);
  if (acc.ll > 0) {
    return -1;
  }
  double prob_sum = 0;
  memory_pool_fold(amb_test->pool, (void *) &prob_sum, &fold_prob);
  return exp(acc.ll) / prob_sum;
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
static void fold_mle(void *x, element_t *elem)
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
  mle.num_dds = CLAMP_DIFF(amb_test->sats.num_sats, 1);
  memory_pool_fold(amb_test->pool, (void *) &mle, &fold_mle);
  memcpy(ambs, mle.N, mle.num_dds * sizeof(s32));
}

/** Updates the IAR process with new measurements.
 *
 * Updates the satellites being tested, adding and removing hypotheses as needed.
 *
 * Assumes that the float estimates include all of the sats in the sdiffs.
 * Assumes that the dimension of the float estimates is one less than the number of sats.
 * Assumes that the sdiffs are arranged entirely in increasing order by PRN.
 *
 * \todo Return error codes?
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
  DEBUG_ENTRY();

  u8 num_sdiffs = state_dim + 1;

  if (amb_test->sats.num_sats < 5) {
    DEBUG_EXIT();
    return;
  }

  sdiff_t ambiguity_sdiffs[amb_test->sats.num_sats];
  double ambiguity_dd_measurements[2*(amb_test->sats.num_sats-1)];
  s8 valid_sdiffs = make_ambiguity_dd_measurements_and_sdiffs(
      amb_test, num_sdiffs, sdiffs, ambiguity_dd_measurements, ambiguity_sdiffs);

  /* Error */
  if (valid_sdiffs != 0) {
    log_error("update_ambiguity_test: invalid sdiffs. return code: %i\n", valid_sdiffs);
    DEBUG_EXIT();
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

  DEBUG_EXIT();
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
 * If a single observation was sufficiently unlikely to come from this hypothesis, we reject
 * the hypothesis. (In addition to the accumulated relative likelihood that is filtered upon later).
 *
 * To be given to memory_pool_fold().
 *
 * \param x     Points to a hyp_filter_t containing the accumulator and everything needed for the map.
 * \param elem  The hypothesis to be mapAccum'd.
 */
static s8 update_and_get_max_ll(void *x_, element_t *elem) {
  hyp_filter_t *x = (hyp_filter_t *) x_;
  hypothesis_t *hyp = (hypothesis_t *) elem;
  double hypothesis_N[x->num_dds];

  for (u8 i=0; i < x->num_dds; i++) {
    hypothesis_N[i] = hyp->N[i];
  }
  double q = get_quadratic_term(x->res_mtxs, x->num_dds, hypothesis_N, x->r_vec);
  hyp->ll += q;
  x->max_ll = MAX(x->max_ll, hyp->ll);
  return (abs(q) < SINGLE_OBS_CHISQ_THRESHOLD);
  /* Doesn't appear to need a dependence on d.o.f. to be effective.
   * We should revisit SINGLE_OBS_CHISQ_THRESHOLD when our noise model is tighter. */
}

/** Keeps track of which integer ambiguities are uninimously agreed upon in the pool.
 * \param num_dds   The number of DDs in each hypothesis. (Used to initialize amb_check).
 * \param hyp       The hypothesis to be checked against.
 * \param amb_check Keeps track of which ambs are still unanimous and their values.
 */
static void check_unanimous_ambs(u8 num_dds, hypothesis_t *hyp,
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
static s8 filter_and_renormalize(void *arg, element_t *elem) {
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
  DEBUG_ENTRY();

  hyp_filter_t x;
  x.num_dds = amb_test->sats.num_sats-1;
  assign_r_vec(&amb_test->res_mtxs, x.num_dds, dd_measurements, x.r_vec);
  x.max_ll = -1e20; // TODO get the first element, or use this as threshold to restart test
  x.res_mtxs = &amb_test->res_mtxs;
  x.unanimous_amb_check = &amb_test->amb_check;
  x.unanimous_amb_check->initialized = 0;

  memory_pool_filter(amb_test->pool, (void *) &x, &update_and_get_max_ll);
  memory_pool_filter(amb_test->pool, (void *) &x, &filter_and_renormalize);
  if (memory_pool_empty(amb_test->pool)) {
    log_debug("Ambiguity pool empty\n");
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
  if (DEBUG) {
    memory_pool_map(amb_test->pool, &x.num_dds, &print_hyp);
    printf("num_unanimous_ndxs=%u\n", x.unanimous_amb_check->num_matching_ndxs);
  }

  DEBUG_EXIT();
}

/* This says whether we can use the ambiguity_test to resolve a position in 3-space.
 */
u8 ambiguity_iar_can_solve(ambiguity_test_t *amb_test)
{
  return amb_test->amb_check.initialized &&
         amb_test->amb_check.num_matching_ndxs >= 3;
}

static bool is_prn_set(u8 len, u8 *prns)
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
  DEBUG_ENTRY();

  if (DEBUG) {
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
    log_error("There is disorder in the amb_test sats.\n");
    printf("amb_test sat prns = {%u, ", ref_prn);
    for (u8 k=0; k < num_dds; k++) {
      printf("%u, ", non_ref_prns[k]);
    }
    printf("}\n");
    DEBUG_EXIT();
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
      DEBUG_EXIT();
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
    DEBUG_EXIT();
    return -1;
  }
  for (i=0; i < num_dds; i++) {
    ambiguity_dd_measurements[i] -= ref_phase;
    ambiguity_dd_measurements[i+num_dds] -= ref_pseudorange;
  }
  if (DEBUG) {
    printf("amb_sdiff_prns = {");
    for (i = 0; i < num_dds+1; i++) {
      printf("%u, ", amb_sdiffs[i].prn);
    }
    printf("}\ndd_measurements = {");
    for (i=0; i < 2 * num_dds; i++) {
      printf("%f, \t", ambiguity_dd_measurements[i]);
    }
    printf("}\n");
  }

  DEBUG_EXIT();
  return 0;
}


/** Make the DD measurements and sdiffs that correspond to the resolved DDs in amb_test.
 * Assuming that amb_test has not been modified since the last check to see
 * whether it can resolve a fixed position (that is, amb_test->amb_check is up
 * to date, this will make a set of sdiffs and DD measurements that correspond
 * to the resolved sats.
 *
 * \todo If the input sdiffs are a subset of the resolved IAR sats, but still
 *       enough to compute a solution, do it.
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
  DEBUG_ENTRY();

  if (DEBUG) {
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
    log_debug("non_ref_prns[%u] = %u, \t (ndx=%u) \t amb[%u] = %"PRId32"\n",
              i, non_ref_prns[i], amb_test->amb_check.matching_ndxs[i],
              i, amb_test->amb_check.ambs[i]);
  }
  s8 valid_sdiffs =
    make_dd_measurements_and_sdiffs(ref_prn, non_ref_prns, num_dds, num_sdiffs,
                                    sdiffs, ambiguity_dd_measurements,
                                    amb_sdiffs);

  DEBUG_EXIT();
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
  DEBUG_ENTRY();

  u8 ref_prn = amb_test->sats.prns[0];
  u8 *non_ref_prns = &amb_test->sats.prns[1];
  u8 num_dds = CLAMP_DIFF(amb_test->sats.num_sats, 1);
  s8 valid_sdiffs = make_dd_measurements_and_sdiffs(ref_prn, non_ref_prns,
                                  num_dds, num_sdiffs, sdiffs,
                                  ambiguity_dd_measurements, amb_sdiffs);

  DEBUG_EXIT();
  return valid_sdiffs;
}


s8 sats_match(const ambiguity_test_t *amb_test, const u8 num_sdiffs, const sdiff_t *sdiffs)
{
  DEBUG_ENTRY();

  if (DEBUG) {
    printf("amb_test.sats.prns = {");
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
    log_debug("sats don't match (different length)\n");
    DEBUG_EXIT();
    return 0;
  }
  const u8 *prns = amb_test->sats.prns;
  const u8 amb_ref = amb_test->sats.prns[0];
  u8 j=0;
  for (u8 i = 1; i<amb_test->sats.num_sats; i++) { //TODO will not having a j condition cause le fault du seg?
    if (j >= num_sdiffs) {
      log_debug("sats don't match\n");
      DEBUG_EXIT();
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
      log_debug("sats don't match\n");
      DEBUG_EXIT();
      return 0;
    }
  }

  DEBUG_EXIT();
  return 1;
}

typedef struct {
  u8 num_sats;
  u8 old_prns[MAX_CHANNELS];
  u8 new_prns[MAX_CHANNELS];
} rebase_prns_t;

static void rebase_hypothesis(void *arg, element_t *elem) //TODO make it so it doesn't have to do all these lookups every time
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

/** Update an ambiguity test's reference satellite.
 * Given a set of sdiffs, choose a new reference that is hopefully already
 * tracked. If that's impossible, just choose a reference.
 * If the ambiguity test has the new reference, rebase the old hypotheses.
 * Otherwise trash the test and start over.
 * On return, the reference sat should be in the sdiffs.
 *
 * \param amb_test              The ambiguity test to update
 * \param num_sdiffs            The length of the sdiffs array.
 * \param sdiffs                sdiffs, sorted by PRN.
 * \param sdiffs_with_ref_first sdiffs, sorted by PRN after the first element,
 *                              which is the reference.
 * \return Whether the reference sat has changed.
 */
u8 ambiguity_update_reference(ambiguity_test_t *amb_test, const u8 num_sdiffs, const sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  DEBUG_ENTRY();

  u8 changed_ref = 0;
  u8 old_prns[amb_test->sats.num_sats];
  memcpy(old_prns, amb_test->sats.prns, amb_test->sats.num_sats * sizeof(u8));

  s8 sats_management_code = rebase_sats_management(&amb_test->sats, num_sdiffs, sdiffs, sdiffs_with_ref_first);
  if (sats_management_code != OLD_REF) {
    log_debug("updating iar reference sat\n");
    changed_ref = 1;
    if (sats_management_code == NEW_REF_START_OVER) {
      create_ambiguity_test(amb_test);
    }
    else {
      u8 new_prns[amb_test->sats.num_sats];
      memcpy(new_prns, amb_test->sats.prns, amb_test->sats.num_sats * sizeof(u8));

      rebase_prns_t prns = {.num_sats = amb_test->sats.num_sats};
      memcpy(prns.old_prns, old_prns, amb_test->sats.num_sats * sizeof(u8));
      memcpy(prns.new_prns, new_prns, amb_test->sats.num_sats * sizeof(u8));
      memory_pool_map(amb_test->pool, &prns, &rebase_hypothesis);
    }
  }

  DEBUG_EXIT();
  return changed_ref;
}

typedef struct {
  u8 num_ndxs;
  u8 intersection_ndxs[MAX_CHANNELS - 1];
} intersection_ndxs_t;

static s32 projection_comparator(void *arg, element_t *a, element_t *b)
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

static void projection_aggregator(element_t *new_, void *x_, u32 n, element_t *elem_)
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

u8 ambiguity_sat_projection(ambiguity_test_t *amb_test, const u8 num_dds_in_intersection, const u8 *dd_intersection_ndxs)
{
  DEBUG_ENTRY();

  u8 num_dds_before_proj = CLAMP_DIFF(amb_test->sats.num_sats, 1);
  if (num_dds_before_proj == num_dds_in_intersection) {
    log_debug("no need for projection\n");
    DEBUG_EXIT();
    return 0;
  }

  intersection_ndxs_t intersection = {.num_ndxs = num_dds_in_intersection};
  memcpy(intersection.intersection_ndxs, dd_intersection_ndxs, num_dds_in_intersection * sizeof(u8));


  log_info("IAR: %"PRIu32" hypotheses before projection\n", memory_pool_n_allocated(amb_test->pool));
  memory_pool_group_by(amb_test->pool,
                       &intersection, &projection_comparator,
                       &intersection, sizeof(intersection),
                       &projection_aggregator);
  log_info("IAR: updates to %"PRIu32"\n", memory_pool_n_allocated(amb_test->pool));
  log_info("After projection, num_sats = %d\n", num_dds_in_intersection + 1);
  u8 work_prns[MAX_CHANNELS];
  memcpy(work_prns, amb_test->sats.prns, amb_test->sats.num_sats * sizeof(u8));
  for (u8 i=0; i<num_dds_in_intersection; i++) {
    amb_test->sats.prns[i+1] = work_prns[dd_intersection_ndxs[i]+1];
  }
  amb_test->sats.num_sats = num_dds_in_intersection+1;

  DEBUG_EXIT();
  return 1;
}

static void vec_plus(u8 cols, u8 rows, z_t *v, z_t *Z, z_t mult, u8 column)
{
  for(u8 i = 0; i < rows; i++) {
    v[i] += Z[i * cols + column] * mult;
  }
}
static s8 increment_matrix_product(u8 len, z_t *counter, u8 vlen, z_t *Z,
                                   z_t *v, z_t *lower_bounds, z_t *upper_bounds) {
  if (memcmp(upper_bounds, counter, len * sizeof(z_t)) == 0) {
    /* counter has reached upper_bound, terminate iteration. */
    return 0;
  }
  for (u8 i=0; i < len; i++) {
    u8 column = i;
    counter[i]++;
    vec_plus(len, vlen, v, Z, 1, column);
    if (counter[i] > upper_bounds[i]) {
      /* counter[i] has reached maximum, reset counter[i]
       * to lower[i] and 'carry' to next 'digit' */
      vec_plus(len, vlen, v, Z, -counter[i] + lower_bounds[i], column);
      counter[i] = lower_bounds[i];
    } else {
      /* Incremented, so now we have the next counter value. */
      break;
    }
  }
  return 1;
}

static bool inside(u32 dim, z_t *point, z_t *lower_bounds, z_t *upper_bounds)
{
  /* Without a slight tolerance, some cases have strange behavior because
   * points on the boundary may not be included; e.g. even with an identity Z
   * matrix some points might not be yielded */
  for (u8 i = 0; i < dim; i++) {
    if (point[i] < lower_bounds[i] || point[i] > upper_bounds[i]) {
      return false;
    }
  }
  return true;
}

/* Initializes x->zimage */
static void init_intersection_count_vector(intersection_count_t *x, hypothesis_t *hyp)
{
  u8 full_dim = x->old_dim + x->new_dim;
  /* Initialize counter using lower bounds. */
  memcpy(x->counter, x->itr_lower_bounds, x->new_dim * sizeof(z_t));
  z_t v0[full_dim];
  /* Map the lower bound vector using Z2_inverse into the second half of v0. */
  matrix_multiply_z_t(x->new_dim, x->new_dim, 1, x->Z2_inv, x->counter, v0 + x->old_dim);
  /* Map the old hypothesis values identically into the first half of v0. */
  for (u8 i = 0; i < x->old_dim; i++) {
    v0[i] = hyp->N[i];
  }
  /* Decorrelate the joint vector. */
  matrix_multiply_z_t(full_dim, full_dim, 1, x->Z1, v0, x->zimage);
}

static void fold_intersection_count(void *arg, element_t *elem)
{
  intersection_count_t *x = (intersection_count_t *) arg;
  hypothesis_t *hyp = (hypothesis_t *)elem;
  u8 full_dim = x->old_dim + x->new_dim;

  /* Set initial image in decorrelated space */
  init_intersection_count_vector(x, hyp);

  do {
    if (inside(full_dim, x->zimage, x->box_lower_bounds, x->box_upper_bounds)) {
      x->intersection_size++;
    }
  } while (0 != increment_matrix_product(
                  x->new_dim, x->counter,
                  full_dim, x->Z, x->zimage,
                  x->itr_lower_bounds, x->itr_upper_bounds));
}

static void round_matrix(u32 rows, u32 cols, const double *A, z_t *B)
{
  for (u8 i=0; i < rows; i++) {
    for (u8 j=0; j < cols; j++) {
      B[i * cols + j] = lround(A[i * cols + j]);
    }
  }
}

/* Computes Z1' * Z2^-1, where Z1' is the rightmost (new_dim) columns of Z1. */
/* Assumes Z1 is written in a basis where the old dds come first. */
static void compute_Z(u8 old_dim, u8 new_dim, const z_t *Z1, const z_t * Z2_inv, z_t *transform)
{
  u8 full_dim = old_dim + new_dim;
  z_t Z1_right[full_dim * new_dim];
  /* Take the right columns of Z1 */
  for (u8 i = 0; i < full_dim; i++) {
    memcpy(&Z1_right[i*new_dim], &Z1[i*full_dim + old_dim], new_dim * sizeof(z_t));
  }
  matrix_multiply_z_t(full_dim, new_dim, new_dim, Z1_right, Z2_inv, transform);
}

/* TODO(dsk) Use submatrix for this instead? */
static void remap_prns(ambiguity_test_t *amb_test, u8 ref_prn,
                       u32 num_added_dds, u8 *added_prns,
                       generate_hypothesis_state_t2 *s)
{
  intersection_count_t *x = s->x;
  u8 i = 0;
  u8 j = 0;
  u8 k = 0;
  u8 old_prns[x->old_dim];
  memcpy(old_prns, &amb_test->sats.prns[1], x->old_dim * sizeof(u8));
  while (k < x->old_dim + num_added_dds) {
    if (j == x->new_dim || (old_prns[i] < added_prns[j] && i != x->old_dim)) {
      s->ndxs_of_old_in_new[i] = k;
      amb_test->sats.prns[k+1] = old_prns[i];
      i++;
      k++;
    } else if (i == x->old_dim || old_prns[i] > added_prns[j]) {
      s->ndxs_of_added_in_new[j] = k;
      amb_test->sats.prns[k+1] = added_prns[j];
      j++;
      k++;
    } else {
      log_error("remap_prns: impossible condition reached.\n");
      printf("old_prns = [");
      for (u8 ii=0; ii < x->old_dim; ii++) {
        printf("%d, ",old_prns[ii]);
      }
      printf("]\n");
      printf("added_prns = [");
      for (u8 jj=0; jj < x->old_dim; jj++) {
        printf("%d, ", added_prns[jj]);
      }
      printf("]\n");
      break;
    }
  }
  amb_test->sats.prns[0] = ref_prn;
  amb_test->sats.num_sats = k+1;
}

static s8 intersection_generate_next_hypothesis0(void *x_, u32 n)
{
  (void) n;
  generate_hypothesis_state_t2 *g = (generate_hypothesis_state_t2 *) x_;
  intersection_count_t *x = g->x;
  u8 full_dim = x->old_dim + x->new_dim;

  do {
    if (inside(full_dim, x->zimage, x->box_lower_bounds, x->box_upper_bounds)) {
      /* Yield current point. */
      return 1;
    }
  } while (0 != increment_matrix_product(
                  x->new_dim, x->counter,
                  full_dim, x->Z, x->zimage,
                  x->itr_lower_bounds, x->itr_upper_bounds));
  return 0;
}

/* Increment the iterator, then continue 0 or more times until it's valid. */
static s8 intersection_generate_next_hypothesis1(void *x_, u32 n)
{
  generate_hypothesis_state_t2 *g = (generate_hypothesis_state_t2 *) x_;
  intersection_count_t *x = g->x;
  u8 full_dim = x->old_dim + x->new_dim;

  if (0 == increment_matrix_product(x->new_dim, x->counter, full_dim, x->Z, x->zimage,
             x->itr_lower_bounds, x->itr_upper_bounds)) {
    return 0;
  }

  return intersection_generate_next_hypothesis0(x_, n);
}

static s8 intersection_init(void *x, element_t *elem)
{
  generate_hypothesis_state_t2 *g = (generate_hypothesis_state_t2 *) x;
  hypothesis_t *hyp = (hypothesis_t *)elem;

  init_intersection_count_vector(g->x, hyp);
  /* Find a valid first point. */
  return intersection_generate_next_hypothesis0(x, 0);
}

static void intersection_hypothesis_prod(element_t *new_, void *x_, u32 n, element_t *elem_)
{
  (void) elem_, (void) n;
  generate_hypothesis_state_t2 *s = (generate_hypothesis_state_t2 *) x_;
  intersection_count_t *x = s->x;
  hypothesis_t *new = (hypothesis_t *) new_;
  u8 *ndxs_of_old_in_new   = s->ndxs_of_old_in_new;
  u8 *ndxs_of_added_in_new = s->ndxs_of_added_in_new;

  s32 old_N[MAX_CHANNELS-1];
  memcpy(old_N, new->N, x->old_dim * sizeof(s32));

  for (u8 i=0; i < x->old_dim; i++) {
    new->N[ndxs_of_old_in_new[i]] = old_N[i];
  }
  for (u8 i=0; i < x->new_dim; i++) {
    new->N[ndxs_of_added_in_new[i]] = 0;
    for (u8 j=0; j < x->new_dim; j++) {
      new->N[ndxs_of_added_in_new[i]] += s->Z_new_inv[i*x->new_dim + j] * x->counter[j];
    }
  }
}

static s32 add_sats(ambiguity_test_t *amb_test,
                    u8 ref_prn, u8 *added_prns,
                    intersection_count_t *x)
{
  generate_hypothesis_state_t2 s;
  s.x = x;
  s.Z_new_inv = x->Z2_inv;
  remap_prns(amb_test, ref_prn, x->new_dim, added_prns, &s);
  s32 count = memory_pool_product_generator(amb_test->pool, &s, MAX_HYPOTHESES, sizeof(s),
                  &intersection_init,
                  &intersection_generate_next_hypothesis1,
                  &intersection_hypothesis_prod);
  (void) count;
  s32 num_hyps = memory_pool_n_allocated(amb_test->pool);
  log_info("IAR: updates to %"PRIu32"\n", num_hyps);
  log_info("add_sats. num sats: %i\n", amb_test->sats.num_sats);
  return num_hyps;
}


/*
 * The satellite inclusion algorithm considers three important vector spaces:
 *  - The correlated space of integer ambiguities considered by the float filter (V0)
 *  - The decorrelated space of all integer ambiguities (V1)
 *  - The decorrelated space of integer ambiguities for only those satellites
 *    not currently considered by the IAR hypothesis test (V2). These sats are
 *    called "new", and the others are called "old."
 *
 *  A decorrelation matrix exists for mapping V0 to V1 (called Z1) and a
 *  similar one exists for mapping the subspace of new ambiguities in V0 to V2
 *  (called Z2).
 *
 *  These decorrelation matrices determine a certain range of likely
 *  ambiguities, lying inside a box in V2 (likely values for new sats) and V1
 *  (likely values for all sats considered jointly). The V1 box considers all
 *  joint correlations and gives a more accurate range of possibilities.
 *  However, we do not want to add any hypotheses which conflict on old
 *  satellites with the existing hypothesis set, so we iterate through the
 *  current hypothesis set, combining each hypothesis with a hypothesis from
 *  the V2 (new sat) box, and map it back to V1, keeping only those that lie
 *  inside the V1 box.
 *
 *  If too many hypotheses result to fit in memory, we repeat the calculation
 *  using fewer new sats. If it is impossible to add a sufficient number of
 *  sats to make progress towards an RTK solution (< 4 double differences
 *  total) we return without adding any.
 */
static u8 inclusion_loop_body(
       u8 num_dds_to_add,
       memory_pool_t *pool, u8 state_dim, u8 num_addible_dds,
       const double *ordered_N_cov, const double *ordered_N_mean,
       const double *addible_cov, const double *addible_mean,
       intersection_count_t *x, u32 *full_size_return)
{
  x->new_dim = num_dds_to_add;
  s32 current_num_hyps = memory_pool_n_allocated(pool);
  u32 max_num_hyps = memory_pool_n_elements(pool);

  u8 num_current_dds = x->old_dim;
  u8 full_dim = num_current_dds + num_dds_to_add;

  /* TODO(dsk) tune this constant.
   * This determines how many hypotheses will be examined by the intersection
   * memory_pool fold below. */
  u32 max_iteration_size = 10000;

  /* Calculate the two decorrelation matrices and their related matrices. */
  u32 full_size =
    float_to_decor(ordered_N_cov, ordered_N_mean,
        state_dim, full_dim,
        x->box_lower_bounds, x->box_upper_bounds, x->Z1, x->Z1_inv);


  /* Useful for debugging. */
  *full_size_return = full_size;

  u32 box_size =
    float_to_decor(addible_cov, addible_mean,
      num_addible_dds, num_dds_to_add,
      x->itr_lower_bounds, x->itr_upper_bounds, x->Z2, x->Z2_inv);

  compute_Z(num_current_dds, num_dds_to_add, x->Z1, x->Z2_inv, x->Z);

  if (full_size <= max_num_hyps) {
    log_debug("BRANCH 1: num dds: %i. full size: %"PRIu32", itr size: %"PRIu32"\n", num_dds_to_add, full_size, box_size);
    /* The hypotheses generated for these double-differences fit. */
    return 1;
  } else if (box_size * current_num_hyps <= max_iteration_size) {
    log_debug("BRANCH 2: num dds: %i. full size: %"PRIu32", itr size: %"PRIu32"\n", num_dds_to_add, full_size, box_size);

    x->intersection_size = 0;

    /* Do intersection */
    memory_pool_fold(pool, (void *)x, &fold_intersection_count);

    log_debug("intersection size: %i\n", x->intersection_size);

    if (x->intersection_size < max_num_hyps) {
      /* The hypotheses generated for these double-differences fit. */
      return 1;
    }
  }

  /* Can't add sats for this value of num_dds_to_add. */
  return 0;
}

/** Perform the inclusion step of adding satellites to the amb test (as we can).
 * Using the Kalman filter's state estimates, add new satellites to the
 * hypotheses being tracked. There's a chance the KF has drifted away from
 * EVERYTHING in the hypothesis pool, so then we have to return an indication
 * of that.
 *
 * \param amb_test                The amb_test struct whose sats we are updating.
 * \param num_dds_in_intersection The number of DD measurements common between
 *                                the float filter and the amb_test last timestep.
 * \param float_sats              The sats for the KF.
 * \param float_mean              The KF estimate
 * \param float_cov_U             The KF covariance U (from UDU decomposition)
 * \param float_cov_D             The KF covariance D (from UDU decomposition)
 * \returns 0 if we didn't change amb_test's sats
 *          1 if we changed the sats, but don't need to start over.
 *          2 if we need to start over (e.g. we have no hypotheses left).
 */
u8 ambiguity_sat_inclusion(ambiguity_test_t *amb_test, const u8 num_dds_in_intersection,
                           const sats_management_t *float_sats, const double *float_mean,
                           const double *float_cov_U, const double *float_cov_D)
{
  if (float_sats->num_sats <= num_dds_in_intersection + 1 || float_sats->num_sats < 5) {
    /* Nothing added if we alread have all the sats or the KF has too few sats
     * such that we couldn't test anyways. Changing the < 5 can allow code to
     * run below that needs to add to no less than 4 DD's.  */
    return 0;
  }

  u8 state_dim = float_sats->num_sats-1;
  double float_cov[state_dim * state_dim];
  u8 float_prns[float_sats->num_sats];
  double N_mean[state_dim];

  matrix_reconstruct_udu(state_dim, float_cov_U, float_cov_D, float_cov);
  memcpy(float_prns, float_sats->prns, float_sats->num_sats * sizeof(u8));
  memcpy(N_mean, float_mean, (float_sats->num_sats-1) * sizeof(double));

  /* After this block, float_prns will have the correct reference,
   * as will N_cov and N_mean */
  if (amb_test->sats.num_sats >= 2 &&
      amb_test->sats.prns[0] != float_sats->prns[0]) {
    u8 old_prns[float_sats->num_sats];
    memcpy(old_prns, float_sats->prns, float_sats->num_sats * sizeof(u8));
    set_reference_sat_of_prns(amb_test->sats.prns[0], float_sats->num_sats, float_prns);
    rebase_mean_N(N_mean, float_sats->num_sats, old_prns, float_prns);
    rebase_covariance_sigma(float_cov, float_sats->num_sats, old_prns, float_prns);
  }
  u8 ref_prn = float_prns[0];

  double N_cov[state_dim * state_dim];
  memcpy(N_cov, float_cov, state_dim * state_dim * sizeof(double));

  /* Find the locations of new prns and old prns so we can reorder our matrices. */
  u8 i = 1;
  u8 j = 1;
  u8 num_addible_dds = 0;
  u32 ndxs_of_new_dds_in_float[MAX_CHANNELS-1];
  u8 num_old_dds = 0;
  u32 ndxs_of_old_dds_in_float[MAX_CHANNELS-1];
  u8 new_dd_prns[MAX_CHANNELS-1];
  while (j < float_sats->num_sats) {
    if (i < amb_test->sats.num_sats && amb_test->sats.prns[i] == float_prns[j]) {
      ndxs_of_old_dds_in_float[num_old_dds++] = j-1;
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
  u32 row_map[1] = {0};
  /* Take just the new dds. */
  submatrix(num_addible_dds, num_addible_dds, state_dim, N_cov,
      ndxs_of_new_dds_in_float, ndxs_of_new_dds_in_float, addible_float_cov);
  submatrix(1, num_addible_dds, state_dim, N_mean,
      row_map, ndxs_of_new_dds_in_float, addible_float_mean);

  u8 num_current_dds = MAX(amb_test->sats.num_sats, 1) - 1;

  assert(num_current_dds == num_old_dds);
  assert(state_dim == num_old_dds + num_addible_dds);

  u8 min_dds_to_add = MAX(1, 4 - num_current_dds);

  if (min_dds_to_add > num_addible_dds) {
    return 0;
  }

  /* Reorder the covariance matrix basis so that old sats come first: */
  /* [ old_sats | new_sats ] */
  u32 reordering[state_dim];
  memcpy(reordering, ndxs_of_old_dds_in_float, num_old_dds * sizeof(u32));
  memcpy(reordering + num_old_dds, ndxs_of_new_dds_in_float,
      num_addible_dds * sizeof(u32));

  /* Rearrange N_cov, N_mean according to reordering. */
  double N_cov_ordered[state_dim * state_dim];
  double N_mean_ordered[state_dim];
  submatrix(state_dim, state_dim, state_dim, N_cov,
      reordering, reordering, N_cov_ordered);
  submatrix(1, state_dim, state_dim, N_mean,
      row_map, reordering, N_mean_ordered);

  /* Initialize intersection struct. */
  z_t counter[num_addible_dds];
  z_t lower_bounds1[state_dim];
  z_t upper_bounds1[state_dim];
  z_t lower_bounds2[num_addible_dds];
  z_t upper_bounds2[num_addible_dds];
  z_t zimage[state_dim];
  z_t Z1[state_dim * state_dim];
  z_t Z1_inv[state_dim * state_dim];
  z_t Z2[num_addible_dds * num_addible_dds];
  z_t Z2_inv[num_addible_dds * num_addible_dds];
  z_t Z1_Z2_inv[state_dim * num_addible_dds];

  intersection_count_t x;
  x.new_dim = num_addible_dds;
  x.old_dim = num_current_dds;
  x.counter = counter;
  // TODO(dsk) clarify name in struct
  x.box_lower_bounds = lower_bounds1;
  x.box_upper_bounds = upper_bounds1;
  x.itr_lower_bounds = lower_bounds2;
  x.itr_upper_bounds = upper_bounds2;
  x.zimage = zimage;
  x.Z = Z1_Z2_inv;
  x.Z1 = Z1;
  x.Z2 = Z2;
  x.Z1_inv = Z1_inv;
  x.Z2_inv = Z2_inv;

  u32 full_size = 0;

  /* Check to see if min_dds_to_add will not fit. If so, don't bother
   * iterating through all the sats below. */
  u8 fits = inclusion_loop_body(
      min_dds_to_add, amb_test->pool, state_dim, num_addible_dds,
      N_cov_ordered, N_mean_ordered, addible_float_cov, addible_float_mean,
      &x, &full_size);
  if (fits == 0) {
    return 0;
  }

  /* Try to add as many dds to the IAR as possible; return 0 if no amount fits. */
  for (u8 num_dds_to_add = num_addible_dds;
       num_dds_to_add >= min_dds_to_add;
       num_dds_to_add--)
  {
    u8 fits = inclusion_loop_body(
        num_dds_to_add, amb_test->pool, state_dim, num_addible_dds,
        N_cov_ordered, N_mean_ordered, addible_float_cov, addible_float_mean,
        &x, &full_size);

    if (fits == 1) {
      /* Sats should be added. The struct x contains new_dim, the correct
       * number to add, along with the matrices needed to do so . */
      s32 num_hyps = add_sats(amb_test, ref_prn, new_dd_prns, &x);
      if (num_hyps == 0) {
        return 2;
      } else {
        return 1;
      }
    }
  }
  log_debug("BRANCH 3: covariance too large. full: %"PRIu32"\n", full_size);
  /* Covariance too large, nothing added. */
  return 0;
}

/* TODO(dsk) remove dead code. */
u8 ambiguity_sat_inclusion_old(ambiguity_test_t *amb_test, u8 num_dds_in_intersection,
                               sats_management_t *float_sats, double *float_mean,
                               double *float_cov_U, double *float_cov_D)
{
  DEBUG_ENTRY();

  if (float_sats->num_sats <= num_dds_in_intersection + 1 || float_sats->num_sats < 2) {
    log_debug("no need for inclusion\n");
    DEBUG_EXIT();
    return 0;
  }

  u32 state_dim = float_sats->num_sats-1;
  double float_cov[state_dim * state_dim];
  matrix_reconstruct_udu(state_dim, float_cov_U, float_cov_D, float_cov);
  u8 float_prns[float_sats->num_sats];
  memcpy(float_prns, float_sats->prns, float_sats->num_sats * sizeof(u8));
  double N_mean[float_sats->num_sats-1];
  memcpy(N_mean, float_mean, (float_sats->num_sats-1) * sizeof(double));
  /* After this block, float_prns will have the correct reference,
   * as will N_cov and N_mean */
  if (amb_test->sats.num_sats >= 2 && amb_test->sats.prns[0] != float_sats->prns[0]) {
    u8 old_prns[float_sats->num_sats];
    memcpy(old_prns, float_sats->prns, float_sats->num_sats * sizeof(u8));
    set_reference_sat_of_prns(amb_test->sats.prns[0], float_sats->num_sats, float_prns);
    rebase_mean_N(N_mean, float_sats->num_sats, old_prns, float_prns);
    rebase_covariance_sigma(float_cov, float_sats->num_sats, old_prns, float_prns);
  }
  double N_cov[(float_sats->num_sats-1) * (float_sats->num_sats-1)];
  memcpy(N_cov, float_cov, state_dim * state_dim * sizeof(double)); //TODO we can just use N_cov throughout

  /* First get all the prns we might add and their covariances. */
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
      addible_float_cov[i*num_addible_dds + j] =
        N_cov[ndxs_of_new_dds_in_float[i]*(float_sats->num_sats-1) + ndxs_of_new_dds_in_float[j]];
    }
    addible_float_mean[i] = N_mean[ndxs_of_new_dds_in_float[i]];
  }

  /* Find the largest set of sats we can add */
  z_t Z_inv[num_addible_dds * num_addible_dds];
  z_t lower_bounds[num_addible_dds];
  z_t upper_bounds[num_addible_dds];
  u8 num_dds_to_add;
  s8 add_any_sats = determine_sats_addition(amb_test,
                                            addible_float_cov, num_addible_dds, addible_float_mean,
                                            lower_bounds, upper_bounds, &num_dds_to_add,
                                            Z_inv);
  if (add_any_sats == 1) {
    add_sats_old(amb_test, float_prns[0], num_dds_to_add, new_dd_prns, lower_bounds, upper_bounds, Z_inv);
    log_debug("adding sats\n");
    DEBUG_EXIT();
    return 1;
  } else {
    log_debug("not adding sats\n");
    DEBUG_EXIT();
    return 0;
  }
}

z_t float_to_decor(const double *addible_float_cov,
                   const double *addible_float_mean,
                   u8 num_addible_dds,
                   u8 num_dds_to_add,
                   z_t *lower_bounds, z_t *upper_bounds,
                   z_t *Z, z_t *Z_inv)
{
  u8 dim = num_dds_to_add;
  double Z_[dim * dim];
  double Z_inv_[dim * dim];

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

  lambda_reduction(num_dds_to_add, added_float_cov, Z_);

  double decor_float_cov_diag[num_dds_to_add];

  memset(decor_float_cov_diag, 0, num_dds_to_add * sizeof(double));

  for (u8 i=0; i < num_dds_to_add; i++) {
    for (u8 j=0; j < num_dds_to_add; j++) {
      for (u8 k=0; k < num_dds_to_add; k++) {
        decor_float_cov_diag[i] += Z_[i*num_dds_to_add + j] * added_float_cov[j * num_dds_to_add + k] * Z_[i*num_dds_to_add + k];
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
      decor_float_mean[i] += Z_[i*num_dds_to_add + j] * addible_float_mean[j];
    }
  }

  u64 new_hyp_set_cardinality = 1;
  for (u8 i=0; i<num_dds_to_add; i++) {
    double search_distance = NUM_SEARCH_STDS * sqrt(decor_float_cov_diag[i]);
    upper_bounds[i] = lround(ceil(decor_float_mean[i] + search_distance));
    lower_bounds[i] = lround(floor(decor_float_mean[i] - search_distance));
    new_hyp_set_cardinality *= upper_bounds[i] - lower_bounds[i] + 1;
  }

  if (Z_inv) {
    round_matrix(dim, dim, Z_, Z);
    matrix_inverse(dim, Z_, Z_inv_);
    round_matrix(dim, dim, Z_inv_, Z_inv);
  }

  return new_hyp_set_cardinality;
}

/* TODO(dsk) remove this function. */
s8 determine_sats_addition(ambiguity_test_t *amb_test,
                           double *float_N_cov, u8 num_float_dds, double *float_N_mean,
                           z_t *lower_bounds, z_t *upper_bounds, u8 *num_dds_to_add,
                           z_t *Z_inv)
{
  u8 num_current_dds = CLAMP_DIFF(amb_test->sats.num_sats, 1);
  /* num_current_dds + min_dds_to_add = 4,
   * so that we have a nullspace projector. */
  u8 min_dds_to_add = MAX(1, 4 - num_current_dds);

  u32 max_new_hyps_cardinality;
  s32 current_num_hyps = memory_pool_n_allocated(amb_test->pool);
  u32 max_num_hyps = memory_pool_n_elements(amb_test->pool);
  if (current_num_hyps <= 0) {
    max_new_hyps_cardinality = max_num_hyps;
  } else {
    max_new_hyps_cardinality = max_num_hyps / current_num_hyps;
  }

  *num_dds_to_add = num_float_dds;
  z_t Z[num_float_dds * num_float_dds];
  while (*num_dds_to_add >= min_dds_to_add) {
    u32 new_hyp_set_cardinality = float_to_decor(float_N_cov,
                                                 float_N_mean,
                                                 num_float_dds,
                                                 *num_dds_to_add,
                                                 lower_bounds, upper_bounds, Z, Z_inv);
    if (new_hyp_set_cardinality <= max_new_hyps_cardinality) {
      return 1;
    }
    else {
      *num_dds_to_add -= 1;
    }
  }
  return -1;
}

/** Add/drop satellites from the ambiguity test, changing reference if needed.
 * INVALIDATES unanimous ambiguities
 * ^ TODO record this in the amb_test state?
 * \param amb_test    An ambiguity test whose tests to update.
 * \param num_sdiffs  The length of the sdiffs array.
 * \param sdiffs      The single differenced observations. Sorted by PRN.
 * \param float_sats  The satellites to correspond to the KF mean and cov.
 * \param float_mean  The KF state estimate.
 * \param float_cov_U The KF state estimate covariance U in UDU decomposiiton.
 * \param float_cov_D The KF state estimate covariance D in UDU decomposiiton.
 * \return  0 if we didn't change the sats
 *          1 if we did change the sats
 *          2 if we need to reset IAR TODO maybe do that in here?
 */
u8 ambiguity_update_sats(ambiguity_test_t *amb_test, const u8 num_sdiffs,
                         const sdiff_t *sdiffs, const sats_management_t *float_sats,
                         const double *float_mean, const double *float_cov_U,
                         const double *float_cov_D)
{
  DEBUG_ENTRY();

  if (num_sdiffs < 2) {
    create_ambiguity_test(amb_test);
    log_debug("< 2 sdiffs, starting over\n");
    DEBUG_EXIT();
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
      create_ambiguity_test(amb_test);
    }

    // u8 num_dds_in_intersection = ambiguity_order_sdiffs_with_intersection(amb_test, sdiffs, float_cov, intersection_ndxs);
    if (ambiguity_sat_projection(amb_test, num_dds_in_intersection, intersection_ndxs)) {
      changed_sats = 1;
    }
    u8 incl = ambiguity_sat_inclusion(amb_test, num_dds_in_intersection,
                float_sats, float_mean, float_cov_U, float_cov_D);
    if (incl == 2) {
      create_ambiguity_test(amb_test);
      changed_sats = 1;
    } else if (incl == 1) {
      changed_sats = 1;
    }
  }
  log_debug("changed_sats = %u\n", changed_sats);

  DEBUG_EXIT();
  return changed_sats;
  /* TODO: Should we order by 'goodness'? Perhaps by volume of hyps? */
}

u8 find_indices_of_intersection_sats(const ambiguity_test_t *amb_test, const u8 num_sdiffs, const sdiff_t *sdiffs_with_ref_first, u8 *intersection_ndxs)
{
  DEBUG_ENTRY();

  if (DEBUG) {
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
      log_debug("(%u, \t%u, \t%u, \t%u, \t%u)\t\t\tamb_test->sats.prns[i] == sdiffs_with_ref_first[j].prn; i,j,k++\n",
                i, j, k, amb_test->sats.prns[i], sdiffs_with_ref_first[j].prn);
      intersection_ndxs[k] = i-1;
      i++;
      j++;
      k++;
    }
    else if (amb_test->sats.prns[i] < sdiffs_with_ref_first[j].prn) {
      log_debug("(%u, \t%u, \t%u, \t%u, \t%u)\t\t\tamb_test->sats.prns[i] <  sdiffs_with_ref_first[j].prn; i++\n",
                i, j, k, amb_test->sats.prns[i], sdiffs_with_ref_first[j].prn);
      i++;
    }
    else {
      log_debug("(%u, \t%u, \t%u, \t%u, \t%u)\t\t\tamb_test->sats.prns[i] >  sdiffs_with_ref_first[j].prn; j++\n",
                i, j, k, amb_test->sats.prns[i], sdiffs_with_ref_first[j].prn);
      j++;
    }
  }

  DEBUG_EXIT();
  return k;
}

/* TODO(dsk) remove dead code. */
typedef struct {
  z_t upper_bounds[MAX_CHANNELS-1];
  z_t lower_bounds[MAX_CHANNELS-1];
  z_t counter[MAX_CHANNELS-1];
  u8 ndxs_of_old_in_new[MAX_CHANNELS-1];
  u8 ndxs_of_added_in_new[MAX_CHANNELS-1];
  u8 num_added_dds;
  u8 num_old_dds;
  z_t Z_inv[(MAX_CHANNELS-1) * (MAX_CHANNELS-1)];
} generate_hypothesis_state_t;

static s8 generate_next_hypothesis(void *x_, u32 n)
{
  (void) n;
  generate_hypothesis_state_t *x = (generate_hypothesis_state_t *)x_;

  if (memcmp(x->upper_bounds, x->counter, x->num_added_dds * sizeof(z_t)) == 0) {
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

  return 1;
}

static void hypothesis_prod(element_t *new_, void *x_, u32 n, element_t *elem_)
{
  (void) elem_;
  (void) n;
  generate_hypothesis_state_t *x = (generate_hypothesis_state_t *)x_;
  hypothesis_t *new = (hypothesis_t *) new_;

  u8 *ndxs_of_old_in_new = x->ndxs_of_old_in_new;
  u8 *ndxs_of_added_in_new = x->ndxs_of_added_in_new;

  s32 old_N[MAX_CHANNELS-1];
  memcpy(old_N, new->N, x->num_old_dds * sizeof(z_t));

  for (u8 i=0; i < x->num_old_dds; i++) {
    new->N[ndxs_of_old_in_new[i]] = old_N[i];
  }
  for (u8 i=0; i<x->num_added_dds; i++) {
    new->N[ndxs_of_added_in_new[i]] = 0;
    for (u8 j=0; j<x->num_added_dds; j++) {
      new->N[ndxs_of_added_in_new[i]] += x->Z_inv[i*x->num_added_dds + j] * x->counter[j];
    }
  }

  /* NOTE: new->ll remains the same as elem->ll as p := exp(ll) is invariant under a
   * constant multiplicative factor common to all hypotheses. TODO: reference^2 document (currently lives in page 3/5.6/2014 of ian's notebook) */
}

/* TODO(dsk) remove dead code. */
static s8 no_init(void *x, element_t *elem) {
  (void) x; (void) elem;
  return 1;
}
void add_sats_old(ambiguity_test_t *amb_test,
                  u8 ref_prn,
                  u32 num_added_dds, u8 *added_prns,
                  z_t *lower_bounds, z_t *upper_bounds,
                  z_t *Z_inv)
{
  /* Make a generator that iterates over the new hypotheses. */
  generate_hypothesis_state_t x0;
  memcpy(x0.upper_bounds, upper_bounds, num_added_dds * sizeof(z_t));
  memcpy(x0.lower_bounds, lower_bounds, num_added_dds * sizeof(z_t));
  memcpy(x0.counter, lower_bounds, num_added_dds * sizeof(z_t));

  x0.num_added_dds = num_added_dds;
  x0.num_old_dds = CLAMP_DIFF(amb_test->sats.num_sats, 1);

  /* Construct the mapping from the old prn indices into the new,
   * and from the added prn indices into the new. */
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
      log_error("add_sats_old: impossible condition reached.\n");
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

  log_info("IAR: %"PRIu32" hypotheses before inclusion\n", memory_pool_n_allocated(amb_test->pool));
  if (DEBUG) {
    memory_pool_map(amb_test->pool, &x0.num_old_dds, &print_hyp);
  }
  memcpy(x0.Z_inv, Z_inv, num_added_dds * num_added_dds * sizeof(z_t));
  /* Take the product of our current hypothesis state with the generator, recorrelating the new ones as we go. */
  memory_pool_product_generator(amb_test->pool, &x0, MAX_HYPOTHESES, sizeof(x0),
                                &no_init, &generate_next_hypothesis, &hypothesis_prod);
  log_info("IAR: updates to %"PRIu32"\n", memory_pool_n_allocated(amb_test->pool));
  if (DEBUG) {
    memory_pool_map(amb_test->pool, &k, &print_hyp);
  }
}

void init_residual_matrices(residual_mtxs_t *res_mtxs, u8 num_dds, double *DE_mtx, double *obs_cov)
{
  res_mtxs->res_dim = num_dds + CLAMP_DIFF(num_dds, 3);
  res_mtxs->null_space_dim = CLAMP_DIFF(num_dds, 3);
  assign_phase_obs_null_basis(num_dds, DE_mtx, res_mtxs->null_projector);
  assign_residual_covariance_inverse(num_dds, obs_cov, res_mtxs->null_projector, res_mtxs->half_res_cov_inv);
}

void assign_residual_covariance_inverse(u8 num_dds, double *obs_cov, double *q, double *r_cov_inv) //TODO make this more efficient (e.g. via page 3/6.2-3/2014 of ian's notebook)
{
  integer dd_dim = 2*num_dds;
  integer res_dim = num_dds + CLAMP_DIFF(num_dds, 3);
  u32 nullspace_dim = CLAMP_DIFF(num_dds, 3);
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
    r_vec[i + res_mtxs->null_space_dim] =
      simple_amb_measurement(dd_measurements[i],
                             dd_measurements[i+num_dds]);
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
