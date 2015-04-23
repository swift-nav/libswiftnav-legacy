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

#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "logging.h"
#include "amb_kf.h"
#include "baseline.h"
#include "single_diff.h"
#include "dgnss_management.h"
#include "linear_algebra.h"
#include "ambiguity_test.h"

nkf_t nkf;
sats_management_t sats_management;
ambiguity_test_t ambiguity_test;

dgnss_settings_t dgnss_settings = {
  .phase_var_test = DEFAULT_PHASE_VAR_TEST,
  .code_var_test = DEFAULT_CODE_VAR_TEST,
  .phase_var_kf = DEFAULT_PHASE_VAR_KF,
  .code_var_kf = DEFAULT_CODE_VAR_KF,
  .amb_drift_var = DEFAULT_AMB_DRIFT_VAR,
  .amb_init_var = DEFAULT_AMB_INIT_VAR,
  .new_int_var = DEFAULT_NEW_INT_VAR,
};

void dgnss_set_settings(double phase_var_test, double code_var_test,
                        double phase_var_kf, double code_var_kf,
                        double amb_drift_var, double amb_init_var,
                        double new_int_var)
{
  dgnss_settings.phase_var_test = phase_var_test;
  dgnss_settings.code_var_test  = code_var_test;
  dgnss_settings.phase_var_kf   = phase_var_kf;
  dgnss_settings.code_var_kf    = code_var_kf;
  dgnss_settings.amb_drift_var  = amb_drift_var;
  dgnss_settings.amb_init_var   = amb_init_var;
  dgnss_settings.new_int_var    = new_int_var;
}

void make_measurements(u8 num_double_diffs, const sdiff_t *sdiffs, double *raw_measurements)
{
  DEBUG_ENTRY();

  double phase0 = sdiffs[0].carrier_phase;
  double code0 = sdiffs[0].pseudorange;
  for (u8 i=0; i<num_double_diffs; i++) {
    raw_measurements[i] = sdiffs[i+1].carrier_phase - phase0;
    raw_measurements[i+num_double_diffs] = sdiffs[i+1].pseudorange - code0;
  }

  DEBUG_EXIT();
}

static bool prns_match(const u8 *old_non_ref_prns, u8 num_non_ref_sdiffs,
                       const sdiff_t *non_ref_sdiffs)
{
  if (sats_management.num_sats-1 != num_non_ref_sdiffs) {
    /* lengths don't match */
    return false;
  }
  for (u8 i=0; i<num_non_ref_sdiffs; i++) {
    /* iterate through the non-reference_sats, checking they match. */
    if (old_non_ref_prns[i] != non_ref_sdiffs[i].prn) {
      return false;
    }
  }
  return true;
}

/** Finds the prns of the intersection between old prns and new measurements.
 * It returns the length of the intersection
 */
static u8 dgnss_intersect_sats(u8 num_old_prns, const u8 *old_prns,
                               u8 num_sdiffs, const sdiff_t *sdiffs,
                               u8 *ndx_of_intersection_in_old,
                               u8 *ndx_of_intersection_in_new)
{
  u8 i, j, n = 0;
  /* Loop over old_prns and sdiffs and check if a PRN is present in both. */
  for (i=0, j=0; i<num_old_prns && j<num_sdiffs; i++, j++) {
    if (old_prns[i] < sdiffs[j].prn)
      j--;
    else if (old_prns[i] > sdiffs[j].prn)
      i--;
    else {
      ndx_of_intersection_in_old[n] = i;
      ndx_of_intersection_in_new[n] = j;
      n++;
    }
  }
  return n;
}

void dgnss_init(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3])
{
  DEBUG_ENTRY();

  sdiff_t corrected_sdiffs[num_sats];
  init_sats_management(&sats_management, num_sats, sdiffs, corrected_sdiffs);

  create_ambiguity_test(&ambiguity_test);

  if (num_sats <= 1) {
    DEBUG_EXIT();
    return;
  }

  double dd_measurements[2*(num_sats-1)];
  make_measurements(num_sats-1, corrected_sdiffs, dd_measurements);

  set_nkf(
    &nkf,
    dgnss_settings.amb_drift_var,
    dgnss_settings.phase_var_kf, dgnss_settings.code_var_kf,
    dgnss_settings.amb_init_var,
    num_sats, corrected_sdiffs, dd_measurements, reciever_ecef
  );

  DEBUG_EXIT();
}

void dgnss_rebase_ref(u8 num_sdiffs, sdiff_t *sdiffs, double reciever_ecef[3], u8 old_prns[MAX_CHANNELS], sdiff_t *corrected_sdiffs)
{
  (void)reciever_ecef;
  /* all the ref sat stuff */
  s8 sats_management_code = rebase_sats_management(&sats_management, num_sdiffs, sdiffs, corrected_sdiffs);
  if (sats_management_code == NEW_REF_START_OVER) {
    log_info("Unable to rebase to new ref, resetting filters and starting over\n");
    dgnss_init(num_sdiffs, sdiffs, reciever_ecef);
    memcpy(old_prns, sats_management.prns, sats_management.num_sats * sizeof(u8));
    if (num_sdiffs >= 1) {
      copy_sdiffs_put_ref_first(old_prns[0], num_sdiffs, sdiffs, corrected_sdiffs);
    }
    /*dgnss_init(num_sdiffs, sdiffs, reciever_ecef); //TODO use current baseline state*/
    return;
  }
  else if (sats_management_code == NEW_REF) {
    /* do everything related to changing the reference sat here */
    rebase_nkf(&nkf, sats_management.num_sats, &old_prns[0], &sats_management.prns[0]);
  }
}


static void sdiffs_to_prns(u8 n, sdiff_t *sdiffs, u8 *prns)
{
  for (u8 i=0; i<n; i++) {
    prns[i] = sdiffs[i].prn;
  }
}

/** Single timestep measurement of the ambiguity vector given sdiffs.
 * Using just the scalar DD carrier phase and pseudoranges for each channel,
 * estimate the integer ambiguities.
 * See docs for simple_amb_measurement.
 *
 * \param num_sdiffs            The number of sdiffs including the reference.
 * \param sdiffs_with_ref_first The sdiffs sorted by prn, but with the ref first.
 */
static void dgnss_simple_amb_meas(const u8 num_sdiffs,
                                  const sdiff_t *sdiffs_with_ref_first,
                                  double *est)
{
  double ref_phase = sdiffs_with_ref_first[0].carrier_phase;
  double ref_code  = sdiffs_with_ref_first[0].pseudorange;
  for (u8 i = 0; i < num_sdiffs - 1; i++) {
    double phase = sdiffs_with_ref_first[i+1].carrier_phase;
    double code  = sdiffs_with_ref_first[i+1].pseudorange;
    phase -= ref_phase;
    code  -= ref_code;
    est[i] = simple_amb_measurement(phase, code);
  }
}

static void dgnss_update_sats(u8 num_sdiffs, double reciever_ecef[3],
                              sdiff_t *sdiffs_with_ref_first,
                              double *dd_measurements)
{
  DEBUG_ENTRY();

  (void)dd_measurements;
  u8 new_prns[num_sdiffs];
  sdiffs_to_prns(num_sdiffs, sdiffs_with_ref_first, new_prns);

  u8 old_prns[MAX_CHANNELS];
  memcpy(old_prns, sats_management.prns, sats_management.num_sats * sizeof(u8));

  if (!prns_match(&old_prns[1], num_sdiffs-1, &sdiffs_with_ref_first[1])) {
    u8 ndx_of_intersection_in_old[sats_management.num_sats];
    u8 ndx_of_intersection_in_new[sats_management.num_sats];
    ndx_of_intersection_in_old[0] = 0;
    ndx_of_intersection_in_new[0] = 0;
    u8 num_intersection_sats = dgnss_intersect_sats(
        sats_management.num_sats-1, &old_prns[1],
        num_sdiffs-1, &sdiffs_with_ref_first[1],
        &ndx_of_intersection_in_old[1],
        &ndx_of_intersection_in_new[1]) + 1;

    set_nkf_matrices(
      &nkf,
      dgnss_settings.phase_var_kf, dgnss_settings.code_var_kf,
      num_sdiffs, sdiffs_with_ref_first, reciever_ecef
    );

    if (num_intersection_sats < sats_management.num_sats) { /* we lost sats */
      nkf_state_projection(&nkf,
                           sats_management.num_sats-1,
                           num_intersection_sats-1,
                           &ndx_of_intersection_in_old[1]);
    }
    if (num_intersection_sats < num_sdiffs) { /* we gained sats */
      double simple_estimates[num_sdiffs-1];
      dgnss_simple_amb_meas(num_sdiffs, sdiffs_with_ref_first,
                            simple_estimates);
      nkf_state_inclusion(&nkf,
                          num_intersection_sats-1,
                          num_sdiffs-1,
                          &ndx_of_intersection_in_new[1],
                          simple_estimates,
                          dgnss_settings.new_int_var);
    }

    update_sats_sats_management(&sats_management, num_sdiffs-1, &sdiffs_with_ref_first[1]);
  }
  else {
    set_nkf_matrices(
      &nkf,
      dgnss_settings.phase_var_kf, dgnss_settings.code_var_kf,
      num_sdiffs, sdiffs_with_ref_first, reciever_ecef
    );
  }

  DEBUG_EXIT();
}

static void dgnss_incorporate_observation(sdiff_t *sdiffs, double * dd_measurements,
                                          double *reciever_ecef)
{
  DEBUG_ENTRY();

  double b2[3];
  least_squares_solve_b(&nkf, sdiffs, dd_measurements, reciever_ecef, b2);

  double ref_ecef[3];

  ref_ecef[0] = reciever_ecef[0] + 0.5 * b2[0];
  ref_ecef[1] = reciever_ecef[1] + 0.5 * b2[0];
  ref_ecef[2] = reciever_ecef[2] + 0.5 * b2[0];

  /* TODO: make a common DE and use it instead. */

  set_nkf_matrices(&nkf,
                   dgnss_settings.phase_var_kf, dgnss_settings.code_var_kf,
                   sats_management.num_sats, sdiffs, ref_ecef);

  nkf_update(&nkf, dd_measurements);
  DEBUG_EXIT();
}

void dgnss_update(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3])
{
  DEBUG_ENTRY();
  if (DEBUG) {
    printf("sdiff[*].prn = {");
    for (u8 i=0; i < num_sats; i++) {
      printf("%u, ",sdiffs[i].prn);
    }
    printf("}\n");
  }

  if (num_sats <= 1) {
    sats_management.num_sats = num_sats;
    if (num_sats == 1) {
      sats_management.prns[0] = sdiffs[0].prn;
    }
    create_ambiguity_test(&ambiguity_test);
    DEBUG_EXIT();
    return;
  }

  if (sats_management.num_sats <= 1) {
    dgnss_init(num_sats, sdiffs, reciever_ecef);
  }

  sdiff_t sdiffs_with_ref_first[num_sats];

  u8 old_prns[MAX_CHANNELS];
  memcpy(old_prns, sats_management.prns, sats_management.num_sats * sizeof(u8));

  /* rebase globals to a new reference sat
   * (permutes sdiffs_with_ref_first accordingly) */
  dgnss_rebase_ref(num_sats, sdiffs, reciever_ecef, old_prns, sdiffs_with_ref_first);

  double dd_measurements[2*(num_sats-1)];
  make_measurements(num_sats-1, sdiffs_with_ref_first, dd_measurements);

  /* all the added/dropped sat stuff */
  dgnss_update_sats(num_sats, reciever_ecef, sdiffs_with_ref_first, dd_measurements);

  double ref_ecef[3];
  if (num_sats >= 5) {
    dgnss_incorporate_observation(sdiffs_with_ref_first, dd_measurements, reciever_ecef);

    double b2[3];
    least_squares_solve_b(&nkf, sdiffs_with_ref_first, dd_measurements, reciever_ecef, b2);

    ref_ecef[0] = reciever_ecef[0] + 0.5 * b2[0];
    ref_ecef[1] = reciever_ecef[1] + 0.5 * b2[1];
    ref_ecef[2] = reciever_ecef[2] + 0.5 * b2[2];
  }

  u8 changed_sats = ambiguity_update_sats(&ambiguity_test, num_sats, sdiffs,
                                          &sats_management, nkf.state_mean,
                                          nkf.state_cov_U, nkf.state_cov_D);

  update_ambiguity_test(ref_ecef,
                        dgnss_settings.phase_var_test,
                        dgnss_settings.code_var_test,
                        &ambiguity_test, nkf.state_dim,
                        sdiffs, changed_sats);

  update_unanimous_ambiguities(&ambiguity_test);

  if (DEBUG) {
    if (num_sats >=4) {
      double b3[3];
      least_squares_solve_b(&nkf, sdiffs_with_ref_first, dd_measurements, reciever_ecef, b3);

      ref_ecef[0] = reciever_ecef[0] + 0.5 * b3[0];
      ref_ecef[1] = reciever_ecef[1] + 0.5 * b3[1];
      ref_ecef[2] = reciever_ecef[2] + 0.5 * b3[2];
      double bb[3];
      u8 num_used;
      dgnss_fixed_baseline(num_sats, sdiffs, ref_ecef,
                           &num_used, bb);
      log_debug("\ndgnss_fixed_baseline:\nb = %f, \t%f, \t%f\nnum_used/num_sats = %u/%u\nusing_iar = %u\n\n",
             bb[0], bb[1], bb[2],
             num_used, num_sats,
             ambiguity_iar_can_solve(&ambiguity_test));
    }
  }
  DEBUG_EXIT();
}

u32 dgnss_iar_num_hyps(void)
{
  if (ambiguity_test.pool == NULL) {
    return 0;
  } else {
    return ambiguity_test_n_hypotheses(&ambiguity_test);
  }
}

u32 dgnss_iar_num_sats(void)
{
  return ambiguity_test.sats.num_sats;
}

s8 dgnss_iar_get_single_hyp(double *dhyp)
{
  u8 num_dds = ambiguity_test.sats.num_sats;
  s32 hyp[num_dds];
  s8 ret = get_single_hypothesis(&ambiguity_test, hyp);
  for (u8 i=0; i<num_dds; i++) {
    dhyp[i] = hyp[i];
  }
  return ret;
}

void dgnss_new_float_baseline(u8 num_sats, sdiff_t *sdiffs, double receiver_ecef[3], u8 *num_used, double b[3])
{
  DEBUG_ENTRY();
  sdiff_t corrected_sdiffs[num_sats];

  u8 old_prns[MAX_CHANNELS];
  memcpy(old_prns, sats_management.prns, sats_management.num_sats * sizeof(u8));
  /* Rebase globals to a new reference sat
   * (permutes corrected_sdiffs accordingly) */
  dgnss_rebase_ref(num_sats, sdiffs, receiver_ecef, old_prns, corrected_sdiffs);

  double dd_measurements[2*(num_sats-1)];
  make_measurements(num_sats-1, corrected_sdiffs, dd_measurements);

  least_squares_solve_b(&nkf, corrected_sdiffs, dd_measurements, receiver_ecef, b);
  *num_used = sats_management.num_sats;
  DEBUG_EXIT();
}

/* Returns the fixed baseline iff there are at least 3 dd ambs unanimously
 * agreed upon in the ambiguity_test.
 * \return 1 If fixed baseline calculation succeeds
 *         0 If iar cannot solve or an error occurs. Signals that float
 *           baseline is needed instead.
 */
s8 dgnss_fixed_baseline(u8 num_sdiffs, sdiff_t *sdiffs, double ref_ecef[3],
                        u8 *num_used, double b[3])
{
  if (!ambiguity_iar_can_solve(&ambiguity_test)) {
    return 0;
  }

  sdiff_t ambiguity_sdiffs[ambiguity_test.amb_check.num_matching_ndxs+1];
  double dd_meas[2 * ambiguity_test.amb_check.num_matching_ndxs];

  s8 valid_sdiffs = make_ambiguity_resolved_dd_measurements_and_sdiffs(
      &ambiguity_test, num_sdiffs, sdiffs, dd_meas, ambiguity_sdiffs);

  /* At this point, sdiffs should be valid due to dgnss_update
   * Return code not equal to 0 signals an error. */
  if (valid_sdiffs != 0) {
    if (valid_sdiffs != -1) {
      log_error("dgnss_fixed_baseline: Invalid sdiffs.");
    }
    return 0;
  }

  double DE[ambiguity_test.amb_check.num_matching_ndxs * 3];
  assign_de_mtx(ambiguity_test.amb_check.num_matching_ndxs + 1,
                ambiguity_sdiffs, ref_ecef, DE);
  *num_used = ambiguity_test.amb_check.num_matching_ndxs + 1;
  s8 ret = lesq_solution_int(ambiguity_test.amb_check.num_matching_ndxs, dd_meas,
                             ambiguity_test.amb_check.ambs, DE, b, 0);
  if (ret) {
    log_error("dgnss_fixed_baseline: "
              "lesq_solution returned error %d\n", ret);
    DEBUG_EXIT();
    return 0;
  }
  return 1;
}

/** Constructs a low latency float baseline measurement.
 * The sdiffs have no particular reason (other than a general tendency
 * brought by hysteresis) to match up with the float filter's sats, so we have
 * to check if we can solve. For now, unless the sdiffs are a superset of the
 * float sats, we don't solve.
 *
 * Requires num_sdiffs >= 4 and (global) sats_management.num_sats >= 4.
 *
 * \TODO solve whenever the information is there.
 *
 * \TODO since we're now using make_dd_measurements_and_sdiffs outside of the
 * amb_test context, pull it into another file.
 *
 * \TODO pull this function into the KF, once we pull the sats_management struct
 *      into the KF too. When we do, do the same for the IAR low lat solution.
 *
 * \param num_sdiffs  The number of sdiffs input. Must be >= 4.
 * \param sdiffs      The sdiffs used to measure. (These should be a superset
 *                    of the float sats).
 * \param ref_ecef    The reference position used for solving, and making
 *                    observation matrices.
 * \param num_used    The number of sats actually used to compute the baseline.
 * \param b           The baseline computed.
 * \return -1 if it can't solve.
 *          0 If it can solve.
 */
s8 _dgnss_low_latency_float_baseline(u8 num_sdiffs, sdiff_t *sdiffs,
                                     double ref_ecef[3], u8 *num_used, double b[3])
{
  DEBUG_ENTRY();
  if (num_sdiffs < 4 || sats_management.num_sats < 4) {
    /* For a position solution, we need at least 4 sats. That means we must
     * have at least 4 sats in common between what the KF is tracking and
     * the sdiffs we give this function. If either is less than 4,
     * this criterion cannot be satisfied. */
    log_debug("Low latency solution can't be computed. Too few observations"
              " or too few sats in the current filter.\n");
    DEBUG_EXIT();
    return -1;
  }

  double float_dd_measurements[2 * (sats_management.num_sats - 1)];
  sdiff_t float_sdiffs[sats_management.num_sats];
  s8 can_make_obs = make_dd_measurements_and_sdiffs(sats_management.prns[0],
             &sats_management.prns[1], sats_management.num_sats - 1,
             num_sdiffs, sdiffs,
             float_dd_measurements, float_sdiffs);
  if (can_make_obs == -1) {
    log_debug("make_float_dd_measurements_and_sdiffs has error code -1\n");
    DEBUG_EXIT();
    return -1;
  }
  least_squares_solve_b(&nkf, float_sdiffs, float_dd_measurements,
                        ref_ecef, b);
  *num_used = sats_management.num_sats;
  DEBUG_EXIT();
  return 0;
}

/** Constructs a low latency IAR resolved baseline measurement.
 * The sdiffs have no particular reason (other than a general tendency
 * brought by hysteresis) to match up with the IAR sats, so we have
 * to check if we can solve. For now, unless the sdiffs are a superset of the
 * IAR sats, we don't solve.
 *
 * Requires num_sdiffs >= 4.
 *
 * \TODO, solve whenever we can
 *
 * \TODO since we're now using make_dd_measurements_and_sdiffs outside of the
 * amb_test context, pull it into another file.
 *
 * \TODO pull this into the IAR file when we do the same for the float low lat
 *      solution.
 *
 * \todo This function is identical to dgnss_fixed_baseline()?
 *
 * \param num_sdiffs  The number of sdiffs input.
 * \param sdiffs      The sdiffs used to measure. (These should be a superset
 *                    of the float sats).
 * \param ref_ecef    The reference position used for solving, and making
 *                    observation matrices.
 * \param num_used    The number of sats actually used to compute the baseline.
 * \param b           The baseline computed.
 * \return -1 if it can't solve.
 *          0 If it can solve.
 */
s8 _dgnss_low_latency_IAR_baseline(u8 num_sdiffs, sdiff_t *sdiffs,
                                   double ref_ecef[3], u8 *num_used, double b[3])
{
  DEBUG_ENTRY();
  assert(num_sdiffs >= 4);
  if (!ambiguity_iar_can_solve(&ambiguity_test)) {
    DEBUG_EXIT();
    return -1;
  }

  sdiff_t ambiguity_sdiffs[ambiguity_test.amb_check.num_matching_ndxs+1];
  double dd_meas[2 * ambiguity_test.amb_check.num_matching_ndxs];

  s8 valid_sdiffs = make_ambiguity_resolved_dd_measurements_and_sdiffs(
      &ambiguity_test, num_sdiffs, sdiffs, dd_meas, ambiguity_sdiffs);

  if (valid_sdiffs != 0) {
    if (valid_sdiffs != -1) {
      log_error("_dgnss_low_latency_IAR_baseline: Invalid sdiffs.");
    }
    DEBUG_EXIT();
    return -1;
  }

  /* TODO: check internals of this if's content and abstract it from the KF */
  double DE[ambiguity_test.amb_check.num_matching_ndxs * 3];
  assign_de_mtx(ambiguity_test.amb_check.num_matching_ndxs + 1,
                ambiguity_sdiffs, ref_ecef, DE);
  *num_used = ambiguity_test.amb_check.num_matching_ndxs + 1;
  s8 ret = lesq_solution_int(ambiguity_test.amb_check.num_matching_ndxs,
                             dd_meas, ambiguity_test.amb_check.ambs, DE, b, 0);
  if (ret) {
    log_error("_dgnss_low_latency_IAR_baseline: "
              "lesq_solution returned error %d\n", ret);
    DEBUG_EXIT();
    return -1;
  }

  DEBUG_EXIT();
  return 0;
}

/** Finds the baseline using low latency sdiffs.
 * The low latency sdiffs are not guaranteed to match up with either the
 * amb_test's or the float sdiffs, and thus care must be taken to transform them
 * accordingly and check their availability in those sat sets.
 *
 * \param num_sdiffs  The number of low-latency sdiffs provided.
 * \param sdiffs      The low-latency sdiffs.
 * \param ref_ecef    The referece position for the baseline.
 *                    (TODO is this the local or remote receiver position?)
 * \param num_used    Output number of sdiffs actually used in the baseline
 *                    estimate.
 * \param b           Output baseline.
 * \return  1 if we are using an IAR resolved baseline.
 *          2 if we are using a float baseline.
 *         -1 if we can't give a baseline.
 */
s8 dgnss_low_latency_baseline(u8 num_sdiffs, sdiff_t *sdiffs,
                              double ref_ecef[3], u8 *num_used, double b[3])
{
  DEBUG_ENTRY();
  if (num_sdiffs < 4 || sats_management.num_sats < 4) {
    /* For a position solution, we need at least 4 sats. That means we must
     * have at least 4 sats in common between what the filters are tracking and
     * the sdiffs we give this function. If num_sdiffs, or the number of KF
     * sats is less than 4, this criterion cannot be satisfied. */
    log_debug("Low latency solution can't be computed. Too few observations"
              " or too few sats in the current filter.\n");
    DEBUG_EXIT();
    return -1;
  }
  if (0 == _dgnss_low_latency_IAR_baseline(num_sdiffs, sdiffs,
                                  ref_ecef, num_used, b)) {
    log_debug("low latency IAR solution\n");
    DEBUG_EXIT();
    return 1;
  }
  /* if we get here, we weren't able to get an IAR resolved baseline.
   * Check if we can get a float baseline. */
  s8 float_ret_code = _dgnss_low_latency_float_baseline(num_sdiffs, sdiffs,
                                              ref_ecef, num_used, b);
  if (float_ret_code == 0) {
    log_debug("low latency float solution\n");
    DEBUG_EXIT();
    return 2;
  }
  log_debug("no low latency solution\n");
  DEBUG_EXIT();
  return -1;
}

void dgnss_reset_iar()
{
  create_ambiguity_test(&ambiguity_test);
}

void dgnss_init_known_baseline(u8 num_sats, sdiff_t *sdiffs,
                               double receiver_ecef[3], double b[3])
{
  double ref_ecef[3];
  ref_ecef[0] = receiver_ecef[0] + 0.5 * b[0];
  ref_ecef[1] = receiver_ecef[1] + 0.5 * b[1];
  ref_ecef[2] = receiver_ecef[2] + 0.5 * b[2];

  sdiff_t corrected_sdiffs[num_sats];

  u8 old_prns[MAX_CHANNELS];
  memcpy(old_prns, sats_management.prns, sats_management.num_sats * sizeof(u8));
  /* rebase globals to a new reference sat
   * (permutes corrected_sdiffs accordingly) */
  dgnss_rebase_ref(num_sats, sdiffs, ref_ecef, old_prns, corrected_sdiffs);

  double dds[2*(num_sats-1)];
  make_measurements(num_sats-1, corrected_sdiffs, dds);

  double DE[(num_sats-1)*3];
  assign_de_mtx(num_sats, corrected_sdiffs, ref_ecef, DE);

  dgnss_reset_iar();

  memcpy(&ambiguity_test.sats, &sats_management, sizeof(sats_management));
  hypothesis_t *hyp = (hypothesis_t *)memory_pool_add(ambiguity_test.pool);
  hyp->ll = 0;
  amb_from_baseline(num_sats-1, DE, dds, b, hyp->N);

  double obs_cov[(num_sats-1) * (num_sats-1) * 4];
  memset(obs_cov, 0, (num_sats-1) * (num_sats-1) * 4 * sizeof(double));
  u8 num_dds = num_sats-1;
  for (u8 i=0; i<num_dds; i++) {
    for (u8 j=0; j<num_dds; j++) {
      u8 i_ = i+num_dds;
      u8 j_ = j+num_dds;
      if (i==j) {
        obs_cov[i*2*num_dds + j] = dgnss_settings.phase_var_test * 2;
        obs_cov[i_*2*num_dds + j_] = dgnss_settings.code_var_test * 2;
      }
      else {
        obs_cov[i*2*num_dds + j] = dgnss_settings.phase_var_test;
        obs_cov[i_*2*num_dds + j_] = dgnss_settings.code_var_test;
      }
    }
  }

  init_residual_matrices(&ambiguity_test.res_mtxs, num_sats-1, DE, obs_cov);
}

static void measure_b(u8 state_dim, const double *state_mean,
                u8 num_sdiffs, sdiff_t *sdiffs_with_ref_first,
                const double receiver_ecef[3], double *b)
{
  double dd_measurements[2*(num_sdiffs-1)];
  make_measurements(num_sdiffs - 1, sdiffs_with_ref_first, dd_measurements);
  double b_old[3] = {0, 0, 0};
  double ref_ecef[3];
  ref_ecef[0] = receiver_ecef[0];
  ref_ecef[1] = receiver_ecef[1];
  ref_ecef[2] = receiver_ecef[2];

  least_squares_solve_b_external_ambs(state_dim, state_mean,
      sdiffs_with_ref_first, dd_measurements, ref_ecef, b);

  while (vector_distance(3, b_old, b) > 1e-4) {
    memcpy(b_old, b, sizeof(double)*3);
    ref_ecef[0] = receiver_ecef[0] + 0.5 * b_old[0];
    ref_ecef[1] = receiver_ecef[1] + 0.5 * b_old[1];
    ref_ecef[2] = receiver_ecef[2] + 0.5 * b_old[2];
    least_squares_solve_b_external_ambs(state_dim, state_mean,
        sdiffs_with_ref_first, dd_measurements, ref_ecef, b);
  }
}


void measure_b_with_external_ambs(u8 state_dim, const double *state_mean,
                                  u8 num_sdiffs, sdiff_t *sdiffs,
                                  const double receiver_ecef[3], double *b)
{
  DEBUG_ENTRY();

  sdiff_t sdiffs_with_ref_first[num_sdiffs];
  /* We require the sats updating has already been done with these sdiffs */
  u8 ref_prn = sats_management.prns[0];
  copy_sdiffs_put_ref_first(ref_prn, num_sdiffs, sdiffs, sdiffs_with_ref_first);

  measure_b(state_dim, state_mean, num_sdiffs, sdiffs_with_ref_first, receiver_ecef, b);

  DEBUG_EXIT();
}

void measure_amb_kf_b(u8 num_sdiffs, sdiff_t *sdiffs,
                      const double receiver_ecef[3], double *b)
{
  DEBUG_ENTRY();

  sdiff_t sdiffs_with_ref_first[num_sdiffs];
  /* We require the sats updating has already been done with these sdiffs */
  u8 ref_prn = sats_management.prns[0];
  copy_sdiffs_put_ref_first(ref_prn, num_sdiffs, sdiffs, sdiffs_with_ref_first);

  measure_b( nkf.state_dim, nkf.state_mean,
      num_sdiffs, sdiffs_with_ref_first, receiver_ecef, b);

  DEBUG_EXIT();
}

void measure_iar_b_with_external_ambs(double *state_mean,
                                      u8 num_sdiffs, sdiff_t *sdiffs,
                                      double receiver_ecef[3],
                                      double *b)
{
  DEBUG_ENTRY();

  sdiff_t sdiffs_with_ref_first[num_sdiffs];
  match_sdiffs_to_sats_man(&ambiguity_test.sats, num_sdiffs, sdiffs, sdiffs_with_ref_first);

  measure_b(CLAMP_DIFF(ambiguity_test.sats.num_sats, 1), state_mean,
      num_sdiffs, sdiffs_with_ref_first, receiver_ecef, b);

  DEBUG_EXIT();
}

static u8 get_de_and_phase(sats_management_t *sats_man,
                           u8 num_sdiffs, sdiff_t *sdiffs,
                           double ref_ecef[3],
                           double *de, double *phase)
{
  u8 ref_prn = sats_man->prns[0];
  u8 num_sats = sats_man->num_sats;
  double e0[3];
  double phi0 = 0;
  /* TODO: Detect if ref_prn is not in prns and return error? */
  u8 i;
  for (i=0; i<num_sdiffs; i++) {
    if (sdiffs[i].prn == ref_prn) {
      e0[0] = sdiffs[i].sat_pos[0] - ref_ecef[0];
      e0[1] = sdiffs[i].sat_pos[1] - ref_ecef[1];
      e0[2] = sdiffs[i].sat_pos[2] - ref_ecef[2];
      vector_normalize(3, e0);
      phi0 = sdiffs[i].carrier_phase;
      break;
    }
  }
  i=1;
  u8 j = 0;
  while (i < num_sats) {
    if (sdiffs[j].prn < sats_man->prns[i]) {
      j++;
    }
    else if (sdiffs[j].prn > sats_man->prns[i]) {
      /* This should never happen. */
      log_warn("sdiffs should be a super set of sats_man prns\n");
      i++;
    }
    else {  /* else they match */
      double e[3];
      e[0] = sdiffs[j].sat_pos[0] - ref_ecef[0];
      e[1] = sdiffs[j].sat_pos[1] - ref_ecef[1];
      e[2] = sdiffs[j].sat_pos[2] - ref_ecef[2];
      vector_normalize(3, e);
      de[(i-1)*3    ] = e[0] - e0[0];
      de[(i-1)*3 + 1] = e[1] - e0[1];
      de[(i-1)*3 + 2] = e[2] - e0[2];
      phase[i-1] = sdiffs[j].carrier_phase - phi0;
      i++;
      j++;
    }
  }
  return num_sats;
}

u8 get_amb_kf_de_and_phase(u8 num_sdiffs, sdiff_t *sdiffs,
                           double ref_ecef[3],
                           double *de, double *phase)
{
  return get_de_and_phase(&sats_management,
                          num_sdiffs, sdiffs,
                          ref_ecef,
                          de, phase);
}

u8 get_iar_de_and_phase(u8 num_sdiffs, sdiff_t *sdiffs,
                        double ref_ecef[3],
                        double *de, double *phase)
{
  return get_de_and_phase(&ambiguity_test.sats,
                          num_sdiffs, sdiffs,
                          ref_ecef,
                          de, phase);
}

u8 get_amb_kf_mean(double *ambs)
{
  u8 num_dds = CLAMP_DIFF(sats_management.num_sats, 1);
  memcpy(ambs, nkf.state_mean, num_dds * sizeof(double));
  return num_dds;
}

u8 get_amb_kf_cov(double *cov)
{
  u8 num_dds = CLAMP_DIFF(sats_management.num_sats, 1);
  matrix_reconstruct_udu(num_dds, nkf.state_cov_U, nkf.state_cov_D, cov);
  return num_dds;
}

u8 get_amb_kf_prns(u8 *prns)
{
  memcpy(prns, sats_management.prns, sats_management.num_sats * sizeof(u8));
  return sats_management.num_sats;
}

u8 get_amb_test_prns(u8 *prns)
{
  memcpy(prns, ambiguity_test.sats.prns, ambiguity_test.sats.num_sats * sizeof(u8));
  return ambiguity_test.sats.num_sats;
}

s8 dgnss_iar_resolved()
{
  return ambiguity_iar_can_solve(&ambiguity_test);
}

u8 dgnss_iar_pool_contains(double *ambs)
{
  return ambiguity_test_pool_contains(&ambiguity_test, ambs);
}

double dgnss_iar_pool_ll(u8 num_ambs, double *ambs)
{
  return ambiguity_test_pool_ll(&ambiguity_test, num_ambs, ambs);
}

double dgnss_iar_pool_prob(u8 num_ambs, double *ambs)
{
  return ambiguity_test_pool_prob(&ambiguity_test, num_ambs, ambs);
}

u8 dgnss_iar_MLE_ambs(s32 *ambs)
{
  ambiguity_test_MLE_ambs(&ambiguity_test, ambs);
  return CLAMP_DIFF(ambiguity_test.sats.num_sats, 1);
}

nkf_t* get_dgnss_nkf(void)
{
  return &nkf;
}

sats_management_t* get_sats_management(void)
{
  return &sats_management;
}

ambiguity_test_t* get_ambiguity_test(void)
{
  return &ambiguity_test;
}
