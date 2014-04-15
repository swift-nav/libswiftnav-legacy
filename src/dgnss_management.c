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
#include "float_kf.h"
#include "stupid_filter.h"
#include "single_diff.h"
#include "dgnss_management.h"
#include "linear_algebra.h"
#include "ambiguity_test.h"

kf_t kf;
stupid_filter_state_t stupid_state;
sats_management_t sats_management;
ambiguity_test_t ambiguity_test;

void make_measurements(u8 num_double_diffs, sdiff_t *sdiffs, double *raw_measurements)
{
  double phase0 = sdiffs[0].carrier_phase;
  double code0 = sdiffs[0].pseudorange;
  for (u8 i=0; i<num_double_diffs; i++) {
    raw_measurements[i] = sdiffs[i+1].carrier_phase - phase0;
    raw_measurements[i+num_double_diffs] = sdiffs[i+1].pseudorange - code0;
  }
}

bool prns_match(u8 *old_non_ref_prns, u8 num_non_ref_sdiffs, sdiff_t *non_ref_sdiffs)
{
  if (sats_management.num_sats-1 != num_non_ref_sdiffs) {
    // printf("len_doesn't match\n");
    return false;
  }
  for (u8 i=0; i<num_non_ref_sdiffs; i++) {
    //iterate through the non-reference_sats
    /*printf("old[%u]=%u, new[%u]=%u\n", i+1, old_non_ref_prns[i], i+1, non_ref_sdiffs[i].prn);*/
    if (old_non_ref_prns[i] != non_ref_sdiffs[i].prn) {
      return false;
    }
  }
  /*printf("prns_match\n");*/
  return true;
}

/** Finds the prns of the intersection between old prns and new measurements.
 * It returns the length of the intersection
 */
u8 dgnss_intersect_sats(u8 num_old_prns, u8 *old_prns,
                  u8 num_sdiffs, sdiff_t *sdiffs,
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

void dgnss_init(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3], double dt)
{
  sdiff_t corrected_sdiffs[num_sats];
  init_sats_management(&sats_management, num_sats, sdiffs, corrected_sdiffs);

  double dd_measurements[2*(num_sats-1)];
  make_measurements(num_sats-1, corrected_sdiffs, dd_measurements);

  get_kf(&kf, PHASE_VAR, CODE_VAR,
         POS_TRANS_VAR, VEL_TRANS_VAR, INT_TRANS_VAR,
         POS_INIT_VAR,  VEL_INIT_VAR,  INT_INIT_VAR,
         num_sats, corrected_sdiffs, dd_measurements, reciever_ecef, dt);

  create_ambiguity_test(&ambiguity_test);
}

void dgnss_rebase_ref(u8 num_sdiffs, sdiff_t *sdiffs, double reciever_ecef[3], double dt, u8 old_prns[MAX_CHANNELS], sdiff_t *corrected_sdiffs)
{
  (void)dt; (void)reciever_ecef;
  //all the ref sat stuff
  s8 sats_management_code = rebase_sats_management(&sats_management, num_sdiffs, sdiffs, corrected_sdiffs);
  if (sats_management_code == NEW_REF_START_OVER) {
    printf("====== START OVER =======\n");
    /*dgnss_init(num_sdiffs, sdiffs, reciever_ecef, dt); //TODO use current baseline state*/
    return;
  }
  else if (sats_management_code == NEW_REF) {
    // do everything related to changing the reference sat here
    rebase_kf(&kf, sats_management.num_sats, &old_prns[0], &sats_management.prns[0]);
  }
}

void sdiffs_to_prns(u8 n, sdiff_t *sdiffs, u8 *prns)
{
  for (u8 i=0; i<n; i++) {
    prns[i] = sdiffs[i].prn;
  }
}

void dgnss_update_sats(u8 num_sdiffs, double reciever_ecef[3], sdiff_t *corrected_sdiffs,
                       double * dd_measurements, double dt)
{
  u8 new_prns[num_sdiffs];
  sdiffs_to_prns(num_sdiffs, corrected_sdiffs, new_prns);

  u8 old_prns[MAX_CHANNELS];
  memcpy(old_prns, sats_management.prns, sats_management.num_sats * sizeof(u8));

  if (!prns_match(&old_prns[1], num_sdiffs-1, &corrected_sdiffs[1])) {
    u8 ndx_of_intersection_in_old[sats_management.num_sats];
    u8 ndx_of_intersection_in_new[sats_management.num_sats];
    ndx_of_intersection_in_old[0] = 0;
    ndx_of_intersection_in_new[0] = 0;
    u8 num_intersection_sats = dgnss_intersect_sats(sats_management.num_sats-1, &old_prns[1],
                                                    num_sdiffs-1, &corrected_sdiffs[1],
                                                    &ndx_of_intersection_in_old[1],
                                                    &ndx_of_intersection_in_new[1]) + 1;
    reset_kf_except_state(&kf,
                          PHASE_VAR, CODE_VAR,
                          POS_TRANS_VAR, VEL_TRANS_VAR, INT_TRANS_VAR,
                          num_sdiffs, corrected_sdiffs, reciever_ecef, dt);
    if (num_intersection_sats < sats_management.num_sats) { //lost sats
      kalman_filter_state_projection(&kf,
                                     sats_management.num_sats-1,
                                     num_intersection_sats-1,
                                     &ndx_of_intersection_in_old[1]);
    }
    if (num_intersection_sats < num_sdiffs) { //gained sats
      kalman_filter_state_inclusion(&kf,
                                    num_intersection_sats-1,
                                    num_sdiffs-1,
                                    &ndx_of_intersection_in_new[1],
                                    NEW_INT_VAR);
    }

    /*print_sats_management(&sats_management);*/
    update_sats_sats_management(&sats_management, num_sdiffs-1, &corrected_sdiffs[1]);
    /*print_sats_management(&sats_management);*/

  }
  else {
    reset_kf_except_state(&kf,
                          PHASE_VAR, CODE_VAR,
                          POS_TRANS_VAR, VEL_TRANS_VAR, INT_TRANS_VAR,
                          num_sdiffs, corrected_sdiffs, reciever_ecef, dt);
  }

}

void dgnss_incorporate_observation(sdiff_t *sdiffs, double * dd_measurements,
                                   double *reciever_ecef, double dt)
{
  (void) dt;

  double ref_ecef[3];
  ref_ecef[0] = reciever_ecef[0] + kf.state_mean[0]*0.5;
  ref_ecef[1] = reciever_ecef[1] + kf.state_mean[1]*0.5;
  ref_ecef[2] = reciever_ecef[2] + kf.state_mean[2]*0.5;
  /* TODO: make a common DE and use it instead. */
  assign_decor_obs_mtx(sats_management.num_sats, sdiffs, ref_ecef,
                       kf.decor_mtx, kf.decor_obs_mtx);
  kalman_filter_update(&kf, dd_measurements);
}

void dvec_printf(double *v, u32 n)
{
  for (u32 i = 0; i < n; i++)
    printf(", %f", v[i]);
}

void dgnss_update(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3], double dt)
{
  sdiff_t corrected_sdiffs[num_sats];

  u8 old_prns[MAX_CHANNELS];
  memcpy(old_prns, sats_management.prns, sats_management.num_sats * sizeof(u8));
  //rebase globals to a new reference sat (permutes corrected_sdiffs accordingly)
  dgnss_rebase_ref(num_sats, sdiffs, reciever_ecef, dt, old_prns, corrected_sdiffs);

  double dd_measurements[2*(num_sats-1)];
  make_measurements(num_sats-1, corrected_sdiffs, dd_measurements);

  //all the added/dropped sat stuff
  dgnss_update_sats(num_sats, reciever_ecef, corrected_sdiffs, dd_measurements, dt);
  /*printf("done updating sats\n");*/
  /*MAT_PRINTF(kf.decor_obs_mtx, kf.obs_dim, kf.state_dim);*/

  // update for observation
  dgnss_incorporate_observation(corrected_sdiffs, dd_measurements, reciever_ecef, dt);

  double ref_ecef[3];
  ref_ecef[0] = reciever_ecef[0] + 0.5 * kf.state_mean[0];
  ref_ecef[1] = reciever_ecef[1] + 0.5 * kf.state_mean[1];
  ref_ecef[2] = reciever_ecef[2] + 0.5 * kf.state_mean[2];

  update_ambiguity_test(ref_ecef, PHASE_VAR, CODE_VAR, &ambiguity_test,
                        kf.state_dim, &sats_management, sdiffs,
                        kf.state_mean, kf.state_cov_U, kf.state_cov_D);

  u32 n_hyps = ambiguity_test_n_hypotheses(&ambiguity_test);
  if (n_hyps > 0) {
    printf("Num hyps: %d\n", n_hyps);
  }
}

void dgnss_float_baseline(u8 *num_used, double b[3])
{
  memcpy(b, kf.state_mean, 3 * sizeof(double));
  *num_used = sats_management.num_sats;
}

void dgnss_fixed_baseline(u8 n, sdiff_t *sdiffs, double ref_ecef[3],
                          u8 *num_used, double b[3])
{
  if (dgnss_iar_resolved()) {
    sdiff_t ambiguity_sdiffs[ambiguity_test.sats.num_sats];
    double dd_meas[2*(ambiguity_test.sats.num_sats-1)];
    make_ambiguity_dd_measurements_and_sdiffs(&ambiguity_test, n, sdiffs,
        dd_meas, ambiguity_sdiffs);
    double DE[(ambiguity_test.sats.num_sats-1) * 3];
    assign_de_mtx(ambiguity_test.sats.num_sats, ambiguity_sdiffs, ref_ecef, DE);
    hypothesis_t *hyp = (hypothesis_t*)ambiguity_test.pool->allocated_nodes_head->elem;
    print_hyp(ambiguity_test.sats.num_sats-1, hyp);
    *num_used = ambiguity_test.sats.num_sats;
    lesq_solution(n, dd_meas, hyp->N, DE, b, 0);
  } else {
    memcpy(b, kf.state_mean, 3 * sizeof(double));
    *num_used = sats_management.num_sats;
  }
}

void dgnss_reset_iar()
{
  create_ambiguity_test(&ambiguity_test);
}

s8 dgnss_iar_resolved()
{
  return ambiguity_test_n_hypotheses(&ambiguity_test) == 1;
}

kf_t * get_dgnss_kf()
{
  return &kf;
}

s32 * get_stupid_filter_ints()
{
  return stupid_state.N;
}

sats_management_t * get_sats_management()
{
  return &sats_management;
}


