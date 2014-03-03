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


#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

kf_t kf;
stupid_filter_state_t stupid_state;
sats_management_t sats_management;

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
    printf("len_doesn't match\n");
    return false;
  }
  for (u8 i=0; i<num_non_ref_sdiffs; i++) {
    //iterate through the non-reference_sats
    printf("old[%u]=%u, new[%u]=%u\n", i+1, old_non_ref_prns[i], i+1, non_ref_sdiffs[i].prn);
    if (old_non_ref_prns[i] != non_ref_sdiffs[i].prn) {
      return false;
    }
  }
  printf("prns_match\n");
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

  kf = get_kf(PHASE_VAR, CODE_VAR,
              POS_TRANS_VAR, VEL_TRANS_VAR, INT_TRANS_VAR,
              POS_INIT_VAR,  VEL_INIT_VAR,  INT_INIT_VAR,
              num_sats, corrected_sdiffs, dd_measurements, reciever_ecef, dt);

  /*double b_init[3] = {0, 0, 0}; // Zero baseline*/
  double b_init[3] = {-1.4861289 ,  0.84761746, -0.01029364}; // colin's piksi data
  // double b_init[3] = {1.02571973, -0.15447333, 0.81029273}; // The antenna tree
  // double b_init[3] = {-1.02571973, 0.15447333, -0.81029273}; // The antenna tree, switched
  init_stupid_filter(&stupid_state, num_sats, corrected_sdiffs, dd_measurements, b_init, reciever_ecef);

}

void dgnss_rebase_ref(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3], double dt, u8 old_prns[MAX_CHANNELS], sdiff_t *corrected_sdiffs)
{
  //all the ref sat stuff
  s8 sats_management_code = rebase_sats_management(&sats_management, num_sats, sdiffs, corrected_sdiffs);
  if (sats_management_code == NEW_REF_START_OVER) {
    dgnss_init(num_sats, sdiffs, reciever_ecef, dt); //TODO use current baseline state
    return;
  }
  else if (sats_management_code == NEW_REF) {
    // do everything related to changing the reference sat here
    rebase_stupid_filter(&stupid_state, sats_management.num_sats, &old_prns[0], &sats_management.prns[0]);
    // rebase_kf(&kf, sats_management.num_sats, &old_prns[0], &sats_management.prns[0]); //TODO implement correctly
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
    MAT_PRINTF(kf.decor_obs_mtx, kf.obs_dim, kf.state_dim);
    if (num_intersection_sats < sats_management.num_sats) { //lost sats
      printf("kf_projection\n");
      VEC_PRINTF(kf.state_mean, kf.state_dim);
      kalman_filter_state_projection(&kf,
                                     sats_management.num_sats-1,
                                     num_intersection_sats-1,
                                     &ndx_of_intersection_in_old[1]);
      VEC_PRINTF(kf.state_mean, kf.state_dim);
    }
    if (num_intersection_sats < num_sdiffs) { //gained sats
      printf("kf_inclusion\n");
      kalman_filter_state_inclusion(&kf,
                                    num_intersection_sats-1,
                                    num_sdiffs-1,
                                    &ndx_of_intersection_in_new[1],
                                    INT_INIT_VAR,
                                    dd_measurements);
    }

    update_sats_stupid_filter(&stupid_state, sats_management.num_sats, old_prns, num_sdiffs,
                            corrected_sdiffs, dd_measurements, reciever_ecef);
    print_sats_management(&sats_management);
    update_sats_sats_management(&sats_management, num_sdiffs-1, &corrected_sdiffs[1]);
    print_sats_management(&sats_management);
  }
MAT_PRINTF(kf.decor_obs_mtx, kf.obs_dim, kf.state_dim);
}

void dgnss_incorporate_observation(sdiff_t *sdiffs, double * dd_measurements, double *reciever_ecef)
{
  double ref_ecef[3];
  ref_ecef[0] = reciever_ecef[0] + kf.state_mean[0]*0.5;
  ref_ecef[1] = reciever_ecef[1] + kf.state_mean[1]*0.5;
  ref_ecef[2] = reciever_ecef[2] + kf.state_mean[2]*0.5;
  assign_decor_obs_mtx(sats_management.num_sats, sdiffs, ref_ecef, kf.decor_mtx, kf.decor_obs_mtx); //TODO make a common DE and use it instead
  kalman_filter_update(&kf, dd_measurements);

  double b[3];
  update_stupid_filter(&stupid_state, sats_management.num_sats, sdiffs,
                        dd_measurements, b, ref_ecef);
  printf("b: %.3f %.3f %.3f\n", b[0], b[1], b[2]);
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
  printf("done updating sats\n");
  MAT_PRINTF(kf.decor_obs_mtx, kf.obs_dim, kf.state_dim);
  // update for observation
  dgnss_incorporate_observation(corrected_sdiffs, dd_measurements, reciever_ecef);
  
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




