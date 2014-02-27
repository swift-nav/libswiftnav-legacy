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


#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

kf_t kf;
stupid_filter_state_t stupid_state;
sats_management_t sats_management;

void make_measurements(u8 num_diffs, sdiff_t *sdiffs, double *raw_measurements)
{
  double phase0 = sdiffs[0].carrier_phase;
  double code0 = sdiffs[0].pseudorange;
  for (u8 i=0; i<num_diffs; i++) {
    raw_measurements[i] = sdiffs[i+1].carrier_phase - phase0;
    raw_measurements[i+num_diffs] = sdiffs[i+1].pseudorange - code0;
  }
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

void dgnss_update_sats(u8 num_sats, double reciever_ecef[3], sdiff_t *corrected_sdiffs,
                       double * dd_measurements)
{
  u8 old_prns[MAX_CHANNELS];
  memcpy(old_prns, sats_management.prns, sats_management.num_sats * sizeof(u8));

  update_sats_stupid_filter(&stupid_state, sats_management.num_sats, old_prns, num_sats,
                            corrected_sdiffs, dd_measurements, reciever_ecef);
}

void dgnss_incorporate_observation(sdiff_t *sdiffs, double * dd_measurements)
{
  
  (void) sdiffs;
  (void) dd_measurements;
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
  dgnss_update_sats(num_sats, reciever_ecef, corrected_sdiffs, dd_measurements);

  // update for observation
  dgnss_incorporate_observation(corrected_sdiffs, dd_measurements);
  double b[3];
  update_stupid_filter(&stupid_state, sats_management.num_sats, corrected_sdiffs,
                        dd_measurements, b, reciever_ecef); //todo use midpoint for reference

  printf("b: %.3f %.3f %.3f\n", b[0], b[1], b[2]);
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




