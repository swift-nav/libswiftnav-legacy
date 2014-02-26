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
#include "sats_management.h"


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
  double dd_measurements[2*(num_sats-1)];
  make_measurements(num_sats-1, sdiffs, dd_measurements);

  sdiff_t sdiffs_with_ref_first[num_sats];
  init_sats_management(&sats_management, num_sats, sdiffs, sdiffs_with_ref_first);
/*
  kf = get_kf(PHASE_VAR, CODE_VAR,
              POS_TRANS_VAR, VEL_TRANS_VAR, INT_TRANS_VAR,
              POS_INIT_VAR,  VEL_INIT_VAR,  INT_INIT_VAR,
              num_sats, sdiffs_with_ref_first, dd_measurements, reciever_ecef, dt);
*/
  (void) dt;
  double b_init[3] = {1.02571973, -0.15447333, 0.81029273}; // The antenna tree
  init_stupid_filter(&stupid_state, num_sats, sdiffs_with_ref_first, dd_measurements, b_init, reciever_ecef);
}

void dgnss_update(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3], double dt)
{
  sdiff_t sdiffs_with_ref_first[num_sats];
  
  //all the ref sat stuff
  u8 old_prns[MAX_CHANNELS];
  memcpy(old_prns, sats_management.prns, sats_management.num_sats * sizeof(u8));
  s8 sats_management_code = update_sats_management(&sats_management, num_sats, sdiffs, sdiffs_with_ref_first);
  if (sats_management_code == NEW_REF_START_OVER) {
    dgnss_init(num_sats, sdiffs, reciever_ecef, dt);
  }
  else if (sats_management_code == NEW_REF) {
    // do everything related to changing the reference sat here
    rebase_stupid_filter(&stupid_state, sats_management.num_sats, &old_prns[0], &sats_management.prns[0]);
    // rebase_kf(&kf, sats_management.num_sats, &old_prns[0], &sats_management.prns[0]);
  }

  double dd_measurements[2*(num_sats-1)];
  make_measurements(num_sats-1, sdiffs_with_ref_first, dd_measurements);

  //all the changed sat stuff
  update_sats_stupid_filter(&stupid_state, sats_management.num_sats, old_prns, num_sats,
                            sdiffs_with_ref_first, dd_measurements, reciever_ecef);

  // update for observation
  double b[3];
  update_stupid_filter(&stupid_state, sats_management.num_sats, sdiffs_with_ref_first,
                        dd_measurements, b, reciever_ecef); //todo use midpoint for reference
}




