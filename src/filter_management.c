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
#include <
#include "float_kf.h"
#include "single_diff.h"


#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

void make_measurements(u8 num_diffs, sdiff_t *sdiffs, double *raw_measurements)
{
  double phase0 = sdiffs[0].carrier_phase;
  double code0 = sdiffs[0].pseudorange;
  for (u8 i=0; i<num_diffs; i++) {
    raw_measurements[i] = sdiffs[i+1].carrier_phase - phase0;
    raw_measurements[i+num_sats] = sdiffs[i+1].pseudorange - code0;
  }
}

void initialize_state(kf_t *kf, double *dd_measurements,
                      double pos_init_var, double vel_init_var, double int_init_var,
                      state_t *state)
{
  double lsq_solution[MAX(kf->obs_dim,kf->state_dim)];
  least_squares_solve(kf, dd_measurements, lsq_solution);
  memcpy(state->mean, lsq_solution, kf->state_dim * sizeof(double));
  eye(kf->state_dim, state->cov_U);
  state->cov_D[0] = pos_init_var;
  state->cov_D[1] = pos_init_var;
  state->cov_D[2] = pos_init_var;
  state->cov_D[3] = vel_init_var;
  state->cov_D[4] = vel_init_var;
  state->cov_D[5] = vel_init_var;
  for (u32 i=6; i<kf->state_dim; i++) {
    state->cov_D[i] = int_init_var;
  }
}

kf_t start_filter(double phase_var, double code_var,
                  double pos_trans_var, double vel_trans_var, double int_trans_var,
                  double pos_init_var, double vel_init_var, double int_init_var,
                  u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3], dt,
                  state_t *state)
{
  kf_t kf = get_kf(phase_var, code_var, pos_trans_var, vel_trans_var, int_trans_var,
                num_sats, sdiffs, reciever_ecef, dt);
  double dd_measurements[kf.obs_dim];
  make_measurements(num_sats-1, sdiffs, dd_measurements);
  initialize_state(&kf, &dd_measurements[0],
                   pos_init_var, vel_init_var, int_init_var,
                   state);
  return kf;
}


bool contains(u8 num_sats, u8 ref_prn, sdiff_t *sdiffs) 
{
  //todo make and use assumptions about the ordering of the prns
  for (u8 i=0; i<num_sats; i++) {
    if (ref_prn == sdiffs[i].prn) {
      return true;
    }
  }
  return false;
}

u8 intersect_sats(u8 num_sats1, u8 num_sats2, u8 *sats1, sdiff_t *sdiffs)
{
  //todo make and use assumptions about the ordering of the prns
  sdiff_t intersection_sats[MAX(num_sats1, num_sats2)];
  u8 i=0;
  for (u8 j=0; j<num_sats1; j++) {
    u8 prn1 = sats1[j];
    for (u8 k=0; k<num_sats2; k++) {
      u8 prn2 = sdiffs[k].prn;
      if (prn1 == prn2) {
        intersect_sats[i] = sdiffs[k];
        k = num_sats2;
        i +=1;
      }
    }
  }
  return i;
}

u8 choose_reference_sat(u8 num_sats, sdiff_t *sats)
{
  double best_snr=sats[0].snr;
  u8 best_prn=sats[0].prn;
  for (u8 i=1; i<num_sats; i++) {
    if (sats[i].snr > best_snr) {
      best_snr = sats[i].snr;
      best_prn = sats[i].prn;
    }
  }
  return best_prn;
}

s8 check_and_rebase_reference_sat(kf_t *kf, u8 num_sats, sdiff_t *sdiffs, state_t *state)
{
  if !contains(num_sats, kf->prns_with_ref_first[0], sdiffs) {
    u8 common_sats[num_sats];
    u8 intersection_size = intersect_sats(num_sats, kf->prns_with_ref_first, sdiffs, common_sats);
    if (intersection_size>=4) {
      u8 new_ref_prn = choose_reference_sat(common_sats);
    }
  }
}

u8 update_filter(kf_t *kf, u8 num_sats, sdiff_t *sdiffs, state_t *state)
{
  //TODO add update_E logic

  u8 control_logic = check_and_rebase_reference_sat(kf, num_sats, sdiffs, state);
}


