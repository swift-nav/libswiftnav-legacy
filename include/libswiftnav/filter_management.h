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


#include "float_kf.h"

typedef struct {
  u32 state_dim;
  double mean[MAX_STATE_DIM];
  double cov_U[MAX_STATE_DIM * MAX_STATE_DIM];
  double cov_D[MAX_STATE_DIM];
} kf_state_t;

void make_measurements(u8 num_diffs, sdiff_t *sdiffs, double *raw_measurements);
void initialize_state(kf_t *kf, double *dd_measurements,
                      double pos_init_var, double vel_init_var, double int_init_var,
                      kf_state_t *state);
kf_t start_filter(double phase_var, double code_var,
                  double pos_trans_var, double vel_trans_var, double int_trans_var,
                  double pos_init_var, double vel_init_var, double int_init_var,
                  u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3], double dt,
                  kf_state_t *state);
u8 update_filter(kf_t *kf, u8 num_sats, sdiff_t *sdiffs, kf_state_t *state);
