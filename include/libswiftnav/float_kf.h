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

#ifndef LIBSWIFTNAV_FLOAT_KF_H
#define LIBSWIFTNAV_FLOAT_KF_H

#include "common.h"

#define MAX_SATS 15
#define MAX_STATE_DIM (MAX_SATS + 6)
#define MAX_OBS_DIM (2 * MAX_SATS)

s8 udu(u32 n, double *M, double *U, double *D);

void triu(u32 n, double *M);

void eye(u32 n, double *M);

void reconstruct_udu(u32 n, double *U, double *D, double *M);

typedef struct {
  u32 state_dim;
  u32 obs_dim;
  double transition_mtx[MAX_STATE_DIM * MAX_STATE_DIM];
  double transition_cov[MAX_STATE_DIM * MAX_STATE_DIM]; 
  double obs_cov_root_inv[MAX_OBS_DIM * MAX_OBS_DIM]; //the decorrelation matrix. takes raw measurements and decorrelates them
  double decor_obs_mtx[MAX_OBS_DIM * MAX_OBS_DIM]; //the observation matrix for decorrelated measurements
  double decor_obs_cov[MAX_OBS_DIM]; //the diagonal of the decorrelated observation covariance (for cholesky is ones)
} kf_t;

void predict_forward(u32 state_dim, double *transition_mtx, double *transition_cov,
                     double *state_mean, double *state_cov_U, double *state_cov_D);

void update_for_obs(u32 state_dim, u32 obs_dim, double *decor_obs_mtx, double *decor_obs_cov, 
                    double *intermediate_mean, double *intermediate_cov_U, double *intermediate_cov_D,
                    double *decor_obs);

void update_scalar_measurement(u32 state_dim, double *h, double R,
                               double *U, double *D, double *k);

void filter_update(kf_t *kf,
                   double *state_mean, double *state_cov_U, double *state_cov_D, 
                   double *raw_measurements);

#endif /* LIBSWIFTNAV_FLOAT_KF_H */

