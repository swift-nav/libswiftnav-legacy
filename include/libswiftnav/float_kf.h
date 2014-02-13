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

s8 udu(u32 n, double M[n][n], double U[n][n], double D[n]);

void triu(u32 n, double M[n][n]);

void eye(u32 n, double M[n][n]);

void reconstruct_udu(u32 n, double U[n][n], double D[n], double M[n][n]);

// typedef struct {
//   u32 state_dim;
//   u32 obs_dim;
//   double transition_mtx;
//   double transition_cov; 
//   double obs_cov_root_inv; //the decorrelation matrix. takes raw measurements and decorrelates them
//   double decor_obs_mtx; //the observation matrix for decorrelated measurements
//   double decor_obs_cov; //the diagonal of the decorrelated observation covariance (for cholesky is ones)
// } kf_t;

// void filter_update(kf_t kf,
//                    double *state_mean, double *state_cov_U, double *state_cov_D, 
//                    double *raw_measurements);

#endif /* LIBSWIFTNAV_FLOAT_KF_H */

