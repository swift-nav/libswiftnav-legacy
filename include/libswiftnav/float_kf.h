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
#include "track.h"
#include "almanac.h"
#include "gpstime.h"
#include "common.h"
#include "single_diff.h"
#include "constants.h"

#define MAX_STATE_DIM_KF (MAX_CHANNELS + 6)
#define MAX_OBS_DIM_KF (2 * MAX_CHANNELS)

s8 udu(u32 n, double *M, double *U, double *D);
void triu(u32 n, double *M);

void eye(u32 n, double *M);


typedef struct {
  u32 state_dim;
  u32 obs_dim;
  double transition_mtx[MAX_STATE_DIM_KF * MAX_STATE_DIM_KF];
  double transition_cov[MAX_STATE_DIM_KF * MAX_STATE_DIM_KF];
  double decor_mtx[MAX_OBS_DIM_KF * MAX_OBS_DIM_KF]; //the decorrelation matrix. takes raw measurements and decorrelates them
  double decor_obs_mtx[MAX_STATE_DIM_KF * MAX_OBS_DIM_KF]; //the observation matrix for decorrelated measurements
  double decor_obs_cov[MAX_OBS_DIM_KF]; //the diagonal of the decorrelated observation covariance (for cholesky is ones)
  double state_mean[MAX_STATE_DIM_KF];
  double state_cov_U[MAX_STATE_DIM_KF * MAX_STATE_DIM_KF];
  double state_cov_D[MAX_STATE_DIM_KF];
} kf_t;



void predict_forward(kf_t *kf);

//void incorporate_obs(kf_t *kf, double *decor_obs);

void update_scalar_measurement(u32 state_dim, double *h, double R,
                               double *U, double *D, double *k);
void kalman_filter_update(kf_t *kf, double *measurements);

void assign_transition_mtx(u32 state_dim, double dt, double *transition_mtx);
void assign_d_mtx(u8 num_sats, double *D);
void assign_e_mtx(u8 num_sats, sdiff_t *sats_with_ref_first, double ref_ecef[3], double *E);
void assign_de_mtx(u8 num_sats, sdiff_t *sats_with_ref_first, double ref_ecef[3], double *DE);
void assign_obs_mtx(u8 num_sats, sdiff_t *sats_with_ref_first, double ref_ecef[3], double *obs_mtx);

void assign_e_mtx_from_alms(u8 num_sats, almanac_t *alms, gps_time_t timestamp, double ref_ecef[3], double *E);
void assign_de_mtx_from_alms(u8 num_sats, almanac_t *alms, gps_time_t timestamp, double ref_ecef[3], double *DE);
void assign_obs_mtx_from_alms(u8 num_sats, almanac_t *alms, gps_time_t timestamp, double ref_ecef[3], double *obs_mtx);

void assign_decor_obs_cov(u8 num_diffs, double phase_var, double code_var,
                          double *decor_mtx, double *decor_obs_cov);
void assign_decor_obs_mtx(u8 num_sats, sdiff_t *sats_with_ref_first,
                          double ref_ecef[3], double *decor_mtx, double *obs_mtx);
void assign_decor_obs_mtx_from_alms(u8 num_sats, almanac_t *alms, gps_time_t timestamp,
                                    double ref_ecef[3], double *decor_mtx, double *obs_mtx);
void assign_transition_cov(u32 state_dim, double pos_var, double vel_var, double int_var, double *transition_cov);

void get_kf(kf_t *kf, double phase_var, double code_var,
            double pos_var, double vel_var, double int_var, 
            double pos_init_var, double vel_init_var, double int_init_var,
            u8 num_sats, sdiff_t *sats_with_ref_first, double *dd_measurements, double ref_ecef[3], double dt);
void reset_kf_except_state(kf_t *kf, 
                           double phase_var, double code_var,
                           double pos_var, double vel_var, double int_var,
                           u8 num_sats, sdiff_t *sats_with_ref_first, double ref_ecef[3], double dt);
kf_t get_kf_from_alms(double phase_var, double code_var, double pos_var, double vel_var, double int_var, 
                      u8 num_sats, almanac_t *alms, gps_time_t timestamp, double ref_ecef[3], double dt);
s32 find_index_of_element_in_u8s(u32 num_elements, u8 x, u8 *list);
void rebase_kf(kf_t *kf, u8 num_sats, u8 *old_prns, u8 *new_prns);

void kalman_filter_state_projection(kf_t *kf,
                                    u8 num_old_non_ref_sats,
                                    u8 num_new_non_ref_sats,
                                    u8 *ndx_of_new_sat_in_old);
void kalman_filter_state_inclusion(kf_t *kf,
                                   u8 num_old_non_ref_sats,
                                   u8 num_new_non_ref_sats,
                                   u8 *ndx_of_old_sat_in_new,
                                   double int_init_var);

void assign_state_rebase_mtx(u8 num_sats, u8 *old_prns, u8 *new_prns, double *rebase_mtx);
void rebase_kf(kf_t *kf, u8 num_sats, u8 *old_prns, u8 *new_prns);

#endif /* LIBSWIFTNAV_FLOAT_KF_H */

