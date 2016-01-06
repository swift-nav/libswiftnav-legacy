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

#ifndef LIBSWIFTNAV_AMB_KF_H
#define LIBSWIFTNAV_AMB_KF_H

#include <libswiftnav/track.h>
#include <libswiftnav/almanac.h>
#include <libswiftnav/time.h>
#include <libswiftnav/common.h>
#include <libswiftnav/observation.h>
#include <libswiftnav/constants.h>

/** \addtogroup amb_kf
 * \{ */

#define MAX_STATE_DIM (MAX_CHANNELS - 1)
#define MAX_OBS_DIM (2 * MAX_CHANNELS - 5)
/** The timescale for smoothing the innovation weighted sum of squares. */
#define KF_SOS_TIMESCALE 7.0f
/** The outlier cutoff for the highpassed innovation weighted sum of squares. */
#define SOS_SWITCH 10.0f

typedef struct {
  /** The dimension of the state vector. */
  u32 state_dim;
  /** The dimension of the observation vector. */
  u32 obs_dim;
  /** The variance to use for the prediction update step (diffusion). */
  double amb_drift_var;
  /** The observation decorrelation matrix. Takes raw measurements and
   * decorrelates them. */
  double decor_mtx[MAX_OBS_DIM * MAX_OBS_DIM];
  /** The observation matrix for decorrelated measurements. */
  double decor_obs_mtx[MAX_STATE_DIM * MAX_OBS_DIM];
  /** The diagonal of the decorrelated observation covariance (for cholesky it's
   * ones). */
  double decor_obs_cov[MAX_OBS_DIM];
  /** A basis for the left nullspace usual DGNSS observatio matrix, made of
   * the DD line of sight vectors from the receivers to the sats. Used to
   * project out the baseline's influence from the observations. */
  double null_basis_Q[(MAX_STATE_DIM - 3) * MAX_OBS_DIM];
  /** The current state estimate. */
  double state_mean[MAX_STATE_DIM];
  /** The upper unit triangular U matrix of the UDU decomposition of the
   * covariance of the current state estimate. Stored dense. */
  double state_cov_U[MAX_STATE_DIM * MAX_STATE_DIM];
  /** The diagonal D matrix of the UDU decomposition of the covariance of the current
   * state estimate. Stored as a vector. */
  double state_cov_D[MAX_STATE_DIM];
  /** A moving average of the log of the weighted sum of squares innovations. */
  double l_sos_avg;
} nkf_t;

/** \} */

bool nkf_update(nkf_t *kf, const double *measurements);

void assign_phase_obs_null_basis(u8 num_dds, double *DE_mtx, double *q);
void set_nkf(nkf_t *kf, double amb_drift_var, double phase_var, double code_var, double amb_init_var,
            u8 num_sdiffs, sdiff_t *sdiffs_with_ref_first, double *dd_measurements, double ref_ecef[3]);
void set_nkf_matrices(nkf_t *kf, double phase_var, double code_var,
                     u8 num_sdiffs, sdiff_t *sdiffs_with_ref_first, double ref_ecef[3]);
s32 find_index_of_signal(const u32 num_elements, const gnss_signal_t x, const gnss_signal_t *list);
void rebase_nkf(nkf_t *kf, u8 num_sats, const gnss_signal_t *old_sids, const gnss_signal_t *new_sids);

void nkf_state_projection(nkf_t *kf,
                                    u8 num_old_non_ref_sats,
                                    u8 num_new_non_ref_sats,
                                    u8 *ndx_of_new_sat_in_old);
void nkf_state_inclusion(nkf_t *kf,
                         u8 num_old_non_ref_sats,
                         u8 num_new_non_ref_sats,
                         u8 *ndx_of_old_sat_in_new,
                         double *estimates,
                         double int_init_var);

void rebase_covariance_udu(double *state_cov_U, double *state_cov_D, u8 num_sats, const gnss_signal_t *old_sids, const gnss_signal_t *new_sids);

void rebase_mean_N(double *mean, const u8 num_sats, const gnss_signal_t *old_sids, const gnss_signal_t *new_sids);
void rebase_covariance_sigma(double *state_cov, const u8 num_sats, const gnss_signal_t *old_sids, const gnss_signal_t *new_sids);

double get_sos_innov(const nkf_t *kf, const double *decor_obs);
double compute_innovation_terms(u32 state_dim, const double *h,
                                double R, const double *U,
                                const double *D, double *f, double *g);
bool outlier_check(nkf_t *kf, const double *decor_obs, double *k_scalar);
void update_kf_state(nkf_t *kf, double R, const double *f, const double *g,
                   double alpha, double k_scalar,
                   double innov);

#endif /* LIBSWIFTNAV_AMB_KF_H */
