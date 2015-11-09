# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from almanac cimport *
from common cimport *
from gpstime cimport *
from libcpp cimport bool as bool_
from observation cimport *

cdef extern from "libswiftnav/amb_kf.h":
  enum:
    MAX_STATE_DIM
    MAX_OBS_DIM
    KF_SOS_TIMESCALE
    SOS_SWITCH

  ctypedef struct nkf_t:
    u32 state_dim
    u32 obs_dim
    double amb_drift_var
    double decor_mtx[MAX_OBS_DIM * MAX_OBS_DIM]
    double decor_obs_mtx[MAX_STATE_DIM * MAX_OBS_DIM]
    double decor_obs_cov[MAX_OBS_DIM]
    double null_basis_Q[(MAX_STATE_DIM - 3) * MAX_OBS_DIM]
    double state_mean[MAX_STATE_DIM]
    double state_cov_U[MAX_STATE_DIM * MAX_STATE_DIM]
    double state_cov_D[MAX_STATE_DIM]
    double l_sos_avg

  # ???
  void assign_phase_obs_null_basis(u8 num_dds, double *DE_mtx, double *q)
  s32 find_index_of_element_in_u8s(const u32 num_elements, const u8 x, const u8 *list)
  void rebase_covariance_udu(double *state_cov_U, double *state_cov_D, u8 num_sats, u8 *old_prns, u8 *new_prns)
  void rebase_mean_N(double *mean, const u8 num_sats, const u8 *old_prns, const u8 *new_prns)
  void rebase_covariance_sigma(double *state_cov, const u8 num_sats, const u8 *old_prns, const u8 *new_prns)
  double compute_innovation_terms(u32 state_dim, const double *h,
                                  double R, const double *U,
                                  const double *D, double *f, double *g)

  # Kalman filter stuff
  bool_ nkf_update(nkf_t *kf, const double *measurements)
  void set_nkf(nkf_t *kf, double amb_drift_var, double phase_var, double code_var, double amb_init_var,
               u8 num_sdiffs, sdiff_t *sdiffs_with_ref_first, double *dd_measurements, double ref_ecef[3])
  void set_nkf_matrices(nkf_t *kf, double phase_var, double code_var,
                        u8 num_sdiffs, sdiff_t *sdiffs_with_ref_first, double ref_ecef[3])
  void rebase_nkf(nkf_t *kf, u8 num_sats, u8 *old_prns, u8 *new_prns)
  void nkf_state_projection(nkf_t *kf,
                            u8 num_old_non_ref_sats,
                            u8 num_new_non_ref_sats,
                            u8 *ndx_of_new_sat_in_old)
  void nkf_state_inclusion(nkf_t *kf,
                           u8 num_old_non_ref_sats,
                           u8 num_new_non_ref_sats,
                           u8 *ndx_of_old_sat_in_new,
                           double *estimates,
                           double int_init_var)
  double get_sos_innov(const nkf_t *kf, const double *decor_obs)
  bool_ outlier_check(nkf_t *kf, const double *decor_obs, double *k_scalar)
  void update_kf_state(nkf_t *kf, double R, const double *f, const double *g, double alpha, double k_scalar, double innov)

cdef class KalmanFilter:
  cdef nkf_t _thisptr
