# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from almanac_c cimport *
from gpstime_c cimport *
from single_diff_c cimport *

cdef extern from "libswiftnav/float_kf.h":
  ctypedef struct kf_t:
    u32 state_dim
    u32 obs_dim
    double *transition_mtx
    double *transition_cov
    double *decor_mtx
    double *decor_obs_mtx
    double *decor_obs_cov
    double *state_mean
    double *state_cov_U
    double *state_cov_D
  s8 udu(u32 n, double *M, double *U, double *D)
  void triu(u32 n, double *M)
  void eye(u32 n, double *M)
  void reconstruct_udu(u32 n, double *U, double *D, double *M)
  void predict_forward(kf_t *kf)
  void update_scalar_measurement(u32 state_dim, double *h, double R,
                                 double *U, double *D, double *k)
  void incorporate_obs(kf_t *kf, double *decor_obs)
  void decorrelate(kf_t *kf, double *measurements, double * decor_measurements)
  void kalman_filter_update(kf_t *kf, double *measurements)

  void assign_transition_mtx(u32 state_dim, double dt, double *transition_mtx)
  void assign_d_mtx(u8 num_sats, double *D)

  void assign_obs_mtx(u8 num_sats, sdiff_t *sats_with_ref_first, double ref_ecef[3], double *obs_mtx)
  void assign_decor_obs_mtx(u8 num_sats, sdiff_t *sats_with_ref_first,
                          double ref_ecef[3], double *decor_mtx, double *obs_mtx)

  void assign_e_mtx_from_alms(u8 num_sats, almanac_t *alms, gps_time_t timestamp, double ref_ecef[3], double *E)
  void assign_de_mtx_from_alms(u8 num_sats, almanac_t *alms, gps_time_t timestamp, double ref_ecef[3], double *E)
  void assign_obs_mtx_from_alms(u8 num_sats, almanac_t *alms, gps_time_t timestamp, double ref_ecef[3], double *obs_mtx)
  void assign_decor_obs_cov(u8 num_diffs, double phase_var, double code_var,
                            double *decor_mtx, double *decor_obs_cov)
  void assign_decor_obs_mtx_from_alms(u8 num_sats, almanac_t *alms, gps_time_t timestamp,
                                    double ref_ecef[3], double *decor_mtx, double *obs_mtx)

  kf_t get_kf(double phase_var, double code_var, double pos_var, double vel_var, double int_var, 
            double pos_init_var, double vel_init_var, double int_init_var,
            u8 num_sats, sdiff_t *sats_with_ref_first, double *dd_measurements, double ref_ecef[3], double dt)
  kf_t get_kf_from_alms(double phase_var, double code_var, double pos_var, double vel_var, double int_var, 
                      u8 num_sats, almanac_t *alms, gps_time_t timestamp, double ref_ecef[3], double dt)

  void least_squares_solve(kf_t *kf, double *measurements, double *lsq_state)

