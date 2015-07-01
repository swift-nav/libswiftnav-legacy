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
from observation_c cimport *

cdef extern from "libswiftnav/amb_kf.h":
  ctypedef struct nkf_t:
    u32 state_dim
    u32 obs_dim
    double amb_drift_var
    double *decor_mtx
    double *decor_obs_mtx
    double *decor_obs_cov
    double *null_basis_Q
    double *state_mean
    double *state_cov_U
    double *state_cov_D
  void predict_forward(nkf_t *kf)
  void nkf_update(nkf_t *kf, double *measurements)
  void set_nkf(nkf_t *kf, double amb_drift_var, double phase_var, double code_var, double amb_init_var,
               u8 num_sdiffs, sdiff_t *sdiffs_with_ref_first, double *dd_measurements, double ref_ecef[3])