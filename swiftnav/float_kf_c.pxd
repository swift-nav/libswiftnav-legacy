# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *

cdef extern from "libswiftnav/float_kf.h":
  s8 udu(u32 n, double *M, double *U, double *D)
  void triu(u32 n, double *M)
  void eye(u32 n, double *M)
  void reconstruct_udu(u32 n, double *U, double *D, double *M)
  void predict_forward(u32 state_dim, double *transition_mtx, double *transition_cov, double *state_mean, double *state_cov_U, double *state_cov_D)
  void update_scalar_measurement(u32 state_dim, double *h, double R,
                               double *U, double *D, double *k)
  void update_for_obs(u32 state_dim, u32 obs_dim, double *decor_obs_mtx, double *decor_obs_cov, 
                    double *intermediate_mean, double *intermediate_cov_U, double *intermediate_cov_D,
                    double *decor_obs)

