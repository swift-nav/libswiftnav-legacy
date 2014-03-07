# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.


from common cimport *

cdef extern from "libswiftnav/ambiguity_test.h":
  ctypedef struct residual_mtxs_t:
    u32 res_dim
    u8 null_space_dim
    double *null_projector
    double *res_cov_inverse

  void assign_phase_obs_null_basis(u8 num_dds, double *DE_mtx, double *q)
