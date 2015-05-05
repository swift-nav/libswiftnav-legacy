# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.


from common cimport *
cimport memory_pool_c
cimport sats_management_c

cdef extern from "libswiftnav/ambiguity_test.h":

  ctypedef struct hypothesis_t:
    s32 N[22]
    float ll

  ctypedef struct residual_mtxs_t:
    u32 res_dim
    u8 null_space_dim
    double *null_projector
    double *half_res_cov_inv

  ctypedef struct ambiguity_test_t:
    u8 num_dds
    memory_pool_c.memory_pool_t *pool
    sats_management_c.sats_management_t sats
    residual_mtxs_t res_mtxs

  void assign_phase_obs_null_basis(u8 num_dds, double *DE_mtx, double *q)
  void assign_residual_covariance_inverse(u8 num_dds, double *obs_cov, double *q, double *r_cov_inv)
  void init_residual_matrices(residual_mtxs_t *res_mtxs, u8 num_dds, double *DE_mtx, double *obs_cov)
  void assign_r_vec(residual_mtxs_t *res_mtxs, u8 num_dds, double *dd_measurements, double *r_vec)
  void assign_r_mean(residual_mtxs_t *res_mtxs, u8 num_dds, double *hypothesis, double *r_mean)
  double get_quadratic_term(residual_mtxs_t *res_mtxs, u8 num_dds, double *hypothesis, double *r_vec)
