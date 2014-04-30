# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

cimport ambiguity_test_c
import numpy as np
cimport numpy as np
from common cimport *

def get_phase_obs_null_basis(DE):
  num_dds = DE.shape[0]
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] DE_ = \
    np.array(DE, dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] q_ = \
    np.empty((max(num_dds-3,0),num_dds), dtype=np.double)

  ambiguity_test_c.assign_phase_obs_null_basis(DE.shape[0], &DE_[0,0], &q_[0,0])
  return q_

def get_residual_covariance_inverse(obs_cov, q):
  cdef u8 num_dds = obs_cov.shape[0]/2

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] q_ = \
    np.array(q, dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] obs_cov_ = \
    np.array(obs_cov, dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] r_cov_inv_ = \
    np.empty((2*num_dds-3, 2*num_dds-3), dtype=np.double)

  ambiguity_test_c.assign_residual_covariance_inverse(num_dds, &obs_cov_[0,0], &q_[0,0], &r_cov_inv_[0,0])
  return r_cov_inv_*2

cdef class ResidualMtxs:
  cdef ambiguity_test_c.residual_mtxs_t residual_mtxs
  cdef u8 num_dds

  def __init__(self,
               np.ndarray[np.double_t, ndim=2, mode="c"] DE_mtx, 
               np.ndarray[np.double_t, ndim=2, mode="c"] obs_cov):

    self.num_dds = obs_cov.shape[0] / 2

    ambiguity_test_c.init_residual_matrices(&(self.residual_mtxs),
                                             self.num_dds,
                                             &DE_mtx[0,0],
                                             &obs_cov[0,0])


  def get_r_vec(self, np.ndarray[np.double_t, ndim=1, mode="c"] dd_measurements):
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] r_vec = \
      np.empty(2*self.num_dds-3, dtype=np.double)
    ambiguity_test_c.assign_r_vec(&(self.residual_mtxs), self.num_dds, &dd_measurements[0], &r_vec[0])
    return r_vec

  def get_r_mean(self, np.ndarray[np.double_t, ndim=1, mode="c"] hypothesis):
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] r_mean = \
      np.empty(2*self.num_dds-3, dtype=np.double)
    ambiguity_test_c.assign_r_mean(&(self.residual_mtxs), self.num_dds, &hypothesis[0], &r_mean[0])
    return r_mean

  def get_quadratic_term(self, np.ndarray[np.double_t, ndim=1, mode="c"] hypothesis,
                         np.ndarray[np.double_t, ndim=1, mode="c"] r_vec):
    return ambiguity_test_c.get_quadratic_term(&(self.residual_mtxs), self.num_dds, &hypothesis[0], &r_vec[0])

# void assign_r_vec(residual_mtxs_t *res_mtxs, u8 num_dds, double *dd_measurements, double *r_vec)
#   void assign_r_mean(residual_mtxs_t *res_mtxs, u8 num_dds, double *hypothesis, double *r_mean)
#   double get_quadratic_term(residual_mtxs_t *res_mtxs, u8 num_dds, double *hypothesis, double *r_vec)