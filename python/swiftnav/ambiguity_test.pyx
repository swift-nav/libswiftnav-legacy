# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

cimport memory_pool
cimport numpy as np
from common cimport *
import numpy as np

cdef hypothesis_to_tuple(u8 num_dds, hypothesis_t *h):
  ambs = [ h.N[i] for i in range(num_dds) ]
  return (h.ll, ambs)

cdef class Hypothesis:

  def __cinit__(self):
    raise NotImplementedError

cdef class UnanimousAmbiguityCheck:

  def __init__(self):
    raise NotImplementedError

cdef class IntersectionCount:

  def __init__(self):
    raise NotImplementedError

  def print_state(self):
    raise NotImplementedError

cdef class GenerateHypothesisState_t2:

  def __init__(self):
    raise NotImplementedError

cdef class AmbiguityTest:

  def __cinit__(self):
    raise NotImplementedError
    #void create_empty_ambiguity_test(ambiguity_test_t *amb_test)
    #void create_ambiguity_test(ambiguity_test_t *amb_test)

  def __dealloc__(self):
    raise NotImplementedError

  def __iter__(self):
    raise NotImplementedError
    # cdef memory_pool.node_t *currp = self.test.pool.allocated_nodes_head
    # while currp is not NULL:
    #   yield hypothesis_to_tuple(self.test.sats.num_sats - 1, <hypothesis_t *>currp.elem)
    #   currp = currp.hdr.next

  def reset(self):
    raise NotImplementedError

  def sets_match(self):
    raise NotImplementedError

  def update_reference(self):
    raise NotImplementedError

  def update_uninimous_ambiguities(self):
    raise NotImplementedError

  def n_hypotheses(self):
    raise NotImplementedError

  def pool_contains(self, ambs):
    raise NotImplementedError

  def pool_ll(self, ambs):
    raise NotImplementedError

  def pool_prob(self, ambs):
    raise NotImplementedError

  def MLE_ambs(self, ambs):
    raise NotImplementedError

  def test_ambiguities(self, ambiguity_dd_measurements):
    raise NotImplementedError

  def ambiguity_update_sats(self, sdiffs, float_sats, float_mean, float_cov_U, float_cov_D, is_bad_measurement):
    raise NotImplementedError

  def find_indices_of_intersection_sats(self, sdiffs_with_ref_first, intersection_ndxs):
    raise NotImplementedError

  def iar_can_solve(self):
    raise NotImplementedError

  def make_dd_measurements_and_sdiffs(self, sdiffs, ambiguity_dd_measurements, amb_sdiffs):
    raise NotImplementedError

  def sat_projection(self, dd_intersection_ndxs):
    raise NotImplementedError

  def sat_inclusion(self, num_dds_in_intersection, float_sats, float_mean, float_cov_U, float_cov_D):
    raise NotImplementedError

def float_to_decor_(addible_float_cov, addible_float_mean, num_addible_dds, num_dds_to_add,
                    lower_bounds, upper_bounds, Z, Z_inv):
  raise NotImplementedError

cdef class ResidualMatrices:

  def __init__(self,
               np.ndarray[np.double_t, ndim=2, mode="c"] DE_mtx,
               np.ndarray[np.double_t, ndim=2, mode="c"] obs_cov):
    self.num_dds = obs_cov.shape[0] / 2
    init_residual_matrices(&(self._thisptr), self.num_dds, &DE_mtx[0,0], &obs_cov[0,0])


  def get_r_vec(self, np.ndarray[np.double_t, ndim=1, mode="c"] dd_measurements):
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] r_vec = np.empty(2*self.num_dds-3, dtype=np.double)
    assign_r_vec(&(self._thisptr), self.num_dds, &dd_measurements[0], &r_vec[0])
    return r_vec

  def get_r_mean(self, np.ndarray[np.double_t, ndim=1, mode="c"] hypothesis):
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] r_mean = np.empty(2*self.num_dds-3, dtype=np.double)
    assign_r_mean(&(self._thisptr), self.num_dds, &hypothesis[0], &r_mean[0])
    return r_mean

  def get_quadratic_term(self,
                         np.ndarray[np.double_t, ndim=1, mode="c"] hypothesis,
                         np.ndarray[np.double_t, ndim=1, mode="c"] r_vec):
    return get_quadratic_term(&(self._thisptr), self.num_dds, &hypothesis[0], &r_vec[0])


# def get_phase_obs_null_basis(DE):
#   num_dds = DE.shape[0]
#   cdef np.ndarray[np.double_t, ndim=2, mode="c"] DE_ = np.array(DE, dtype=np.double)
#   cdef np.ndarray[np.double_t, ndim=2, mode="c"] q_ = np.empty((max(num_dds-3,0),num_dds), dtype=np.double)
#   assign_phase_obs_null_basis(DE.shape[0], &DE_[0,0], &q_[0,0])
#   return q_

# def get_residual_covariance_inverse(obs_cov, q):
#   cdef u8 num_dds = obs_cov.shape[0]/2
#   cdef np.ndarray[np.double_t, ndim=2, mode="c"] q_ = np.array(q, dtype=np.double)
#   cdef np.ndarray[np.double_t, ndim=2, mode="c"] obs_cov_ = np.array(obs_cov, dtype=np.double)
#   cdef np.ndarray[np.double_t, ndim=2, mode="c"] r_cov_inv_ = np.empty((2*num_dds-3, 2*num_dds-3), dtype=np.double)
#   assign_residual_covariance_inverse(num_dds, &obs_cov_[0,0], &q_[0,0], &r_cov_inv_[0,0])
#   return r_cov_inv_*2
