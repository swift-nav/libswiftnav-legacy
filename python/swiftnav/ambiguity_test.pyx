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
from fmt_utils import fmt_repr
from libc.string cimport memcpy, memset
from observation cimport *
from sats_management cimport *
import numpy as np

cdef hypothesis_to_tuple(u8 num_dds, hypothesis_t *h):
  ambs = [ h.N[i] for i in range(num_dds) ]
  return (h.ll, ambs)

cdef class Hypothesis:

  def __cinit__(self,
                np.ndarray[np.int32_t, ndim=1, mode="c"] N,
                float ll):
    memset(&self._thisptr, 0, sizeof(hypothesis_t))
    self._thisptr.N = N
    self._thisptr.ll = ll

  def __getattr__(self, k):
    return self._thisptr.get(k)

  def __repr__(self):
    return fmt_repr(self)

  def to_dict(self):
    return self._thisptr

  def from_dict(self, d):
    self._thisptr = d

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

  def __init__(self,
               SatsManagement sats or None):
    memset(&self._thisptr, 0, sizeof(ambiguity_test_t))
    create_ambiguity_test(&self._thisptr)
    if sats:
      self._thisptr.sats = sats._thisptr

  # def __dealloc__(self):
  #   destroy_ambiguity_test(&self._thisptr)

  property sats:
    def __get__(self):
     return self._thisptr.sats

  def reset(self):
    raise NotImplementedError("No implementation in C file")

  def sets_match(self, sdiffs):
    num_sdiffs = len(sdiffs)
    cdef sdiff_t sdiffs_[32]
    mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
    return sats_match(&self._thisptr, num_sdiffs, &sdiffs_[0])

  def update_reference(self, sdiffs, sdiffs_with_ref_first):
    num_sdiffs = len(sdiffs)
    cdef sdiff_t sdiffs_[32], sdiffs_with_ref_first_[32]
    mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
    mk_sdiff_array(sdiffs_with_ref_first, 32, &sdiffs_with_ref_first_[0])
    return ambiguity_update_reference(&self._thisptr, num_sdiffs,
                                      &sdiffs_[0],
                                      &sdiffs_with_ref_first_[0])

  def update_uninimous_ambiguities(self):
    update_unanimous_ambiguities(&self._thisptr)

  def n_hypotheses(self):
    return ambiguity_test_n_hypotheses(&self._thisptr)

  def pool_contains(self, ambs):
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] ambs_ = np.array(ambs, dtype=np.double)
    return ambiguity_test_pool_contains(&self._thisptr, &ambs_[0])

  def pool_ll(self, ambs):
    num_ambs = len(ambs)
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] ambs_ = np.array(ambs, dtype=np.double)
    return ambiguity_test_pool_ll(&self._thisptr, num_ambs, &ambs_[0])

  def pool_prob(self, ambs):
    num_ambs = len(ambs)
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] ambs_ = np.array(ambs, dtype=np.double)
    return ambiguity_test_pool_prob(&self._thisptr, num_ambs, &ambs_[0])

  def MLE_ambs(self):
    cdef np.ndarray[np.int32_t, ndim=1, mode="c"] ambs_ = np.empty(self._thisptr.num_dds, dtype=np.int32)
    ambiguity_test_MLE_ambs(&self._thisptr, &ambs_[0])
    return ambs_

  def test_ambiguities(self, ambiguity_dd_measurements):
    """Updates the IAR hypothesis pool log likelihood ratios and filters
   them.  It assumes that the observations are structured to match the
   amb_test sats.  INVALIDATES unanimous ambiguities.

    """
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] ambs_ \
      = np.array(ambiguity_dd_measurements, dtype=np.double)
    test_ambiguities(&self._thisptr, &ambs_[0])

  def ambiguity_update_sats(self,
                            sdiffs,
                            SatsManagement float_sats or None,
                            np.ndarray[np.double_t, ndim=1, mode="c"] float_mean or None,
                            np.ndarray[np.double_t, ndim=2, mode="c"] float_cov_U or None,
                            np.ndarray[np.double_t, ndim=2, mode="c"] float_cov_D or None,
                            is_bad_measurement):
    num_sdiffs = len(sdiffs)
    cdef sdiff_t sdiffs_[32]
    print "here"
    mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
    print sdiffs_
    return ambiguity_update_sats(&self._thisptr,
                                 num_sdiffs,
                                 &sdiffs_[0],
                                 &float_sats._thisptr,
                                 &float_mean[0] if float_mean else NULL,
                                 &float_cov_U[0, 0] if float_cov_U else NULL,
                                 &float_cov_D[0, 0] if float_cov_D else NULL,
                                 is_bad_measurement)

  def find_indices_of_intersection_sats(self, sdiffs_with_ref_first):
    num_sdiffs = len(sdiffs_with_ref_first)
    cdef sdiff_t sdiffs_with_ref_first_[32]
    mk_sdiff_array(sdiffs_with_ref_first, 32, &sdiffs_with_ref_first_[0])
    cdef np.ndarray[np.uint8_t, ndim=1, mode="c"] intersection_ndxs = np.empty(num_sdiffs, dtype=np.uint8)
    cdef u8 k = find_indices_of_intersection_sats(&self._thisptr, num_sdiffs,
                                                  &sdiffs_with_ref_first_[0], &intersection_ndxs[0])
    return (k, intersection_ndxs)

  def iar_can_solve(self):
    return ambiguity_iar_can_solve(&self._thisptr)

  # def make_ambiguity_dd_measurements_and_sdiffs(self, sdiffs, ambiguity_dd_measurements, amb_sdiffs):
  #   num_sdiffs = len(sdiffs)
  #   cdef sdiff_t sdiffs_[32], amb_sdiffs[32]
  #   mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
  #   cdef c8 code = make_ambiguity_dd_measurements_and_sdiffs(&self._thisptr,
  #                                                            num_sdiffs,
  #                                                            &sdiffs_[0],
  #                                                            &ambiguity_dd_measurements[0],
  #                                                            &amb_sdiffs[0])
  #   return (code, ambiguity_dd_measurements, amb_sdiffs)

  def sat_projection(self, dd_intersection_ndxs):
    num_dds_intersection = len(dd_intersection_ndxs)
    cdef np.ndarray[np.uint8_t, ndim=1, mode="c"] dd_intersection_ndxs_ = np.array(dd_intersection_ndxs, dtype=np.uint8)
    return ambiguity_sat_projection(&self._thisptr, num_dds_intersection, &dd_intersection_ndxs_[0])

  def sat_inclusion(self, num_dds_in_intersection, SatsManagement float_sats, float_mean, float_cov_U, float_cov_D):
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] float_mean_ = np.array(float_mean, dtype=np.double)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] float_cov_U_ = np.matrix(float_cov_U, dtype=np.double)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] float_cov_D_ = np.matrix(float_cov_U, dtype=np.double)
    cdef u8 code = ambiguity_sat_inclusion(&self._thisptr,
                                           num_dds_in_intersection,
                                           &float_sats._thisptr,
                                           &float_mean_[0],
                                           &float_cov_U_[0, 0],
                                           &float_cov_D_[0, 0])
    return code


def float_to_decor_(addible_float_cov, addible_float_mean, num_addible_dds, num_dds_to_add,
                    lower_bounds, upper_bounds, Z, Z_inv):
  raise NotImplementedError

cdef class ResidualMatrices:

  def __init__(self,
               np.ndarray[np.double_t, ndim=2, mode="c"] DE_mtx,
               np.ndarray[np.double_t, ndim=2, mode="c"] obs_cov):
    self.num_dds = obs_cov.shape[0] / 2
    init_residual_matrices(&self._thisptr, self.num_dds, &DE_mtx[0,0], &obs_cov[0,0])

  def get_r_vec(self, np.ndarray[np.double_t, ndim=1, mode="c"] dd_measurements):
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] r_vec = np.empty(2*self.num_dds-3, dtype=np.double)
    assign_r_vec(&self._thisptr, self.num_dds, &dd_measurements[0], &r_vec[0])
    return r_vec

  def get_r_mean(self, np.ndarray[np.double_t, ndim=1, mode="c"] hypothesis):
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] r_mean = np.empty(2*self.num_dds-3, dtype=np.double)
    assign_r_mean(&self._thisptr, self.num_dds, &hypothesis[0], &r_mean[0])
    return r_mean

  def get_quadratic_term(self,
                         np.ndarray[np.double_t, ndim=1, mode="c"] hypothesis,
                         np.ndarray[np.double_t, ndim=1, mode="c"] r_vec):
    return get_quadratic_term(&self._thisptr, self.num_dds, &hypothesis[0], &r_vec[0])


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
