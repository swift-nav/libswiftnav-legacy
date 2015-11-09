# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

"""Baseline calculations

Functions for relating the baseline vector with carrier phase
observations and ambiguities.
"""

cimport numpy as np
from fmt_utils import fmt_repr
from libc.string cimport memcpy, memset
from observation cimport *
from observation import SingleDiff
import numpy as np

DEFAULT_RAIM_THRESHOLD_ = DEFAULT_RAIM_THRESHOLD

cdef class Ambiguity:

  def __init__(self, **kwargs):
    memset(&self._thisptr, 0, sizeof(ambiguities_t))
    if kwargs:
      self._thisptr = kwargs

  def __getattr__(self, k):
    return self._thisptr.get(k)

  def __repr__(self):
    return fmt_repr(self)

  def to_dict(self):
    return self._thisptr

  def from_dict(self, d):
    self._thisptr = d


cdef class Ambiguities:

  def __init__(self, **kwargs):
    memset(&self._thisptr, 0, sizeof(ambiguities_t))
    if kwargs:
      self._thisptr = kwargs

  def __getattr__(self, k):
    return self._thisptr.get(k)

  def __repr__(self):
    return fmt_repr(self)

  def to_dict(self):
    return self._thisptr

  def from_dict(self, d):
    self._thisptr = d


def predict_carrier_obs_(N, DE, b):
  num_dds = N.shape[0]
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] N_ = np.array(N, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] DE_ = np.array(DE, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b_ = np.array(b, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] dd_obs = np.empty(num_dds, dtype=np.double)
  predict_carrier_obs(num_dds, &N_[0], &DE_[0,0], &b_[0], &dd_obs[0])
  return dd_obs

def amb_from_baseline_(DE, dd_obs, b):
  num_dds = dd_obs.shape[0]
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] DE_ = np.array(DE, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] dd_obs_ = np.array(dd_obs, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b_ = np.array(b, dtype=np.double)
  cdef np.ndarray[np.int32_t, ndim=1, mode="c"] N = np.empty(num_dds, dtype=np.int32)
  amb_from_baseline(num_dds, &DE_[0,0], &dd_obs_[0], &b_[0], &N[0])
  return N

def lesq_solution_float_(dd_obs, N, DE):
  num_dds = len(dd_obs)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] dd_obs_ = np.array(dd_obs, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] N_ = np.array(N, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] DE_ = np.array(DE, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = np.empty(3, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] resid = np.array(dd_obs, dtype=np.double)
  cdef s8 ret = lesq_solution_float(num_dds, &dd_obs_[0], &N_[0], &DE_[0, 0], &b[0], &resid[0])
  return (ret, b, resid)

def least_squares_solve_b_external_ambs_(ambs, sdiffs_with_ref_first,
                                         dd_measurements, ref_ecef,
                                         disable_raim, raim_threshold):
  num_dds = len(ambs)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ambs_ = np.array(ambs, dtype=np.double)
  cdef sdiff_t sdiffs_[32]
  mk_sdiff_array(sdiffs_with_ref_first, 32, &sdiffs_[0])
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] dd_measurements_ = np.array(dd_measurements, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(ref_ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = np.empty(3, dtype=np.double)
  cdef s8 ret = least_squares_solve_b_external_ambs(num_dds,
                                                    &ambs_[0],
                                                    &sdiffs_[0],
                                                    &dd_measurements_[0],
                                                    &ref_ecef_[0],
                                                    &b[0],
                                                    disable_raim,
                                                    raim_threshold)
  return (ret, b)

cdef mk_ambiguity_array(py_ambs, u8 n_c_ambs, ambiguity_t *c_ambs):
  if n_c_ambs < len(py_ambs):
    raise ValueError("The length of the c ambs array (" + str(n_c_ambs) + \
                      ") must be at least the length of the python ambs array " + \
                      str(len(py_ambs)) + ").")
  cdef ambiguity_t sd_
  for (i,amb) in enumerate(py_ambs):
    sd_ = (<Ambiguity> amb)._thisptr
    memcpy(&c_ambs[i], &sd_, sizeof(ambiguity_t))

def diff_ambs_(GNSSSignal ref_sid, amb_set):
  cdef u8 num_ambs = len(amb_set)
  cdef ambiguity_t amb_set_[32]
  mk_ambiguity_array(amb_set, num_ambs, &amb_set_[0])
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] dd_ambs = np.empty(num_ambs, dtype=np.double)
  diff_ambs(ref_sid._thisptr, num_ambs, &amb_set_[0], &dd_ambs[0])
  return dd_ambs

def _baseline_(sdiffs, ref_ecef, single_ambs, disable_raim, raim_threshold):
  cdef u8 num_ambs = len(single_ambs)
  cdef u8 num_sdiffs = len(sdiffs)
  cdef sdiff_t sdiffs_[32]
  mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
  cdef ambiguity_t ambs_[32]
  mk_ambiguity_array(single_ambs, 32, &ambs_[0])
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(ref_ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = np.empty(3, dtype=np.double)
  cdef u8 num_used
  cdef s8 code = baseline_(num_sdiffs, &sdiffs_[0], &ref_ecef_[0],
                           num_ambs, &ambs_[0], &num_used, &b[0],
                           disable_raim, raim_threshold)
  return (code, num_used, b)

# TODO (Buro): Seems to segfault!
def _baseline(sdiffs, ref_ecef, Ambiguities ambs, disable_raim, raim_threshold):
  cdef u8 num_sdiffs = len(sdiffs)
  cdef sdiff_t sdiffs_[32]
  mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(ref_ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = np.empty(3, dtype=np.double)
  cdef u8 num_used
  cdef s8 code = baseline(num_sdiffs, &sdiffs_[0], &ref_ecef_[0],
                          &ambs._thisptr, &num_used, &b[0],
                          disable_raim, raim_threshold)
  return (code, num_used, b)

def ambiguities_init_(Ambiguities ambs):
  ambiguities_init(&ambs._thisptr)

def lesq_solve_raim_(dd_obs, N, DE, disable_raim, raim_threshold):
  cdef u8 num_dds = len(dd_obs)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] dd_obs_ = np.array(dd_obs, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] N_ = np.array(N, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] DE_ = np.array(DE, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = np.empty(3, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] residuals = np.empty(num_dds, dtype=np.double)
  cdef u8 num_used, removed_obs
  cdef s8 code = lesq_solve_raim(num_dds, &dd_obs_[0],
                                 &N_[0], &DE_[0,0], &b[0],
                                 disable_raim, raim_threshold,
                                 &num_used, &residuals[0], &removed_obs)
  return (code, num_used, residuals, removed_obs, b)
