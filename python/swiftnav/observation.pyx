# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

# cython: embedsignature=True

cimport numpy as np
from fmt_utils import fmt_repr
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy, memset
from time cimport *
from ephemeris cimport *
from signal cimport *
from track cimport *
import numpy as np

"""Single Difference Observations

Functions for storing and manipulating single difference observations.
"""

cdef class SingleDiff:

  def __init__(self, **kwargs):
    memset(&self._thisptr, 0, sizeof(sdiff_t))
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

def single_diff_(m_a, m_b):
  """Calculate single differences from two sets of observations.
  Undifferenced input observations are assumed to be both taken at the
  same time.

  Parameters
  ----------

  Returns
  -------

  """
  n_a = len(m_a)
  n_b = len(m_b)
  cdef navigation_measurement_t nav_msg_a[32], nav_msg_b[32]
  mk_nav_meas_array(m_a, n_a, &nav_msg_a[0])
  mk_nav_meas_array(m_b, n_b, &nav_msg_a[0])
  n_sdiffs = min(n_a, n_b)
  cdef sdiff_t sdiffs[32]
  cdef u8 count = single_diff(n_a, nav_msg_a, n_b, nav_msg_b, &sdiffs[0])
  return read_sdiff_array(count, &sdiffs[0])

# TODO (Buro): special note that there are a bunch of statically
# declared things in the .c file for observation.c

def make_propagated_sdiffs_(m_local, m_remote, remote_dists, remote_pos_ecef,
                            es, GpsTime t):
  """Propagates remote measurements to a local time and makes
  sdiffs. When we get two sets of observations that aren't time
  matched to each other (but are internally time matched within each
  set), we need to adjust one set of measurements to be our best guess
  of what it would have been had we measured it at the other set's
  time. This function does that and differences those measurements
  from sats present in both sets.

  Parameters
  ----------

  Returns
  -------

  """
  n_local = len(m_local)
  n_remote = len(m_remote)
  cdef navigation_measurement_t m_local_[32], m_remote_[32]
  mk_nav_meas_array(m_local, n_local, &m_local_[0])
  mk_nav_meas_array(m_remote, n_remote, &m_remote_[0])
  cdef ephemeris_t** es_ = <ephemeris_t **>malloc(len(es)*sizeof(ephemeris_t *))
  mk_ephemeris_array(es, len(es), es_)
  # TODO (Buro): These can be moved up into type hints.
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] remote_dists_ \
    = np.array(remote_dists, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] remote_pos_ecef_ \
    = np.array(remote_pos_ecef, dtype=np.double)
  cdef sdiff_t sdiffs_[32]
  n_sdiffs = min(n_local, n_remote)
  cdef u8 n_sats =  make_propagated_sdiffs(n_local, &m_local_[0],
                                           n_remote, &m_remote_[0],
                                           &remote_dists_[0],
                                           &remote_pos_ecef_[0],
                                           es_,
                                           &t._thisptr,
                                           &sdiffs_[0])
  free(es_)
  return read_sdiff_array(n_sats, &sdiffs_[0])

def make_dd_measurements_and_sdiffs_(GNSSSignal ref_sid, non_ref_sids, sdiffs_in):
  """Constructs the double differenced measurements and sdiffs needed
  for IAR. This requires that all IAR prns are in the sdiffs used
  (sdiffs is a superset of IAR's PRNs), that the sdiffs are ordered by
  prn, and that the IAR non_ref_prns are also ordered (both ascending).

  Parameters
  ----------

  Returns
  -------

  """
  num_sdiffs = len(sdiffs_in)
  num_dds = len(non_ref_sids)
  cdef sdiff_t sdiffs_in_[32], sdiffs_out_[32]
  cdef gnss_signal_t non_ref_sids_[32]
  mk_signal_array(non_ref_sids, num_dds, &non_ref_sids_[0])
  mk_sdiff_array(sdiffs_in, num_sdiffs, &sdiffs_in_[0])
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] dd_meas = np.empty(num_dds, dtype=np.double)
  cdef s8 code = make_dd_measurements_and_sdiffs(ref_sid._thisptr,
                                                 &non_ref_sids_[0],
                                                 num_dds, num_sdiffs,
                                                 &sdiffs_in_[0],
                                                 &dd_meas[0], &sdiffs_out_[0])
  return (code, dd_meas, read_sdiff_array(num_dds + 1, &sdiffs_out_[0]))

def check_lock_counters_(sds, lock_counters):
  """Checks to see if any satellites have had their lock counter values
  have changed.

  If the lock counter changes it indicates that the satellite should
  be reinitialized in the filter. This function checks the current
  lock counter values in a set of single difference observations with
  a stored state (`lock_counters`).

  Parameters
  ----------

  Returns
  -------

  """
  n_sdiffs = len(sds)
  cdef sdiff_t sdiffs_[32]
  mk_sdiff_array(sds, n_sdiffs, &sdiffs_[0])
  cdef np.ndarray[np.uint16_t, ndim=1, mode="c"] lock_counters_ \
    = np.array(lock_counters, dtype=np.uint16)
  cdef gnss_signal_t sats_to_drop[32]
  cdef u16 *locks = &lock_counters_[0]
  cdef s8 dropped = check_lock_counters(n_sdiffs, &sdiffs_[0], &locks, &sats_to_drop[0])
  return (lock_counters_, read_signal_array(dropped, sats_to_drop))

def copy_sdiffs_put_ref_first_(GNSSSignal ref_sid, sdiffs):
  """Wraps copy_sdiffs_put_ref_first

  Parameters
  ----------

  Returns
  -------
  """
  num_sdiffs = len(sdiffs)
  cdef sdiff_t sdiffs_[32], sdiffs_out[32]
  mk_sdiff_array(sdiffs, num_sdiffs, &sdiffs_[0])
  cdef s8 code = copy_sdiffs_put_ref_first(ref_sid._thisptr, num_sdiffs,&sdiffs_[0], &sdiffs_out[0])
  return (code, read_sdiff_array(num_sdiffs, &sdiffs_out[0]))

def filter_sdiffs_(sdiffs, num_sats_to_drop, sats_to_drop):
  """Wraps filter_sdiffs.

  """
  raise NotImplementedError

def debug_sdiff_(SingleDiff sd):
  """Prints an SingleDiff.
  """
  debug_sdiff(sd._thisptr)

def debug_sdiffs_(sd):
  n = len(sd)
  cdef sdiff_t sd_[32]
  mk_sdiff_array(sd, n, &sd_[0])
  debug_sdiffs(n, &sd_[0])

# TODO (Buro): Look into Numpy typed memory views for dealing with
# this more idiomatically.

cdef mk_sdiff_array(py_sdiffs, u8 n_c_sdiffs, sdiff_t *c_sdiffs):
  """Given an array of python SingleDiffs, copies their contents to an
  array of sdiff_t's.

  """
  if n_c_sdiffs < len(py_sdiffs):
    raise ValueError("The length of the c sdiffs array (" + str(n_c_sdiffs) + \
                      ") must be at least the length of the python sdiffs array " + \
                      str(len(py_sdiffs)) + ").")
  cdef sdiff_t sd_
  for (i, sdiff) in enumerate(py_sdiffs):
    sd_ = (<SingleDiff> sdiff)._thisptr
    memcpy(&c_sdiffs[i], &sd_, sizeof(sdiff_t))

cdef read_sdiff_array(u8 n_c_sdiffs, sdiff_t *c_sdiffs):
  """Given an array of c sdiff_t's, returns an array of SingleDiffs.

  """
  cdef sdiff_t sd_
  py_sdiffs = [SingleDiff()]*n_c_sdiffs
  for (i, sdiff) in enumerate(py_sdiffs):
    sd_ = (<SingleDiff> sdiff)._thisptr
    memcpy(&sd_, &c_sdiffs[i], sizeof(sdiff_t))
  return py_sdiffs
