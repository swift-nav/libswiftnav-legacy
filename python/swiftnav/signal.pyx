# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

# cython: embedsignature=True

from common cimport *
from fmt_utils import fmt_repr
from libcpp cimport bool
from libc.string cimport memset, memcpy

cdef class GNSSSignal:

  def __init__(self, **kwargs):
    memset(&self._thisptr, 0, sizeof(gnss_signal_t))
    if kwargs:
      self._thisptr = kwargs

  def __rich_cmp__(self, GNSSSignal sig, op):
    return sid_compare(self._thisptr, sig._thisptr)

  def __getattr__(self, k):
    return self._thisptr.get(k)

  def __repr__(self):
    return fmt_repr(self)

  def to_dict(self):
    return self._thisptr

  def from_dict(self, d):
    self._thisptr = d

  def is_valid(self):
    return sid_valid(self._thisptr)

  def to_index(self):
    return sid_to_index(self._thisptr)

  def to_string(self):
    cdef u8 n = 255
    cdef char s[255]
    cdef s8 length = sid_to_string(&s[0], n, self._thisptr)
    return s[:length].decode('UTF-8')

def sid_from_index(i):
  sid = GNSSSignal()
  sid._thispr = sid_from_index(i)
  return sid

cdef mk_signal_array(py_signals, u8 n_c_signals, gnss_signal_t *c_signals):
  """Given an array of python SingleDiffs, copies their contents to an
  array of gnss_signal_t's.

  """
  if n_c_signals < len(py_signals):
    raise ValueError("The length of the c signals array (" + str(n_c_signals) + \
                      ") must be at least the length of the python signals array " + \
                      str(len(py_signals)) + ").")
  cdef gnss_signal_t sd_
  for (i, signal) in enumerate(py_signals):
    sd_ = (<GNSSSignal> signal)._thisptr
    memcpy(&c_signals[i], &sd_, sizeof(gnss_signal_t))

cdef read_signal_array(u8 n_c_signals, gnss_signal_t *c_signals):
  """Given an array of c gnss_signal_t's, returns an array of GNSSSignals.

  """
  cdef gnss_signal_t sd_
  py_signals = [GNSSSignal()]*n_c_signals
  for (i, signal) in enumerate(py_signals):
    sd_ = (<GNSSSignal> signal)._thisptr
    memcpy(&sd_, &c_signals[i], sizeof(gnss_signal_t))
  return py_signals
