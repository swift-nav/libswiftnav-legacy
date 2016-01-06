# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

# cython: embedsignature=True

from fmt_utils import fmt_repr
from signal cimport gnss_signal_t

cdef class BitSync:

  def __cinit__(self, sid):
    bit_sync_init(&self._thisptr, sid)

  def __getattr__(self, k):
    return self._thisptr.get(k)

  def __repr__(self):
    return fmt_repr(self)

  def to_dict(self):
    return self._thisptr

  def from_dict(self, d):
    self._thisptr = d

  def update(self, corr_prompt_real, ms):
    cdef bool result
    cdef s32 bit_integrate
    result = bit_sync_update(&self._thisptr, corr_prompt_real, ms, &bit_integrate)
    return (result, bit_integrate)
