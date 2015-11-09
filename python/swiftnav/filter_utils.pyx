# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

"""Utility functions for use in filters.

"""

cimport numpy as np
import numpy as np
from observation cimport *

def simple_amb_measurement_(carrier, code):
  return simple_amb_measurement(carrier, code)

def assign_de_mtx_(sats_with_ref_first, ref_ecef):
  num_sats = len(sats_with_ref_first)
  cdef sdiff_t sats_with_ref_first_[32]
  mk_sdiff_array(sats_with_ref_first, 32, &sats_with_ref_first_[0])
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(ref_ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] DE = np.empty([0,0,0]*(num_sats-1), dtype=np.double)
  cdef s8 code = assign_de_mtx(num_sats, &sats_with_ref_first_[0], &ref_ecef_[0], &DE[0,0])
  return (code, DE)
