# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

import numpy as np
cimport numpy as np
from libc.string cimport memcpy

# TODO (Buro): Finish this!

cdef class SingleDiff:

  def __init__(self, **kwargs):
    if kwargs:
      self._thisptr = kwargs

  def __repr__(self):
    return "<SingleDiff prn=" + str(self.prn) + ", snr=" + str(self.snr) + ">"

def single_diff_(m_a, m_b):
  raise NotImplementedError

def sdiff_search_prn_(a, b):
  raise NotImplementedError

def make_propagated_sdiffs_(m_local, m_remote, remote_dists, remote_pos_ecef, es, t, sds):
  raise NotImplementedError

def make_dd_measurements_and_sdiffs_(ref_prn, non_ref_prns, num_dds, num_sdiffs, sdiffs_in):
  raise NotImplementedError

def check_lock_counters_(n_sds, sds, lock_counters):
  raise NotImplementedError

def copy_sdiffs_put_ref_first_(ref_prn, num_sdiffs, sdiffs):
  raise NotImplementedError

def filter_sdiffs_(num_sdiffs, sdiffs, num_sats_to_drop, sats_to_drop):
  raise NotImplementedError

def debug_sdiff_(sd):
  raise NotImplementedError
