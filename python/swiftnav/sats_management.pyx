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
from libc.string cimport memset, memcpy
from observation cimport *
from signal cimport *

cdef class SatsManagement:

  def __init__(self, sids):
    memset(&self._thisptr, 0, sizeof(sats_management_t))
    num_sids = len(sids)
    cdef gnss_signal_t sids_[11]
    mk_signal_array(sids, num_sids, &sids_[0])
    self._thisptr.num_sats = num_sids
    self._thisptr.sids = sids_

  def __getattr__(self, k):
    return self._thisptr.get(k)

  def __repr__(self):
    return fmt_repr(self)

  def to_dict(self):
    return self._thisptr

  def from_dict(self, d):
    self._thisptr = d

  def init(self, sdiffs):
    """Init sats management?

    Parameters
    ----------

    Returns
    -------

    """
    n_sdiffs = len(sdiffs)
    cdef sdiff_t sdiffs_[32], sdiffs_with_ref_first[32]
    mk_sdiff_array(sdiffs, n_sdiffs, &sdiffs_[0])
    init_sats_management(&self._thisptr, n_sdiffs, sdiffs_, sdiffs_with_ref_first)
    return read_sdiff_array(n_sdiffs, &sdiffs_with_ref_first[0])

  def _print(self):
    print_sats_management(&self._thisptr)

  def _print_short(self):
    print_sats_management_short(&self._thisptr)

  def rebase(self, sdiffs):
    """Updates sats to the new measurements' sat set.

    Parameters
    ----------

    Returns
    -------

    """
    num_sats = len(sdiffs)
    cdef sdiff_t sdiffs_[32], sdiffs_with_ref_first[32]
    mk_sdiff_array(sdiffs, num_sats, &sdiffs_[0])
    cdef u8 code = rebase_sats_management(&self._thisptr, num_sats,
                                          &sdiffs_[0], &sdiffs_with_ref_first[0])
    return (code, read_sdiff_array(num_sats, &sdiffs_with_ref_first[0]))

  def update(self, non_ref_sdiffs):
    """Updates sats to the new measurements' sat set.

    Parameters
    ----------

    Returns
    -------

    """
    num_non_ref_sdiffs = len(non_ref_sdiffs)
    cdef sdiff_t non_ref_sdiffs_[32]
    mk_sdiff_array(non_ref_sdiffs, num_non_ref_sdiffs, &non_ref_sdiffs_[0])
    update_sats_sats_management(&self._thisptr, num_non_ref_sdiffs, &non_ref_sdiffs_[0])

  def match_sdiffs_to_sats_man(self, sdiffs):
    """

    Parameters
    ----------

    Returns
    -------

    """
    num_sdiffs = len(sdiffs)
    cdef sdiff_t sdiffs_[32], sdiffs_with_ref_first[32]
    mk_sdiff_array(sdiffs, num_sdiffs, &sdiffs_[0])
    cdef s8 code \
      = match_sdiffs_to_sats_man(&self._thisptr, num_sdiffs, &sdiffs_[0], &sdiffs_with_ref_first[0])
    return (code, read_sdiff_array(num_sdiffs, &sdiffs_with_ref_first[0]))

def choose_reference_sat_(sats):
  """

  Parameters
  ----------

  Returns
  -------

  """
  num_sats = len(sats)
  cdef sdiff_t sdiffs[32]
  mk_sdiff_array(sats, num_sats, &sdiffs[0])
  best_sid = GNSSSignal()
  best_sid._thisptr = choose_reference_sat(num_sats, &sdiffs[0])
  return best_sid

def set_reference_sat_of_sids_(GNSSSignal ref_sid, sids):
  """

  Parameters
  ----------

  Returns
  -------

  """
  num_sats = len(sids)
  cdef gnss_signal_t sids_[32]
  mk_signal_array(sids, num_sats, &sids_[0])
  set_reference_sat_of_sids(ref_sid._thisptr, num_sats, &sids_[0])
  return read_signal_array(num_sats, &sids_[0])
