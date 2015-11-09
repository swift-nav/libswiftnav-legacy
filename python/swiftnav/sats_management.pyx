# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

cdef class SatsManagement:

  def __cinit__(self, sdiffs, sdiffs_with_ref_first):
    raise NotImplementedError
    #init_sats_management(&self._thisptr, sdiffs, sdiffs_with_ref_first)

  def _print(self):
    print_sats_management(&self._thisptr)

  def _print_short(self):
    print_sats_management_short(&self._thisptr)

  def rebase(self, sdiffs, sdiffs_with_ref_first):
    raise NotImplementedError
    #return rebase_sats_management(&self._thisptr, num_sats, sdiffs, sdiffs_with_ref_first)

  def update(self, non_ref_sdiffs):
    raise NotImplementedError
    #update_sats_sats_management(&self._thisptr, num_non_ref_sdiffs, sdiff_t *non_ref_sdiffs)

  def match_sdiffs_to_sats_man(self, u8 num_sdiffs, sdiffs, sdiffs_with_ref_first):
    raise NotImplementedError
    #match_sdiffs_to_sats_man(sats_management_t *sats, u8 num_sdiffs, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first):

def choose_reference_sat_(num_sats, sats):
  raise NotImplementedError

def set_reference_sat_of_prns_(u8 ref_prn, prns):
  raise NotImplementedError
