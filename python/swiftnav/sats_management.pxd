# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from constants cimport *
from observation cimport *

cdef extern from "libswiftnav/sats_management.h":
  enum:
    OLD_REF
    NEW_REF
    NEW_REF_START_OVER
    INTERSECTION_SATS_THRESHOLD_SIZE

  ctypedef struct sats_management_t:
    u8 num_sats
    u8 prns[MAX_CHANNELS]

  u8 choose_reference_sat(const u8 num_sats, const sdiff_t *sats)

  void init_sats_management(sats_management_t *sats_management, const u8 num_sats,
                            const sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
  void print_sats_management(sats_management_t *sats_management)
  void print_sats_management_short(sats_management_t *sats_management)
  s8 rebase_sats_management(sats_management_t *sats_management,
                            const u8 num_sats, const sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
  void update_sats_sats_management(sats_management_t *sats_management, u8 num_non_ref_sdiffs, sdiff_t *non_ref_sdiffs)

  void set_reference_sat_of_prns(u8 ref_prn, u8 num_sats, u8 *prns)
  s8 match_sdiffs_to_sats_man(sats_management_t *sats, u8 num_sdiffs, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)

cdef class SatsManagement:
  cdef sats_management_t _thisptr
