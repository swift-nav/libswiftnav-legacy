# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from libcpp cimport bool
from signal cimport gnss_signal_t

cdef extern from "libswiftnav/bit_sync.h":
  enum:
    BITSYNC_UNSYNCED
    BIT_LENGTH_MAX

  ctypedef struct bit_sync_t:
    u8 bit_length
    u8 bit_phase
    s8 bit_phase_ref
    s32 bit_integrate
    u8 bitsync_count
    s32 bitsync_prev_corr[BIT_LENGTH_MAX]
    u32 bitsync_histogram[BIT_LENGTH_MAX]

  void bit_sync_init(bit_sync_t *b, gnss_signal_t sid)
  bool bit_sync_update(bit_sync_t *b, s32 corr_prompt_real, u32 ms,
                       s32 *bit_integrate)

cdef class BitSync:
  cdef bit_sync_t _thisptr
