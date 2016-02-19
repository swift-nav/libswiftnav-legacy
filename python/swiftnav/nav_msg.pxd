# Copyright (C) 2012 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from ephemeris cimport *
from libcpp cimport bool
from signal cimport gnss_signal_t

cdef extern from "libswiftnav/nav_msg.h":
  enum:
    NAV_MSG_SUBFRAME_BITS_LEN
    TOW_INVALID
    BITSYNC_UNSYNCED
    BIT_POLARITY_NORMAL
    BIT_POLARITY_INVERTED
    BIT_POLARITY_UNKNOWN

  ctypedef struct nav_msg_t:
    u32 subframe_bits[NAV_MSG_SUBFRAME_BITS_LEN]
    u16 subframe_bit_index
    bool overrun
    s16 subframe_start_index
    u32 frame_words[3][8]
    u8 next_subframe_id
    s8 bit_polarity

  ctypedef struct gps_l1ca_decoded_data_t:
    ephemeris_t ephemeris
    bool ephemeris_upd_flag
    u32 gps_l2c_sv_capability
    bool gps_l2c_sv_capability_upd_flag

  void nav_msg_init(nav_msg_t *n)
  s32 nav_msg_update(nav_msg_t *n, bool bit_val)
  bool subframe_ready(nav_msg_t *n)
  s8 process_subframe(nav_msg_t *n, gnss_signal_t sid,
                      gps_l1ca_decoded_data_t *data)

cdef class NavMsg:
  cdef nav_msg_t _thisptr

