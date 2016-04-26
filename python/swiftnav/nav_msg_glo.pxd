# Copyright (C) 2016 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

# cython: embedsignature=True

from common cimport *
from ephemeris cimport *
from libcpp cimport bool
from signal cimport gnss_signal_t

cdef extern from "libswiftnav/nav_msg_glo.h":
  enum:
    NAV_MSG_GLO_STRING_BITS_LEN
    GLO_TM
    GLO_TM_LEN

  ctypedef enum glo_receive_machine:
    INVALID
    SYNC_TM
    GET_DATA_BIT

  ctypedef struct nav_msg_glo_t:
    u32 string_bits[NAV_MSG_GLO_STRING_BITS_LEN]
    u16 current_head_bit_index
    u16 nt
    u8 next_string_id
    u8 n4
    u8 hrs
    u8 min
    u8 sec
    glo_receive_machine state
    u8 meander_bits_cnt
    u8 manchester
    u8 decode_done

  void nav_msg_init_glo(nav_msg_glo_t *n)
  s8 process_string_glo(nav_msg_glo_t *n, ephemeris_t *e)
  s8 nav_msg_update_glo(nav_msg_glo_t *n, bool bit_val)
  double nav_msg_get_tow_glo(nav_msg_glo_t *n)
  s8 error_detection_glo(nav_msg_glo_t *n)

cdef class NavMsgGlo:
  cdef nav_msg_glo_t _thisptr
