# Copyright (C) 2016 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from libcpp cimport bool

cdef extern from "libswiftnav/cnav_msg.h":

  ctypedef struct cnav_msg_decoder_t:
    pass

  ctypedef struct cnav_msg_t:
    u8   prn    # SV PRN. 0..31
    u8   msg_id # Message id. 0..31
    u32  tow    # GPS ToW in 6-second units. Multiply to 6 to get seconds. 
    bool alert  # CNAV message alert flag

  void cnav_msg_decoder_init(cnav_msg_decoder_t *dec)
  bool cnav_msg_decoder_add_symbol(cnav_msg_decoder_t *dec,
                                   u8 symbol,
                                   cnav_msg_t *msg,
                                   u32 *delay)

cdef extern from "libl2cbitstream/l2cbitstream.h":
  bool get_l2c_message(u8 *au_message, u8 prn, u8 msg_id, u32 tow)

cdef class CNavMsgDecoder:
  cdef cnav_msg_decoder_t _thisptr
cdef class CNavMsg:
  cdef cnav_msg_t _thisptr
