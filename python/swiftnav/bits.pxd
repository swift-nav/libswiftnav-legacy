# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *

cdef extern from "libswiftnav/bits.h":
  u8 parity(u32 x)
  u32 getbitu(const u8 *buff, u32 pos, u8 len)
  s32 getbits(const u8 *buff, u32 pos, u8 len)
  void setbitu(u8 *buff, u32 pos, u32 len, u32 data)
  void setbits(u8 *buff, u32 pos, u32 len, s32 data)
