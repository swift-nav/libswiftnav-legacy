# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from libc.stdint cimport *

cdef extern from "libswiftnav/common.h":
  double ABS(double)
  double MIN(double, double)
  double MAX(double, double)
  double CLAMP_DIFF(double, double)

  ctypedef int8_t s8
  ctypedef int16_t s16
  ctypedef int32_t s32
  ctypedef int64_t s64
  ctypedef uint8_t u8
  ctypedef uint16_t u16
  ctypedef uint32_t u32
  ctypedef uint64_t u64
