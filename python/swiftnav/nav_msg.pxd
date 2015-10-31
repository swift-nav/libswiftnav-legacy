# Copyright (C) 2012 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

cimport nav_msg_c
cimport ephemeris_c
from cpython cimport bool

cdef class NavMsg:
  cdef nav_msg_c.nav_msg_t state
  cdef ephemeris_c.ephemeris_t eph
  cdef readonly bool eph_valid
  cdef readonly bit_phase
  cdef readonly bit_phase_ref
