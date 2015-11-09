# Copyright (C) 2012 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from ephemeris cimport ephemeris_t
from gpstime cimport gps_time_t

cdef class NavMsg:

  def __cinit__(self):
    nav_msg_init(&self._thisptr)

  def update(self, corr_prompt_real, ms):
    return nav_msg_update(&self._thisptr, corr_prompt_real, ms)

  def subframe_ready(self):
    return subframe_ready(&self._thisptr)

  def process_subframe(self, e):
    cdef ephemeris_t tmp = e._thisptr
    return process_subframe(&self._thisptr, &tmp)
