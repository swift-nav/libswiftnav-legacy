# Copyright (C) 2014 Swift Navigation Inc.
# Contact: Fergus Noble <fergus@swiftnav.com>
#          Bhaskar Mookerji <mookerji@swiftnav.com>
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from posix.types cimport time_t

cdef extern from "libswiftnav/time.h":
  u8 GPS_WEEK_CYCLE
  u8 GPS_MINUS_UTC_SECS
  u16 GPS_EPOCH
  double WN_UNKNOWN
  u16 WEEK_SECS
  u16 DAYS_SECS

  ctypedef struct gps_time_t:
    double tow
    s16 wn

  void normalize_gps_time(const gps_time_t *t)
  time_t gps2time(const gps_time_t *t_gps)
  double gpsdifftime(const gps_time_t *end, const gps_time_t *beginning)
  void gps_time_match_weeks(gps_time_t *t, const gps_time_t *ref)

cdef class GpsTime:
  cdef gps_time_t _thisptr
