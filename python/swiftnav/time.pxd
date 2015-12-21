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
  u32 GPS_EPOCH
  u32 MJD_JAN_6_1980
  u32 MJD_JAN_1_1601
  double WN_UNKNOWN
  u32 WEEK_SECS
  u32 DAY_SECS
  u32 HOUR_SECS
  u32 MINUTE_SECS
  u32 COMMON_YEAR_DAYS
  u32 LEAP_YEAR_DAYS
  u32 FOUR_YEARS_DAYS
  u32 HUNDRED_YEARS_DAYS
  u32 FOUR_HUNDRED_YEARS_DAYS

  ctypedef struct gps_time_t:
    double tow
    s16 wn

  void normalize_gps_time(const gps_time_t *t)
  time_t gps2time(const gps_time_t *t_gps)
  double gpsdifftime(const gps_time_t *end, const gps_time_t *beginning)
  void gps_time_match_weeks(gps_time_t *t, const gps_time_t *ref)

cdef class GpsTime:
  cdef gps_time_t _thisptr
