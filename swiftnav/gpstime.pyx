# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

cimport gpstime_c
import datetime

cdef class GpsTime:
  # cdef gpstime_c.gps_time_t gps_time
  def __init__(self,
				       wn,
               tow):
    self.wn = wn
    self.tow = tow

  def __repr__(self):
    return "<GpsTime wn " + str(self.wn) + ", tow " + str(self.tow) + ">"

  property wn:
    def __get__(self):
      return self.gps_time.wn
    def __set__(self, wn):
      self.gps_time.wn = wn

  property tow:
    def __get__(self):
      return self.gps_time.tow
    def __set__(self, tow):
      self.gps_time.tow = tow

def gpst_components2datetime(wn, tow):
  """
  Convert GPS time in terms of week number into GPS time in terms of the
  calendar. IMPORTANTLY, the calendar GPS time is off from UTC by a number
  of leap seconds (currently 17). The UTC to GPS time conversion can change up to twice
  per year.

  Parameters
  ----------
  wn : int
    The GPS week number.
  tow : int
    The GPS time of week.

  Returns
  -------
  datetime
    The GPS time in calendar form.
  """
  dt = datetime.timedelta(weeks=wn, seconds=tow)
  return dt + datetime.datetime(1980, 1, 6, 0, 0, 0)

def gpst2datetime(GpsTime gpst):
"""
Convert a GpsTime into datetime. Both are assumed to be GPS times.
(See note in gpst_components2datetime for explanation).

Parameters
----------
gpst : GpsTime
  GPS time struct storing week number and time of week.

Returns
-------
datetime
  GPS time stored as a datetime.
"""
  return gpst_components2datetime(gpst.wn, gpst.tow)

def datetime2gpst(timestamp):
  """
  Convert a datetime (still GPS time) into a GpsTime.
  The datetime should be GPS time, which is (currently) 17 seconds
  off from UTC. The UTC to GPS time conversion can change up to twice
  per year.

  Parameters
  ----------
  timestamp : datetime
    GPS time in terms of a datetime object.

  Returns
  -------
  GpsTime
    GPS time in terms of week number and time of week.
  """
  dt = timestamp - datetime.datetime(1980, 1, 6, 0, 0, 0)
  wn = dt.days / 7
  tow = (dt - datetime.timedelta(weeks=wn)).total_seconds()
  #wn % 1024 would be standard, but we want a bijective tranformation.
  return GpsTime(wn, tow)
