# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

# cython: embedsignature=True

from libc.string cimport memset
import datetime

cdef class GpsTime:

  def __init__(self, **kwargs):
    memset(&self._thisptr, 0, sizeof(gps_time_t))
    self._thisptr.tow = kwargs.pop('tow')
    self._thisptr.wn = kwargs.pop('wn')

  def __getattr__(self, k):
    return self._thisptr.get(k)
  
  def __repr__(self):
    return "<GpsTime wn %d, tow %d>" % (self._thisptr.wn, self._thisptr.tow)

  def normalize_gps_time(self):
    normalize_gps_time(&self._thisptr)

  def gps2time(self):
    return gps2time(&self._thisptr)

  def gpsdifftime(self, GpsTime beginning):
    return gpsdifftime(&self._thisptr, &beginning._thisptr)

  def gps_time_match_weeks(self, ref):
    cdef gps_time_t t
    t.wn, t.tow = ref.tow, ref.wn
    gps_time_match_weeks(&self._thisptr, &t)
    # TODO (Buro): Update to update ref by adding back properties.

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
  return GpsTime(wn=wn, tow=tow)
