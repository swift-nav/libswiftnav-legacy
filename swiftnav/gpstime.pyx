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
  #cdef gpstime_c.gps_time_t gps_time
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
  dt = datetime.timedelta(weeks=wn, seconds=tow)
  return dt + datetime.datetime(1980, 1, 6, 0, 0, 0) + \
      datetime.timedelta(seconds=16)

def gpst2datetime(GpsTime gpst):
  return gpst_components2datetime(gpst.wn, gpst.tow)
