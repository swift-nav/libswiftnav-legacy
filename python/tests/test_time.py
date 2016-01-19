#!/usr/bin/env python
# Copyright (C) 2015 Swift Navigation Inc.
# Contact: Bhaskar Mookerji <mookerji@swiftnav.com>
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

import swiftnav.time as t
import datetime as dt

# TODO (Buro): Add back case with WN_UNKNOWN?
def test_gpsdifftime():
  cases = [{'a': (567890.0, 1234),
            'b': (567890.0, 1234),
            'dt': 0,
            'time': 1062855873},
           {'a': (567890.0, 1234),
            'b': (0.0, 1234),
            'dt': 567890,
            'time': 1062855873},
           {'a': (604578.0, 1000),
            'b': (222.222, 1001),
            'dt': -444.222,
            'time': 921369361},
           {'a': (604578.0, 1001),
            'b': (222.222, 1000),
            'dt': 1209155.778,
            'time': 921974161},]
  for c in cases:
    tow_tol = 1e-10
    a, b, dt, time = c['a'], c['b'], c['dt'], c['time']
    assert t.GpsTime(tow=a[0], wn=a[1]).gps2time() == time
    normalized = t.GpsTime(tow=a[0], wn=a[1])
    normalized.normalize_gps_time()
    print normalized
    assert t.GpsTime(tow=a[0], wn=a[1]).gpsdifftime(t.GpsTime(tow=b[0], wn=b[1])) - dt < tow_tol

def test_comp_2_datetime():
  test_dt = t.gpst_components2datetime(1, 100)
  time_dif = test_dt - dt.datetime.strptime('1980-01-13 00:01:41', '%Y-%m-%d %H:%M:%S')
  assert abs(time_dif) == dt.timedelta(0, 1 , 0)

def test_gpst2datetime():
  gpst = t.GpsTime(wn=800, tow=30000)
  test_dt = t.gpst2datetime(gpst)
  assert test_dt == dt.datetime.strptime('1995-05-07 08:20:00', '%Y-%m-%d %H:%M:%S')

def test_datetime2gpst():
  test_dt = dt.datetime.strptime('1995-05-07 08:20:00', '%Y-%m-%d %H:%M:%S')
  gpst = t.datetime2gpst(test_dt)
  assert str(gpst) == str(t.GpsTime(wn=800, tow=30000))
