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
    print t.GpsTime(tow=a[0], wn=a[1]).normalize_gps_time()
    assert t.GpsTime(tow=a[0], wn=a[1]).gpsdifftime(t.GpsTime(tow=b[0], wn=b[1])) - dt < tow_tol
