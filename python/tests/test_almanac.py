#!/usr/bin/env python
# Copyright (C) 2016 Swift Navigation Inc.
# Contact: Bhaskar Mookerji <mookerji@swiftnav.com>
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

import numpy as np
import swiftnav.almanac

def test_init():
  alm = {
    'gps': {
        'a': 8,
        'af0': 9,
        'af1': 10,
        'argp': 7,
        'ecc': 8,
        'inc': 3,
        'ma': 8,
        'raaw': 6,
        'rora': 4,
        'toa': 2,
        'week': 11
    },
    'healthy': 1,
    'sid': {
        'code': 0,
        'sat': 1
    },
    'valid': 1,
  }

  satAlmanac = swiftnav.almanac.Almanac(**alm)
  assert np.isclose(alm['gps']['a'], satAlmanac.gps['a'])
  assert np.isclose(alm['gps']['ecc'], satAlmanac.gps['ecc'])

def test_almanac_functions():
  alm = {
    'healthy': 1,
    'gps': {
        'a': 26559810.38052176,
        'week': 814,
        'ecc': 0.004033565521,
        'argp': 0.380143734,
        'af0': -7.629394531e-06,
        'rora': -7.874613724e-09,
        'ma': 2.030394348,
        'toa': 233472.0,
        'inc': 0.9619461694,
        'af1': 0.0,
        'raaw': -0.9244574916
    },
    'valid': 1,
    'sid': {
        'code': 0,
        'sat': 1
    }
  }
  tow = 168214.6
  wn = 814
  pos = np.array([-2704369.61784456, -4263211.09418205,  3884641.21270987])

  satAlmanac = swiftnav.almanac.Almanac(**alm)
  position, velocity = satAlmanac.calc_state(tow, week=wn)
  assert np.allclose(position, np.array([7938925.33852649, -19535003.95771277, -16086183.49111858]))
  assert np.allclose(velocity, np.array([1771.68410121, -1031.78569083,  2151.49172217]))

  doppler = satAlmanac.calc_doppler(tow, pos, week=wn)
  assert np.isclose(doppler, -1607.88747151)

  azimuth, elevation = satAlmanac.calc_az_el(tow, pos, week=wn)
  assert np.isclose(azimuth, 2.43752431894)
  assert np.isclose(elevation, -0.239507449606)
