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
import swiftnav.almanac as a
import swiftnav.time as t

def test_init():
  alm = {
    'sid': {'code': 0, 'sat': 1},
    'toa': {'tow': 2, 'wn': 11,},
    'ura': 10.0,
    'fit_interval': 13,
    'valid': 1,
    'health_bits': 0,
    'kepler': {
        'sqrta': 8,
        'af0': 9,
        'af1': 10,
        'w': 7,
        'ecc': 8,
        'inc': 3,
        'm0': 8,
        'omega0': 6,
        'omegadot': 4,
    }
  }

  satAlmanac = a.Almanac(**alm)
  assert np.isclose(alm['kepler']['sqrta'], satAlmanac.kepler['sqrta'])
  assert np.isclose(alm['kepler']['ecc'], satAlmanac.kepler['ecc'])

def test_almanac_functions():
  alm = {
    'sid': {'sat': 1, 'code': 0},
    'toa': {'tow': 233472.0, 'wn': 814,},
    'ura': 900.0,
    'fit_interval': 144 * 60 * 60,
    'valid': 1,
    'health_bits': 0,
    'kepler': {
        'sqrta': np.sqrt(26559810.38052176),
        'ecc': 0.004033565521,
        'w': 0.380143734,
        'omega0': -0.9244574916,
        'omegadot': -7.874613724e-09,
        'm0': 2.030394348,
        'inc': 0.9619461694,
        'af0': -7.629394531e-06,
        'af1': 0.0,
    }
  }
  time = t.GpsTime(**{ 'wn': 814, 'tow': 168214.0,})
  pos = np.array([-2704369.61784456, -4263211.09418205,  3884641.21270987])

  satAlmanac = a.Almanac(**alm)
  position, velocity, clock_err, clock_rate_err = satAlmanac.calc_state(time)
  assert np.allclose(position, np.array([7937862.2780455 , -19534384.87638699, -16087474.32428426]))
  assert np.allclose(velocity, np.array([1771.85083411, -1031.81872491, 2151.28549392]))
  assert np.isclose(clock_err, -7.629394531e-06)
  assert np.isclose(clock_rate_err, 0.0)

  doppler = satAlmanac.calc_doppler(time, pos)
  assert np.isclose(doppler, -1607.665252834069)

  azimuth, elevation = satAlmanac.calc_az_el(time, pos)
  assert np.isclose(azimuth, 2.437585312667075)
  assert np.isclose(elevation, -0.23953409326103278)
