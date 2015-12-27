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

import numpy as np
import swiftnav.almanac

def test_imports():
  """Verify that distributed packages survive setuptools installation.

  """
  assert True

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
        'band': 0,
        'constellation': 0,
        'sat': 1
    },
    'valid': 1,
  }

  satAlmanac = swiftnav.almanac.Almanac(**alm)
  assert np.isclose(alm['gps']['a'], satAlmanac.gps['a'])
  assert np.isclose(alm['gps']['ecc'], satAlmanac.gps['ecc'])

