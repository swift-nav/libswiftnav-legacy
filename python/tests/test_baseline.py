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

import swiftnav.baseline as base

def test_basic():
  num_dds = 4
  N = base.carray_double(4)
  DE = base.carray_double(12)
  b = base.carray_double(3)
  dd_obs = base.carray_double(3)
  base.predict_carrier_obs(num_dds, N, DE, b, dd_obs)
  return dd_obs[0], dd_obs[1], dd_obs[2]
