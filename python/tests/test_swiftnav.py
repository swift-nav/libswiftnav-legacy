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

def test_imports():
  """Verify that distributed packages survive setuptools installation.

  """
  import swiftnav.almanac
  import swiftnav.amb_kf
  import swiftnav.ambiguity_test
  import swiftnav.baseline
  import swiftnav.bits
  import swiftnav.constants
  import swiftnav.coord_system
  import swiftnav.correlate
  import swiftnav.dgnss_management
  import swiftnav.edc
  import swiftnav.ephemeris
  import swiftnav.gpstime
  import swiftnav.lambda_py
  import swiftnav.linear_algebra
  import swiftnav.logging
  import swiftnav.memory_pool
  import swiftnav.nav_msg
  import swiftnav.printing_utils
  import swiftnav.prns
  import swiftnav.pvt
  import swiftnav.rtcm3
  import swiftnav.sats_management_py
  import swiftnav.single_diff
  import swiftnav.track
  import swiftnav.tropo
