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

import swiftnav.nav_msg as n

def test_basic():
  msg = n.nav_msg_t()
  assert msg.bit_phase == 0
  msg.bit_phase = 1
  msg.subframe_start_index = 1
  x = n.carray_u32(n.NAV_MSG_SUBFRAME_BITS_LEN)
  msg.subframe_bits = n.carray_u32(n.NAV_MSG_SUBFRAME_BITS_LEN)
  assert msg.bit_phase == 1
  #msg.hist = [1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]
  msg.hist = n.carray_float(20)
