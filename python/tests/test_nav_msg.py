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

import swiftnav.nav_msg

def test_imports():
  """Verify that distributed packages survive setuptools installation.

  """
  assert True

def build_nav_msg():
  """ Instantiate NavMsg with known attribute values. """
  nm = swiftnav.nav_msg.NavMsg()
  nm.from_dict(
    {
      'bit_polarity': -1,
      'frame_words': [[0, 0, 0, 0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0, 0, 0, 0]],
      'next_subframe_id': 1,
      'overrun': False,
      'subframe_bit_index': 0,
      'subframe_bits': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
      'subframe_start_index': 0
    }
  )
  return nm

def test_instantiate():
  nm = swiftnav.nav_msg.NavMsg()

def test_build_nav_msg():
  build_nav_msg()

def test_richcmp():
  nm_a = build_nav_msg()
  nm_b = build_nav_msg()

  # Change nm_b's fields slightly.
  nm_b_dict = nm_b.to_dict()
  nm_b_dict['bit_polarity'] = 1
  nm_b.from_dict(nm_b_dict)

  assert nm_a == nm_a
  assert not(nm_a != nm_a)
  assert nm_a != nm_b
  assert not(nm_a == nm_b)

def test_rebuild():
  nm = build_nav_msg()

  rebuild, args = nm.__reduce__()

  assert rebuild(args) == nm
