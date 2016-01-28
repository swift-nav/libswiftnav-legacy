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

import swiftnav.ambiguity_test as a
import swiftnav.observation as o
import swiftnav.sats_management as sat
import swiftnav.signal as s

def test_update_sats_same_sats():
  sids = [s.GNSSSignal(sat=3, code=0),
          s.GNSSSignal(sat=1, code=0),
          s.GNSSSignal(sat=2, code=0),
          s.GNSSSignal(sat=4, code=0)]
  sm = sat.SatsManagement(sids=sids)
  test = a.AmbiguityTest(sats=sm)
  sdiffs = [o.SingleDiff(sid={'sat':1, 'code': 0},
                         pseudorange=0, sat_pos=(0, 0, 0), sat_vel=(0, 0, 0),
                         carrier_phase=0, raw_doppler=0, doppler=0, snr=0),
            o.SingleDiff(sid={'sat':2, 'code': 0},
                         pseudorange=0, sat_pos=(0, 0, 0), sat_vel=(0, 0, 0),
                         carrier_phase=0, raw_doppler=0, doppler=0, snr=0),
            o.SingleDiff(sid={'sat':3, 'code': 0},
                         pseudorange=0, sat_pos=(0, 0, 0), sat_vel=(0, 0, 0),
                         carrier_phase=0, raw_doppler=0, doppler=0, snr=0),
            o.SingleDiff(sid={'sat':4, 'code': 0},
                         pseudorange=0, sat_pos=(0, 0, 0), sat_vel=(0, 0, 0),
                         carrier_phase=0, raw_doppler=0, doppler=0, snr=0)]
  assert not test.ambiguity_update_sats(sdiffs, None, None, None, None, False)
  assert test.sats['sids'][0]['sat'] == 3
  assert test.sats['sids'][1]['sat'] == 1
  assert test.sats['sids'][2]['sat'] == 2
  assert test.sats['sids'][3]['sat'] == 4

def test_bad_measurements():
  sids = [s.GNSSSignal(sat=1, code=0),
          s.GNSSSignal(sat=2, code=0),
          s.GNSSSignal(sat=3, code=0),
          s.GNSSSignal(sat=5, code=0),
          s.GNSSSignal(sat=6, code=0)]
  float_sats = sat.SatsManagement(sids=sids)
  sids = [s.GNSSSignal(sat=3, code=0),
          s.GNSSSignal(sat=1, code=0),
          s.GNSSSignal(sat=2, code=0),
          s.GNSSSignal(sat=5, code=0),
          s.GNSSSignal(sat=6, code=0)]
  amb_sats_init = sat.SatsManagement(sids=sids)
  test = a.AmbiguityTest(sats=float_sats)
  sdiffs = [o.SingleDiff(sid={'sat':1, 'code': 0},
                         pseudorange=0, sat_pos=(0, 0, 0), sat_vel=(0, 0, 0),
                         carrier_phase=0, raw_doppler=0, doppler=0, snr=0),
            o.SingleDiff(sid={'sat':2, 'code': 0},
                         pseudorange=0, sat_pos=(0, 0, 0), sat_vel=(0, 0, 0),
                         carrier_phase=0, raw_doppler=0, doppler=0, snr=0),
            o.SingleDiff(sid={'sat':3, 'code': 0},
                         pseudorange=0, sat_pos=(0, 0, 0), sat_vel=(0, 0, 0),
                         carrier_phase=0, raw_doppler=0, doppler=0, snr=0),
            o.SingleDiff(sid={'sat':4, 'code': 0},
                         pseudorange=0, sat_pos=(0, 0, 0), sat_vel=(0, 0, 0),
                         carrier_phase=0, raw_doppler=0, doppler=0, snr=0),
            o.SingleDiff(sid={'sat':5, 'code': 0},
                         pseudorange=0, sat_pos=(0, 0, 0), sat_vel=(0, 0, 0),
                         carrier_phase=0, raw_doppler=0, doppler=0, snr=0),
            o.SingleDiff(sid={'sat':6, 'code': 0},
                         pseudorange=0, sat_pos=(0, 0, 0), sat_vel=(0, 0, 0),
                         carrier_phase=0, raw_doppler=0, doppler=0, snr=0)]
  #ambiguity_update_sats(&amb_test, num_sdiffs, sdiffs, &float_sats, est, U, D, false);
