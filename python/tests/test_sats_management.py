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

import swiftnav.observation as o
import swiftnav.sats_management as sat
import swiftnav.signal as s

def test_creation():
  new_ref = s.GNSSSignal(sat=8, band=0, constellation=0)
  sids = [s.GNSSSignal(sat=2, band=0, constellation=0),
          s.GNSSSignal(sat=1, band=0, constellation=0),
          s.GNSSSignal(sat=3, band=0, constellation=0),
          s.GNSSSignal(sat=4, band=0, constellation=0)]
  sm = sat.SatsManagement(sids=sids)
  assert sm.num_sats
  assert sm.sids[0]['sat'] == 2
  assert sm.sids[1]['sat'] == 1
  assert sm.sids[2]['sat'] == 3
  assert sm.sids[3]['sat'] == 4
  assert not sm._print()
  assert not sm._print_short()
  sdiffs = [o.SingleDiff(sid={'sat':1, 'band': 0, 'constellation': 0},
                         pseudorange=0, sat_pos=(0, 0, 0), sat_vel=(0, 0, 0),
                         carrier_phase=0, raw_doppler=0, doppler=0, snr=1),
            o.SingleDiff(sid={'sat':2, 'band': 0, 'constellation': 0},
                         pseudorange=0, sat_pos=(0, 0, 0), sat_vel=(0, 0, 0),
                         carrier_phase=0, raw_doppler=0, doppler=0, snr=5),
            o.SingleDiff(sid={'sat':3, 'band': 0, 'constellation': 0},
                         pseudorange=0, sat_pos=(0, 0, 0), sat_vel=(0, 0, 0),
                         carrier_phase=0, raw_doppler=0, doppler=0, snr=10)]
  # TODO (Buro): Check outputs!
  assert sm.rebase(sdiffs)
  assert not sm.update(sdiffs)

def test_choose_ref_sat():
  sdiffs = [o.SingleDiff(sid={'sat':1, 'band': 0, 'constellation': 0},
                         pseudorange=0, sat_pos=(0, 0, 0), sat_vel=(0, 0, 0),
                         carrier_phase=0, raw_doppler=0, doppler=0, snr=0),
            o.SingleDiff(sid={'sat':2, 'band': 0, 'constellation': 0},
                         pseudorange=0, sat_pos=(0, 0, 0), sat_vel=(0, 0, 0),
                         carrier_phase=0, raw_doppler=0, doppler=0, snr=0),
            o.SingleDiff(sid={'sat':3, 'band': 0, 'constellation': 0},
                         pseudorange=0, sat_pos=(0, 0, 0), sat_vel=(0, 0, 0),
                         carrier_phase=0, raw_doppler=0, doppler=0, snr=0)]
  assert isinstance(sat.choose_reference_sat_(sdiffs), s.GNSSSignal)
  assert sat.choose_reference_sat_(sdiffs).sat == 1
  sdiffs = [o.SingleDiff(sid={'sat':1, 'band': 0, 'constellation': 0},
                         pseudorange=0, sat_pos=(0, 0, 0), sat_vel=(0, 0, 0),
                         carrier_phase=0, raw_doppler=0, doppler=0, snr=1),
            o.SingleDiff(sid={'sat':2, 'band': 0, 'constellation': 0},
                         pseudorange=0, sat_pos=(0, 0, 0), sat_vel=(0, 0, 0),
                         carrier_phase=0, raw_doppler=0, doppler=0, snr=5),
            o.SingleDiff(sid={'sat':3, 'band': 0, 'constellation': 0},
                         pseudorange=0, sat_pos=(0, 0, 0), sat_vel=(0, 0, 0),
                         carrier_phase=0, raw_doppler=0, doppler=0, snr=10)]
  assert isinstance(sat.choose_reference_sat_(sdiffs), s.GNSSSignal)
  assert sat.choose_reference_sat_(sdiffs).sat == 3
