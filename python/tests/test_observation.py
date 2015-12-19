#!/usr/bin/env pythons
# Copyright (C) 2015 Swift Navigation Inc.
# Contact: Bhaskar Mookerji <mookerji@swiftnav.com>
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from swiftnav.pvt import calc_PVT_
import numpy as np
import pytest
import swiftnav.time as t
import swiftnav.track as t
import swiftnav.observation as o
import swiftnav.signal as s

nms = \
  [t.NavigationMeasurement(sid=s.GNSSSignal(sat=9, band=0, constellation=0),
                           pseudorange=23946993.888943646,
                           raw_pseudorange=23946993.888943646,
                           sat_pos=(-19477278.087422125, -7649508.9457812719, 16674633.163554827),
                           sat_vel=(0, 0, 0), carrier_phase=0,
                           raw_doppler=0, doppler=0, lock_counter=0,
                           snr=0, lock_time=0,
                           tot=t.GpsTime(tow=0, wn=0)),
   t.NavigationMeasurement(sid=s.GNSSSignal(sat=1, band=0, constellation=0),
                           pseudorange=22932174.156858064,
                           raw_pseudorange=22932174.156858064,
                           sat_pos=(-9680013.5408340245, -15286326.354385279, 19429449.383770257),
                           sat_vel=(0, 0, 0), carrier_phase=0,
                           raw_doppler=0, doppler=0, lock_counter=0,
                           snr=0, lock_time=0,
                           tot=t.GpsTime(tow=0, wn=0)),
   t.NavigationMeasurement(sid=s.GNSSSignal(sat=2, band=0, constellation=0),
                           pseudorange=24373231.648055989,
                           raw_pseudorange=24373231.648055989,
                           sat_pos=(-19858593.085281931, -3109845.8288993631, 17180320.439503901),
                           sat_vel=(0, 0, 0), carrier_phase=0,
                           raw_doppler=0, doppler=0, lock_counter=0,
                           snr=0, lock_time=0,
                           tot=t.GpsTime(tow=0, wn=0)),
   t.NavigationMeasurement(sid=s.GNSSSignal(sat=3, band=0, constellation=0),
                           pseudorange=24779663.252316438,
                           raw_pseudorange=24779663.252316438,
                           sat_pos=(6682497.8716542246, -14006962.389166718, 21410456.275678463),
                           sat_vel=(0, 0, 0), carrier_phase=0,
                           raw_doppler=0, doppler=0, lock_counter=0,
                           snr=0, lock_time=0,
                           tot=t.GpsTime(tow=0, wn=0)),
   t.NavigationMeasurement(sid=s.GNSSSignal(sat=4, band=0, constellation=0),
                           pseudorange=26948717.022331879,
                           raw_pseudorange=26948717.022331879,
                           sat_pos=(7415370.9916331079, -24974079.044485383, -3836019.0262199985),
                           sat_vel=(0, 0, 0), carrier_phase=0,
                           raw_doppler=0, doppler=0, lock_counter=0,
                           snr=0, lock_time=0,
                           tot=t.GpsTime(tow=0, wn=0)),
   t.NavigationMeasurement(sid=s.GNSSSignal(sat=5, band=0, constellation=0),
                           pseudorange=23327405.435463827,
                           raw_pseudorange=23327405.435463827,
                           sat_pos=(-2833466.1648670658, -22755197.793894723, 13160322.082875408),
                           sat_vel=(0, 0, 0), carrier_phase=0,
                           raw_doppler=0, doppler=0, lock_counter=0,
                           snr=0, lock_time=0,
                           tot=t.GpsTime(tow=0, wn=0)),
   t.NavigationMeasurement(sid=s.GNSSSignal(sat=6, band=0, constellation=0),
                           pseudorange=27371419.016328193,
                           raw_pseudorange=27371419.016328193,
                           sat_pos=(14881660.383624561, -5825253.4316490609, 21204679.68313824),
                           sat_vel=(0, 0, 0), carrier_phase=0,
                           raw_doppler=0, doppler=0, lock_counter=0,
                           snr=0, lock_time=0,
                           tot=t.GpsTime(tow=0, wn=0)),
   t.NavigationMeasurement(sid=s.GNSSSignal(sat=7, band=0, constellation=0),
                           pseudorange=26294221.697782904,
                           raw_pseudorange=26294221.697782904,
                           sat_pos=(12246530.477279386, -22184711.955107089, 7739084.2855069181),
                           sat_vel=(0, 0, 0), carrier_phase=0,
                           raw_doppler=0, doppler=0, lock_counter=0,
                           snr=0, lock_time=0,
                           tot=t.GpsTime(tow=0, wn=0)),
   t.NavigationMeasurement(sid=s.GNSSSignal(sat=8, band=0, constellation=0),
                           pseudorange=25781999.479948733,
                           raw_pseudorange=25781999.479948733,
                           sat_pos=(-25360766.249484103, -1659033.490658124, 7821492.0398916304),
                           sat_vel=(0, 0, 0), carrier_phase=0,
                           raw_doppler=0, doppler=0, lock_counter=0,
                           snr=0, lock_time=0,
                           tot=t.GpsTime(tow=0, wn=0))]

TOL = 1e-10

def test_pvt_failed_repair():
  n_used = 5
  ret, soln, dops = calc_PVT_(nms[0:8], False)
  assert ret == 1
  assert np.allclose(soln.pos_llh, [0.6590905981580969, -2.136087375092587, 37.806255077127076])
  assert soln.valid
  assert soln.n_used == 7
  assert soln.time['tow'] - 0.07107547315305934 <= TOL
  assert soln.time['wn'] == 0
  assert soln.clock_offset - 0.008803100339752288 <= TOL
  assert soln.clock_bias - 0 <= TOL
  assert dops.pdop - 2.04639475526 <= TOL
  assert dops.gdop - 2.36782174999 <= TOL
  assert dops.tdop - 1.19115420721 <= TOL
  assert dops.hdop - 1.32894965567 <= TOL
  assert dops.vdop - 1.5561569031 <= TOL

def test_pvt_repair():
  ret, soln, dops = calc_PVT_(nms[0:9], False)
  assert ret == 1, "Return code should be 1 (pvt repair)"
  assert soln.n_used == 8
  assert np.allclose(soln.pos_llh, [0.6590899959250786, -2.1360890086733595, 28.617422790579482])
  assert soln.valid

def test_pvt_disable_pvt_raim():
  ret, soln, dops = calc_PVT_(nms[0:9], True)
  assert soln.n_used == 9
  assert ret == 2, "Return code should be 2 (raim not used)"
  assert soln.valid == 1, "Solution should be valid!"

@pytest.mark.skipif(True, reason="Don't know what's going on here.")
def test_single_diff():
  # Test for when they are interleaved
  assert len(o.single_diff_([nms[0], nms[1]], [nms[2], nms[3]])) == 0
  assert len(o.single_diff_([nms[2], nms[3]], [nms[0], nms[1]])) == 0
  # Test construction from construction
  sdiffs = o.single_diff_([nms[1], nms[0]], [nms[1], nms[0]])
  assert len(sdiffs) == 2
  sdiffs = o.single_diff_([nms[1], nms[0]], [nms[2], nms[0]])
  assert len(sdiffs) == 1
