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

from swiftnav.ephemeris import Ephemeris
from swiftnav.pvt import calc_PVT_
import numpy as np
import pytest
import swiftnav.time as ti
import swiftnav.track as t
import swiftnav.observation as o
import swiftnav.signal as s

nms = \
  [t.NavigationMeasurement(sid=s.GNSSSignal(sat=9, code=0),
                           pseudorange=23946993.888943646,
                           raw_pseudorange=23946993.888943646,
                           sat_pos=(-19477278.087422125, -7649508.9457812719, 16674633.163554827),
                           sat_vel=(0, 0, 0), carrier_phase=0,
                           raw_doppler=0, doppler=0, lock_counter=0,
                           snr=0, lock_time=0,
                           tot=ti.GpsTime(tow=0, wn=0)),
   t.NavigationMeasurement(sid=s.GNSSSignal(sat=1, code=0),
                           pseudorange=22932174.156858064,
                           raw_pseudorange=22932174.156858064,
                           sat_pos=(-9680013.5408340245, -15286326.354385279, 19429449.383770257),
                           sat_vel=(0, 0, 0), carrier_phase=0,
                           raw_doppler=0, doppler=0, lock_counter=0,
                           snr=0, lock_time=0,
                           tot=ti.GpsTime(tow=0, wn=0)),
   t.NavigationMeasurement(sid=s.GNSSSignal(sat=2, code=0),
                           pseudorange=24373231.648055989,
                           raw_pseudorange=24373231.648055989,
                           sat_pos=(-19858593.085281931, -3109845.8288993631, 17180320.439503901),
                           sat_vel=(0, 0, 0), carrier_phase=0,
                           raw_doppler=0, doppler=0, lock_counter=0,
                           snr=0, lock_time=0,
                           tot=ti.GpsTime(tow=0, wn=0)),
   t.NavigationMeasurement(sid=s.GNSSSignal(sat=3, code=0),
                           pseudorange=24779663.252316438,
                           raw_pseudorange=24779663.252316438,
                           sat_pos=(6682497.8716542246, -14006962.389166718, 21410456.275678463),
                           sat_vel=(0, 0, 0), carrier_phase=0,
                           raw_doppler=0, doppler=0, lock_counter=0,
                           snr=0, lock_time=0,
                           tot=ti.GpsTime(tow=0, wn=0)),
   t.NavigationMeasurement(sid=s.GNSSSignal(sat=4, code=0),
                           pseudorange=26948717.022331879,
                           raw_pseudorange=26948717.022331879,
                           sat_pos=(7415370.9916331079, -24974079.044485383, -3836019.0262199985),
                           sat_vel=(0, 0, 0), carrier_phase=0,
                           raw_doppler=0, doppler=0, lock_counter=0,
                           snr=0, lock_time=0,
                           tot=ti.GpsTime(tow=0, wn=0)),
   t.NavigationMeasurement(sid=s.GNSSSignal(sat=5, code=0),
                           pseudorange=23327405.435463827,
                           raw_pseudorange=23327405.435463827,
                           sat_pos=(-2833466.1648670658, -22755197.793894723, 13160322.082875408),
                           sat_vel=(0, 0, 0), carrier_phase=0,
                           raw_doppler=0, doppler=0, lock_counter=0,
                           snr=0, lock_time=0,
                           tot=ti.GpsTime(tow=0, wn=0)),
   t.NavigationMeasurement(sid=s.GNSSSignal(sat=6, code=0),
                           pseudorange=27371419.016328193,
                           raw_pseudorange=27371419.016328193,
                           sat_pos=(14881660.383624561, -5825253.4316490609, 21204679.68313824),
                           sat_vel=(0, 0, 0), carrier_phase=0,
                           raw_doppler=0, doppler=0, lock_counter=0,
                           snr=0, lock_time=0,
                           tot=ti.GpsTime(tow=0, wn=0)),
   t.NavigationMeasurement(sid=s.GNSSSignal(sat=7, code=0),
                           pseudorange=26294221.697782904,
                           raw_pseudorange=26294221.697782904,
                           sat_pos=(12246530.477279386, -22184711.955107089, 7739084.2855069181),
                           sat_vel=(0, 0, 0), carrier_phase=0,
                           raw_doppler=0, doppler=0, lock_counter=0,
                           snr=0, lock_time=0,
                           tot=ti.GpsTime(tow=0, wn=0)),
   t.NavigationMeasurement(sid=s.GNSSSignal(sat=8, code=0),
                           pseudorange=25781999.479948733,
                           raw_pseudorange=25781999.479948733,
                           sat_pos=(-25360766.249484103, -1659033.490658124, 7821492.0398916304),
                           sat_vel=(0, 0, 0), carrier_phase=0,
                           raw_doppler=0, doppler=0, lock_counter=0,
                           snr=0, lock_time=0,
                           tot=ti.GpsTime(tow=0, wn=0))]

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

def test_prop_single_diff():
  rover_nm = [t.NavigationMeasurement(sat_pos=[-19277207.52067729, -8215764.2479763795, 16367744.770246204],
                                      pseudorange=21123480.27105955, doppler=0, raw_doppler=0, carrier_phase=-110993309.26669005,
                                      sat_vel=[-1025.0370403901404, -1821.9217799467374, -2091.6303199092254],
                                      lock_time=0, tot=ti.GpsTime(wn=1876, tow=167131.92954684512),
                                      raw_pseudorange=21121324.476403236,
                                      snr=30.0, sid=s.GNSSSignal(code=0, sat=0), lock_counter=0),
              t.NavigationMeasurement(sat_pos=[-16580310.849158794, 918714.1939749047, 20731444.258332774],
                                      pseudorange=22432049.84763688, doppler=0, raw_doppler=0, carrier_phase=-117882743.21601027,
                                      sat_vel=[1060.9864205977192, -2411.43509917502, 953.6270954519971],
                                      lock_time=0, tot=ti.GpsTime(wn=1876, tow=167131.9251737675),
                                      raw_pseudorange=22432340.166121125,
                                      snr=30.0, sid=s.GNSSSignal(code=0, sat=2), lock_counter=0)]
  remote_dists = np.array([ 21121393.87562408,  22432814.46819838])
  base_pos = np.array( [-2704375, -4263211,  3884637])
  es = [Ephemeris(**{'ura': 0.0, 'healthy': 1,
                     'xyz': {'acc': [-6.51925802230835e-08, 4.718410826464475e-09, 2.6226835183928943],
                             'iod': 0, 'a_gf0': 5153.648313522339, 'a_gf1': 0.5977821557062277,
                             'pos': [5.122274160385132e-09, 19.0625, 259.9375],
                             'rate': [1.0151416063308716e-06, 6.260350346565247e-06, -6.332993507385254e-08],
                             'toa': 4096},
                     'valid': 1,
                     'sid': {'code': 0, 'sat': 0},
                     'toe': {'wn': 1876, 'tow': 172800.0},
                     'kepler': {'inc_dot': 4.4966158735287064e-10, 'tgd': 5.122274160385132e-09, 'omegadot': -8.099980253833005e-09,
                                'sqrta': 5153.648313522339, 'inc': 0.9634151551139846, 'cus': 6.260350346565247e-06,
                                'omega0': 0.5977821557062277, 'cuc': 1.0151416063308716e-06, 'm0': 2.6226835183928943,
                                'toc': {'wn': 1876, 'tow': 172800.0}, 'dn': 4.718410826464475e-09, 'ecc': 0.0049016030970960855,
                                'cic': -6.332993507385254e-08, 'crs': 19.0625, 'iode': 74, 'iodc': 21845, 'cis': -6.51925802230835e-08,
                                'crc': 259.9375, 'w': 0.4885959643259506, 'af0': 7.212162017822266e-06, 'af1': 9.094947017729282e-13, 'af2': 0.0},
                     'fit_interval': 4}),
        Ephemeris(**{'ura': 0.0, 'healthy': 1,
                     'xyz': {'acc': [-6.146728992462158e-08, 4.742340394655613e-09, -0.8114190126645531], 'iod': 0,
                             'a_gf0': 5153.798839569092, 'a_gf1': 1.6390180338742641,
                             'pos': [1.862645149230957e-09, -40.5625, 221.625],
                             'rate': [-2.2239983081817627e-06, 8.001923561096191e-06, 3.725290298461914e-09], 'toa': 0},
                     'valid': 1,
                     'sid': {'code': 0, 'sat': 2},
                     'toe': {'wn': 1876, 'tow': 172800.0},
                     'kepler': {'inc_dot': -4.475186409476941e-10, 'tgd': 1.862645149230957e-09, 'omegadot': -8.103551831174965e-09,
                                'sqrta': 5153.798839569092, 'inc': 0.9587534715647247, 'cus': 8.001923561096191e-06,
                                'omega0': 1.6390180338742641, 'cuc': -2.2239983081817627e-06, 'm0': -0.8114190126645531,
                                'toc': {'wn': 1876, 'tow': 172800.0}, 'dn': 4.742340394655613e-09, 'ecc': 0.00035172002390027046,
                                'cic': 3.725290298461914e-09, 'crs': -40.5625, 'iode': 59, 'iodc': 21845, 'cis': -6.146728992462158e-08,
                                'crc': 221.625, 'w': 2.9037284161724037, 'af0': -9.969808161258698e-07, 'af1': -5.229594535194337e-12, 'af2': 0.0},
                     'fit_interval': 4})]
  gpst = ti.GpsTime(wn=1876, tow=167132)
  sdiffs = o.make_propagated_sdiffs_(rover_nm, rover_nm, remote_dists, base_pos, es, gpst)
  assert len(sdiffs) == 2
  assert sdiffs[0].pseudorange > 0
  assert sdiffs[1].pseudorange > 0
  assert sdiffs[0].carrier_phase < 0
  assert sdiffs[0].doppler > 0
  assert sdiffs[0].sid['sat'] == 0
  assert sdiffs[1].sid['sat'] == 2
