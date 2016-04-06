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
import swiftnav.ephemeris as e
import swiftnav.time as t

tol = 1e-10

def test_sat_state():
  eph = e.Ephemeris(**{'sid': {'sat': 17, 'code': 0,},
                       'toe': {'wn': 1867, 'tow': 518400.0,},
                       'ura': 2.0,
                       'fit_interval': 4 * 60 * 60,
                       'valid': 1,
                       'healthy': 1,
                       'kepler': {
                                  'toc': {'tow': 518400.0, 'wn': 1867,},
                                  'crs': 25.125,
                                  'inc_dot': 3.78944355982045e-10,
                                  'tgd': -1.1175870895385742e-08,
                                  'crc': 106.65625,
                                  'ecc': 0.016364791779778898,
                                  'omegadot': -8.013190924423338e-09,
                                  'inc': 0.9253317285121154,
                                  'cuc': 1.255422830581665e-06,
                                  'omega0': 2.7384009602031045,
                                  'cus': 1.280754804611206e-05,
                                  'm0': -2.057975194561658,
                                  'dn': 5.035924052164783e-09,
                                  'cic': 2.7194619178771973e-07,
                                  'sqrta': 5153.647108078003,
                                  'cis': 9.313225746154785e-09,
                                  'iode': 50,
                                  'iodc': 21845,
                                  'crc': 49614,
                                  'w': -1.9329047030450934,
                                  'af0': 0.0004458986222743988,
                                  'af1': 3.637978807091713e-12,
                                  'af2': 0.0}})
  assert eph.is_valid(t.GpsTime(**{ 'wn': 1867, 'tow': 518400.0,}))
  assert eph.is_healthy()
  time = t.GpsTime(**{ 'wn': 1867, 'tow': 518400.0,})
  pos, vel, clock_err, clock_rate_err =  eph.calc_sat_state(time)
  assert np.allclose(pos, [8939514.51292471, -19122117.65874583, 16447534.13118957])
  assert np.allclose(vel, [1892.86704302, -765.38725981, -1985.61338019])
  assert clock_err - 0.00336435649383 < tol
  assert clock_rate_err - 3.63797880709e-12 < tol
  gnd_pos = np.array([1115054, -4843905,  3983584])
  doppler = eph.calc_doppler(time, gnd_pos)
  assert np.isclose(doppler, 253.83513370122773)
  azimuth, elevation = eph.calc_az_el(time, gnd_pos)
  assert np.isclose(azimuth, 1.5874048317078373)
  assert np.isclose(elevation, 1.3533993525984187)
  assert not eph.is_valid(t.GpsTime(**{ 'wn': 1866, 'tow': 518400.0,}))
  assert not eph.is_valid(t.GpsTime(**{ 'wn': 1867, 'tow': 518400.0 + 3600.*2 + 1}))
  assert not eph.is_valid(t.GpsTime(**{ 'wn': 1867, 'tow': 518400.0 - 3600.*2 - 1}))
  assert eph.is_valid(t.GpsTime(**{ 'wn': 1867, 'tow': 518400.0 + 3600.*1}))
  assert eph.is_valid(t.GpsTime(**{ 'wn': 1867, 'tow': 518400.0 - 3600.*1}))
  assert not eph.is_valid(t.GpsTime(**{ 'wn': 1868, 'tow': 518400.0,}))
  assert len(eph.to_dict()) == 9
  assert repr(eph)
  eph = e.Ephemeris(**{'sid': {'sat': 17, 'code': 0},
                     'toe': { 'wn': 1867, 'tow': 518400.0,},
                     'ura': 2.0,
                     'fit_interval': 4 * 60 * 60,
                     'valid': 1,
                     'healthy': 0,
                     'kepler': {'crs': 25.125,
                                'inc_dot': 3.78944355982045e-10,
                                'tgd': -1.1175870895385742e-08,
                                'crc': 106.65625,
                                'ecc': 0.016364791779778898,
                                'omegadot': -8.013190924423338e-09,
                                'inc': 0.9253317285121154,
                                'cuc': 1.255422830581665e-06,
                                'omega0': 2.7384009602031045,
                                'cus': 1.280754804611206e-05,
                                'm0': -2.057975194561658,
                                'toc': {'tow': 518400.0, 'wn': 1867,},
                                'dn': 5.035924052164783e-09,
                                'cic': 2.7194619178771973e-07,
                                'sqrta': 5153.647108078003,
                                'cis': 9.313225746154785e-09,
                                'iode': 50,
                                'iodc': 21845,
                                'crc': 49614,
                                'w': -1.9329047030450934,
                                'af0': 0.0004458986222743988,
                                'af1': 3.637978807091713e-12,
                                'af2': 0.0}})
  assert eph.is_valid(t.GpsTime(**{ 'wn': 1867, 'tow': 518400.0,}))
  assert not eph.is_healthy()
  eph = e.Ephemeris(**{'sid': {'sat': 17, 'code': 0},
                     'toe': { 'wn': 1867, 'tow': 518400.0,},
                     'ura': 2.0,
                     'fit_interval': 4 * 60 * 60,
                     'valid': 0,
                     'healthy': 0,
                     'kepler': {'crs': 25.125,
                                'inc_dot': 3.78944355982045e-10,
                                'tgd': -1.1175870895385742e-08,
                                'crc': 106.65625,
                                'ecc': 0.016364791779778898,
                                'omegadot': -8.013190924423338e-09,
                                'inc': 0.9253317285121154,
                                'cuc': 1.255422830581665e-06,
                                'omega0': 2.7384009602031045,
                                'cus': 1.280754804611206e-05,
                                'm0': -2.057975194561658,
                                'toc': {'tow': 518400.0, 'wn': 1867,},
                                'dn': 5.035924052164783e-09,
                                'cic': 2.7194619178771973e-07,
                                'sqrta': 5153.647108078003,
                                'cis': 9.313225746154785e-09,
                                'iode': 50,
                                'iodc': 21845,
                                'crc': 49614,
                                'w': -1.9329047030450934,
                                'af0': 0.0004458986222743988,
                                'af1': 3.637978807091713e-12,
                                'af2': 0.0}})
  assert not eph.is_valid(t.GpsTime(**{ 'wn': 1867, 'tow': 518400.0,}))
  assert eph.is_healthy()
  eph = e.Ephemeris(**{'sid': {'sat': 17, 'code': 0},
                         'toe': { 'wn': 1867, 'tow': 518400.0,},
                         'ura': 2.0,
                         'fit_interval': 6 * 60 * 60,
                         'valid': 1,
                         'healthy': 1,
                         'kepler': {'crs': 25.125,
                                    'inc_dot': 3.78944355982045e-10,
                                    'tgd': -1.1175870895385742e-08,
                                    'crc': 106.65625,
                                    'ecc': 0.016364791779778898,
                                    'omegadot': -8.013190924423338e-09,
                                    'inc': 0.9253317285121154,
                                    'cuc': 1.255422830581665e-06,
                                    'omega0': 2.7384009602031045,
                                    'cus': 1.280754804611206e-05,
                                    'm0': -2.057975194561658,
                                    'toc': {'tow': 518400.0, 'wn': 1867,},
                                    'dn': 5.035924052164783e-09,
                                    'cic': 2.7194619178771973e-07,
                                    'sqrta': 5153.647108078003,
                                    'cis': 9.313225746154785e-09,
                                    'iode': 50,
                                    'iodc': 21845,
                                    'crc': 49614,
                                    'w': -1.9329047030450934,
                                    'af0': 0.0004458986222743988,
                                    'af1': 3.637978807091713e-12,
                                    'af2': 0.0}})
  assert eph.is_valid(t.GpsTime(**{ 'wn': 1867, 'tow': 518400.0 + 3600.*2}))
  assert not eph.is_valid(t.GpsTime(**{ 'wn': 1867, 'tow': 518400.0 + 3600.*3 + 1}))
  assert eph.is_healthy()

def test_parameters():
  toe = t.GpsTime(**{ 'wn': 1867, 'tow': 518400.0,})
  time = t.GpsTime(**{ 'wn': 1867, 'tow': 518400.0,})
  valid = 1
  fit_interval = 4 * 60 * 60
  assert e.Ephemeris.is_params_valid(valid, fit_interval, toe, time)

  toe = t.GpsTime(**{ 'wn': 1867, 'tow': 518400.0,})
  time = t.GpsTime(**{ 'wn': 1867, 'tow': 518400.0,})
  valid = 0
  fit_interval = 4 * 60 * 60
  assert not e.Ephemeris.is_params_valid(valid, fit_interval, toe, time)

  toe = t.GpsTime(**{ 'wn': 1867, 'tow': 518400.0,})
  time = t.GpsTime(**{ 'wn': 1867, 'tow': 518400.0,})
  valid = 1
  fit_interval = 0 * 60 * 60
  assert not e.Ephemeris.is_params_valid(valid, fit_interval, toe, time)

  toe = t.GpsTime(**{ 'wn': 1867, 'tow': 518400.0,})
  time = t.GpsTime(**{ 'wn': 1867, 'tow': 518400.0 + 3600.0 * 2,})
  valid = 1
  fit_interval = 4 * 60 * 60
  assert e.Ephemeris.is_params_valid(valid, fit_interval, toe, time)

  toe = t.GpsTime(**{ 'wn': 1867, 'tow': 518400.0,})
  time = t.GpsTime(**{ 'wn': 1867, 'tow': 518400.0 + 3600.0 * 3 ,})
  valid = 1
  fit_interval = 4 * 60 * 60
  assert not e.Ephemeris.is_params_valid(valid, fit_interval, toe, time)

def test_glo_sat_state():
  eph1 = e.Ephemeris(**{'sid': {'sat': 3, 'code': 3,},
                     'toe': {'wn': 1892, 'tow': 303900.0,}, #from GPS
                     'ura': 2.0,
                     'fit_interval': 4 * 60 * 60,
                     'valid': 1,
                     'healthy': 1,
                     'glo': {  #data at 2016-04-13T12:25:00
                              'gamma' : 9.09494701772928238e-13,
                              'tau' : -6.72144815325737000e-05,
                              'pos' : [-1.4379786621093750e+07,-7.6587543945312500e+06,1.9676794921875000e+07],
                              'vel' : [2.6736078262329102e+03,-2.4776077270507813e+02,1.8579368591308594e+03],
                              'acc' : [0.0,0.0,-2.79396772384643555e-06],
                            }})
  eph2 = e.Ephemeris(**{'sid': {'sat': 3, 'code': 3,},
                     'toe': {'wn': 1892, 'tow': 305700.0,}, #from GPS
                     'ura': 2.0,
                     'fit_interval': 4 * 60 * 60,
                     'valid': 1,
                     'healthy': 1,
                     'glo': {  #data at 2016-04-13T12:55:00
                              'gamma' : 9.09494701772928238e-13,
                              'tau' : -6.72172755002975464e-05,
                              'pos' : [-9.2793393554687500e+06,-8.5263706054687500e+06,2.2220868652343750e+07],
                              'vel' : [2.9441843032836914e+03,-7.2328472137451172e+02,9.5054054260253906e+02],
                              'acc' : [0.0,0.0,-2.79396772384643555e-06],
                            }})
  #check that both ephemeris are valid at 2016-04-13T12:40:00
  assert eph1.is_valid(t.GpsTime(**{ 'wn': 1892, 'tow': 304800.0,}))
  assert eph2.is_valid(t.GpsTime(**{ 'wn': 1892, 'tow': 304800.0,}))
  assert eph1.is_healthy()
  assert eph2.is_healthy()
  #set time in between 2 ephemeris
  time = t.GpsTime(**{ 'wn': 1892, 'tow': 304800.0,})
  #calculate SV orbits at the time
  pos1, vel1, clock_err1, clock_rate_err1 =  eph1.calc_sat_state(time)
  pos2, vel2, clock_err2, clock_rate_err2 =  eph2.calc_sat_state(time)
  #check if position difference of the SV in 2.5 m range
  assert np.allclose(pos1, pos2, 0, 2.5)
  #check if velocity difference of the SV in 0.005 m/s range
  assert np.allclose(vel1, vel2, 0, 0.005)
  assert abs(clock_err1 - clock_err2) < 1.0e-8
  assert abs(clock_rate_err1 - clock_rate_err2) < 1.0e-8
  # Now check GLO orbit computation (pos and vel) at particular time
  # set reference ephemeris at 2007-15-11T06.15.00
  eph = e.Ephemeris(**{'sid': {'sat': 3, 'code': 3,},
                    'toe': {'wn': 1453, 'tow': 368100.0,}, #from GPS
                    'ura': 2.0,
                    'fit_interval': 4 * 60 * 60,
                    'valid': 1,
                    'healthy': 1,
                    'glo': {  #input data from GLO ICD, pg. 55
                             'gamma' : 0, #not necessary for the test
                             'tau' : 0, #not necessary for the test
                             'pos' : [-14081752.701,18358958.252,10861302.124],
                             'vel' : [-1025.76358,1086.72147,-3157.32343],
                             'acc' : [1.7156e-6-9.2581e-7,1.0278e-6-1.0343e-6,-1.0368e-6-1.1260e-6],
                          }})
  #set time to 6:30
  time = t.GpsTime(**{ 'wn': 1453, 'tow': 369000.0,})
  #calculate SV orbits at the time
  pos, vel, clock_err, clock_rate_err =  eph.calc_sat_state(time)
  #check calculated vel and pos data at 2007-15-11T06:30:00
  #reference values from GLO ICD, pg 55. Note: it seems there is an typo in
  # result Vz in ICD example
  assert np.allclose(pos, [-14836563.872,19249935.476,7924017.196], 0, 2.5)
  assert np.allclose(vel, [-653.97782,882.62958,-3359.44444], 0, 0.005)
