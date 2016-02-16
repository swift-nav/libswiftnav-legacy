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

from swiftnav.observation import SingleDiff
import numpy as np
import pytest
import swiftnav.baseline as bl
import swiftnav.constants as c

# All these tests are ported from libswiftnav's
# tests/check_baseline.c, unless otherwise noted

def test_predict_carrier_obs():
  N = np.array([22.2, 23.3, 34.4, 123.4])
  DE = np.matrix([[1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1],
                  [1, 1, 1]])
  b = np.array([1, 1, 1])
  dd_obs = bl.predict_carrier_obs_(N, DE, b)
  dd_obs_expected = np.array([1.0/0.19023800915688557 + 22.2,
                              1.0/0.19023800915688557 + 23.3,
                              1.0/0.19023800915688557 + 34.4,
                              3.0/0.19023800915688557 + 123.4])
  assert np.allclose(dd_obs, dd_obs_expected)

def test_predict_carrier_obs2():
  N = np.array([222.2])
  DE = np.matrix([[3, 4, 5]])
  b = np.array([7, 8, 9])
  dd_obs_expected = np.array([(3*7 + 4*8 + 5*9)/0.19023800915688557 + 222.2])
  dd_obs = bl.predict_carrier_obs_(N, DE, b)
  assert np.allclose(dd_obs, dd_obs_expected)

def test_amb_from_baseline():
  N_true = np.array([22, 23, 34, -123])
  DE = np.matrix([[1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1],
                  [1, 1, 1]])
  b = np.array([1, 1, 1])
  dd_obs = bl.predict_carrier_obs_(N_true, DE, b)
  N = bl.amb_from_baseline_(DE, dd_obs, b)
  assert np.allclose(N, N_true)

def test_lesq_solution():
  # Over constrained.
  N = np.array([22, 23, 34, -123])
  DE = np.matrix([[1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1],
                  [1, 1, 1]])
  b_true = np.array([1.234, 1.456, 1.789])
  dd_obs = bl.predict_carrier_obs_(N, DE, b_true)
  code, b, resid = bl.lesq_solution_float_(dd_obs, N, DE)
  assert code == BASELINE_SUCCESS
  assert np.allclose(b, b_true)

def test_lesq_solution2():
  # Exactly constrained.
  N = np.array([22, 23, 34])
  DE = np.matrix([[1, 0, 0],
                  [0, 1, 0],
                  [1, 1, 1]])
  b_true = np.array([1.234, 1.456, 1.789])
  dd_obs = bl.predict_carrier_obs_(N, DE, b_true)
  code, b, resid = bl.lesq_solution_float_(dd_obs, N, DE)
  assert code == BASELINE_SUCCESS
  assert np.allclose(b, b_true)

# TODO (Buro): Add back
# def test_lesq_solution3():
#   # Under constrained, should fail with correct return code.
#   N[2]
#   DE[2*3]
#   b[3]
#   resid[2]
#   dd_obs[2]
#   code, b, resid = bl.lesq_solution_float_(dd_obs, N, DE)
#   assert code == BASELINE_NOT_ENOUGH_SATS_FLOAT

def test_lesq_solution4():
  # Over constrained, integer valued ambiguity
  N = np.array([22, 23, 34, -123])
  DE = np.matrix([[1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1],
                  [1, 1, 1]])
  b_true = np.array([1.234, 1.456, 1.789])
  dd_obs = bl.predict_carrier_obs_(N, DE, b_true)
  N_int = N
  code, num_used, residuals, removed_obs, b \
    = bl.lesq_solve_raim_(dd_obs, N_int, DE, False, bl.DEFAULT_RAIM_THRESHOLD_)
  assert code >= BASELINE_SUCCESS
  assert num_used == 4
  assert np.allclose(residuals,
                     np.array([  1.77635684e-15,  -3.55271368e-15,
                                 -5.32907052e-15,   3.55271368e-15]))
  assert removed_obs == 0
  assert np.allclose(b, np.array([ 1.234,  1.456,  1.789]))
  assert np.allclose(b, b_true)

# TODO (Buro): Crashes!
def test_lesq_solution5():
  # Over constrained with non-zero residuals
  N = np.array([0, 0, 0, 0])
  DE = np.matrix([[1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1],
                  [1, 0, 0]])
  b_true = np.array([0, 0, 0])
  dd_obs = np.array([0, 0, 0, 1])
  code, b, resid = bl.lesq_solution_float_(dd_obs, N, DE)
  b_expected = np.array([0.5*c.GPS_L1_LAMBDA_NO_VAC_, 0, 0])
  resid_expected = np.array([-0.5, 0, 0, 0.5])
  assert code == BASELINE_SUCCESS, "solution returned error %d" % code
  np.allclose(b, b_true)
  np.allclose(resid, resid_expected)

def baseline_setup():
  # Initialize sdiffs used in baseline() tests
  ref_ecef = np.array([0,0,0])
  sdiffs = [SingleDiff(sid={'sat': 1, 'code': 0},
                       pseudorange=0,
                       sat_pos=np.array([0, 0, 0]),
                       sat_vel=np.array([0, 0, 0]),
                       doppler=0,
                       carrier_phase=1,
                       snr=0.0),
            SingleDiff(sid={'sat':2, 'code': 0},
                       pseudorange=0,
                       sat_pos=np.array([1, 0, 0]),
                       sat_vel=np.array([0, 0, 0]),
                       doppler=0,
                       carrier_phase=2,
                       snr=0.0),
            SingleDiff(sid={'sat': 3, 'code': 0},
                       pseudorange=0,
                       sat_pos=np.array([0, 1, 0]),
                       sat_vel=np.array([0, 0, 0]),
                       doppler=0,
                       carrier_phase=3,
                       snr=0.0),
            SingleDiff(sid={'sat': 4, 'code': 0},
                       pseudorange=0,
                       sat_pos=np.array([0, 1, 1]),
                       sat_vel=np.array([0, 0, 0]),
                       doppler=0,
                       carrier_phase=4,
                       snr=0.0),
            SingleDiff(sid={'sat': 5, 'code': 0},
                       pseudorange=0,
                       sat_pos=np.array([0, 0, 1]),
                       sat_vel=np.array([0, 0, 0]),
                       doppler=0,
                       carrier_phase=5,
                       snr=0.0)]
  return (sdiffs, ref_ecef)

@pytest.mark.skipif(True, reason="Segfault!")
def test_baseline_ref_first():
  # Check that it works with the first sdiff as the reference
  # sat. This should verify that the loop can start correctly.
  ambs = bl.Ambiguities(n=4,
                        sids=[{'sat': 1, 'code': 0},
                              {'sat': 2, 'code': 0},
                              {'sat': 3, 'code': 0},
                              {'sat': 4, 'code': 0},
                              {'sat': 5, 'code': 0},
                              # Remaining
                              {'sat': 6, 'code': 0},
                              {'sat': 7, 'code': 0},
                              {'sat': 8, 'code': 0},
                              {'sat': 9, 'code': 0},
                              {'sat': 10, 'code': 0},
                              {'sat': 11, 'code': 0}],
                        ambs=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
  sdiffs, ref_ecef = baseline_setup()
  code, num_used, b = bl._baseline(sdiffs, ref_ecef, ambs, False, bl.DEFAULT_RAIM_THRESHOLD_)
  assert valid == BASELINE_SUCCESS
  assert num_used == 5
  assert np.allclose(b, np.array([-0.742242, -0.492905, -0.0533294]))

@pytest.mark.skipif(True, reason="Segfault!")
def test_baseline_ref_middle():
  # Check that it works with a middle sdiff as the reference sat. This
  # should verify that the induction works.
  sdiffs, ref_ecef = baseline_setup()
  ambs = bl.Ambiguities(n=4,
                        sids=[{'sat': 1, 'code': 0},
                              {'sat': 2, 'code': 0},
                              {'sat': 3, 'code': 0},
                              {'sat': 4, 'code': 0},
                              {'sat': 5, 'code': 0},
                              # Remaining
                              {'sat': 6, 'code': 0},
                              {'sat': 7, 'code': 0},
                              {'sat': 8, 'code': 0},
                              {'sat': 9, 'code': 0},
                              {'sat': 10, 'code': 0},
                              {'sat': 11, 'code': 0}],
                        ambs=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
  old_snr_2 = sdiffs[1].snr
  sdiffs[1] = SingleDiff(sid={'sat':2, 'code': 0},
                         pseudorange=0,
                         sat_pos=np.array([1, 0, 0]),
                         sat_vel=np.array([0, 0, 0]),
                         doppler=0,
                         carrier_phase=2,
                         snr=22.0)
  code, num_used, b = bl._baseline(sdiffs, ref_ecef, ambs, False, bl.DEFAULT_RAIM_THRESHOLD_)
  assert code == BASELINE_SUCCESS
  assert num_used == 5
  assert np.allclose(b, np.array([-0.622609, -0.432371, -0.00461595]))

@pytest.mark.skipif(True, reason="Segfault!")
def test_baseline_ref_end():
  # Check that it works with the last sdiff as the reference sat. This
  # should verify that the loop can terminate correctly
  sdiffs, ref_ecef = baseline_setup()
  ambs = bl.Ambiguities(n=4,
                        sids=[{'sat': 5, 'code': 0},
                              {'sat': 1, 'code': 0},
                              {'sat': 2, 'code': 0},
                              {'sat': 3, 'code': 0},
                              {'sat': 4, 'code': 0},
                              # Remaining
                              {'sat': 6, 'code': 0},
                              {'sat': 7, 'code': 0},
                              {'sat': 8, 'code': 0},
                              {'sat': 9, 'code': 0},
                              {'sat': 10, 'code': 0},
                              {'sat': 11, 'code': 0}],
                        ambs=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
  index = ambs.n
  old_snr_2 = sdiffs[index].snr
  s = sdiffs[index].to_dict()
  s.update({'snr': 22})
  sdiffs[index] = SingleDiff(**s)
  code, num_used, b = bl._baseline(sdiffs, ref_ecef, ambs, False, bl.DEFAULT_RAIM_THRESHOLD_)
  assert code == BASELINE_SUCCESS
  assert num_used == 5
  assert np.allclose(b, np.array([-0.589178, -0.35166, 0.0288157]))

# TODO (Buro): Add this back!
# def test_baseline_fixed_point():
#   # Check that measurements generated from a baseline result in an
#   # estimate matching the baseline.
#   sdiffs, ref_ecef = baseline_setup()
#   ambs = bl.Ambiguities(**{'n': 4,
#                         'sids': [{'sat': 5}, {'sat': 1}, {'sat': 2}, {'sat': 3}, {'sat': 4}],
#                         'ambs': np.array([0, 0, 0, 0])})
#   b_orig = np.array([1, 1, 1])
#   ref_ecef = np.array([0, 0, 0])
#   for (u8 i=0 i<5 i++) {
#     sdiffs[i].carrier_phase = vector_dot(3, b_orig, sdiffs[i].sat_pos) /
#                               vector_norm(3, sdiffs[i].sat_pos) /
#                               GPS_L1_LAMBDA_NO_VAC
#   code, num_used, b = bl._baseline(sdiffs, ref_ecef, ambs, False, bl.DEFAULT_RAIM_THRESHOLD_)
#   assert valid == 0
#   assert num_used == 5
#   assert np.allclose(b, b_orig)

def test_baseline_few_sats():
  sdiffs, ref_ecef = baseline_setup()
  ambs = bl.Ambiguities(n=0,
                        sids=[{'sat': 5, 'code': 0},
                              {'sat': 1, 'code': 0},
                              {'sat': 2, 'code': 0},
                              {'sat': 3, 'code': 0},
                              {'sat': 4, 'code': 0},
                              # Remaining
                              {'sat': 6, 'code': 0},
                              {'sat': 7, 'code': 0},
                              {'sat': 8, 'code': 0},
                              {'sat': 9, 'code': 0},
                              {'sat': 10, 'code': 0},
                              {'sat': 11, 'code': 0}],
                        ambs=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
  code, num_used, b = bl._baseline(sdiffs, ref_ecef, ambs, False, bl.DEFAULT_RAIM_THRESHOLD_)
  assert num_used == 0
  assert np.allclose(b, np.array([0, 0, 0]))
  assert code == BASELINE_NOT_ENOUGH_SATS_ROVER

def test_lesq_repair8():
  # Test raim repair: Over constrained with bad DE row.
  N = np.array([0,0,0,0,0,0,0,0])
  DE = np.matrix([[1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1],
                  [1, 1, 1],
                  [1, 1, 1],
                  [1, 1, 1],
                  [1, 1, 1],
                  [22, 222, 222]])
  dd_obs = np.array([1, 1, 1, 3, 3, 3, 3, 3])
  code, num_used, residuals, removed_obs, b \
    = bl.lesq_solve_raim_(dd_obs, N, DE, False, bl.DEFAULT_RAIM_THRESHOLD_)
  assert num_used == 7
  assert np.allclose(residuals, np.zeros(num_used))
  assert removed_obs == 7
  assert np.allclose(b, np.array([ 0.19023801,  0.19023801,  0.19023801]))
  assert code == BASELINE_SUCCESS_RAIM_REPAIR,
    "Expecting BASELINE_SUCCESS_RAIM_REPAIR for repaired solution, got: %i." % code
  assert removed_obs == 7, \
    "Expecting repaired solution (dropping index 4 of DE), got: %i." % removed_obs

def test_lesq_repair1():
  # Test raim repair: over constrained with bad DE row.
  N = np.array([0, 0, 0, 0, 0])
  DE = np.matrix([[1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1],
                  [1, 1, 1],
                  [22, 222, 222]])
  dd_obs = np.array([1, 1, 1, 3, 1])
  code, num_used, residuals, removed_obs, b \
    = bl.lesq_solve_raim_(dd_obs, N, DE, False, bl.DEFAULT_RAIM_THRESHOLD_)
  assert num_used == 4
  assert np.allclose(residuals, np.zeros(num_used))
  assert removed_obs == 4
  assert np.allclose(b, np.array([ 0.19023801,  0.19023801,  0.19023801]))
  assert code == BASELINE_SUCCESS_RAIM_REPAIR,
    "Expecting BASELINE_SUCCESS_RAIM_REPAIR for repaired solution, got: %i." % code
  assert removed_obs == 4, \
    "Expecting repaired solution (dropping index 4 of DE), got: %i." % removed_obs

def test_lesq_repair_disabled():
  # Test raim disabling flag: Over constrained with bad DE row.
  N = np.array([0, 0, 0, 0, 0])
  DE = np.matrix([[1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1],
                  [1, 1, 1],
                  [22, 222, 222]])
  dd_obs = np.array([1, 1, 1, 3, 1])
  # DISABLE raim
  code, num_used, residuals, removed_obs, b \
    = bl.lesq_solve_raim_(dd_obs, N, DE, True, bl.DEFAULT_RAIM_THRESHOLD_)
  assert np.allclose(b, np.array([ 0.37698481, -0.01824651, -0.01824651]))
  assert num_used == 5
  assert np.allclose(residuals,
                     np.array([-0.9816482, 1.09591413,  1.09591413,  1.21018006, -0.01038781]))
  assert removed_obs == 0
  assert code == BASELINE_SUCCESS_NO_RAIM,
    "Expecting BASELINE_SUCCESS_NO_RAIM for repaired solution, got: %i." % code

def test_lesq_repair2():
  # Bad DE row, not enough rows to repair.
  N = np.array([0, 0, 0, 0])
  DE = np.matrix([[1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1],
                  [22, 222, 222]])
  dd_obs = np.array([1, 1, 1, 1])
  code, num_used, residuals, removed_obs, b \
    = bl.lesq_solve_raim_(dd_obs, N, DE, False, bl.DEFAULT_RAIM_THRESHOLD_)
  assert np.allclose(b, np.array([ 0.1705906,  -0.00802221, -0.00802221]))
  assert num_used == 0
  assert np.allclose(residuals, np.zeros(num_used))
  assert removed_obs == 0
  assert code == BASELINE_NOT_ENOUGH_SATS_RAIM,
    "Expecting BASELINE_NOT_ENOUGH_SATS_RAIM for not enough dds to repair, got: %i." % code
