# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

cimport numpy as np
from almanac cimport *
from amb_kf cimport *
from amb_kf import KalmanFilter
from ambiguity_test cimport *
from ambiguity_test import AmbiguityTest
from baseline cimport *
from constants cimport MAX_SATS
from time cimport *
from libc.stdio cimport printf
from libc.string cimport memcpy
from observation cimport *
from observation cimport SingleDiff
from observation import SingleDiff
from sats_management cimport *
from signal cimport *
import numpy as np


def dgnss_set_settings_(phase_var_test, code_var_test,
                        phase_var_kf, code_var_kf,
                        amb_drift_var, amb_init_var, new_int_var):
  dgnss_set_settings(phase_var_test, code_var_test,
                     phase_var_kf, code_var_kf,
                     amb_drift_var, amb_init_var,new_int_var)

def make_measurements_(sdiffs):
  num_ddiffs = len(sdiffs) - 1
  cdef sdiff_t sdiffs_[32]
  mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] raw_measurements = np.empty(2 * num_ddiffs, dtype=np.double)
  make_measurements(num_ddiffs, &sdiffs_[0], &raw_measurements[0])
  return raw_measurements

def dgnss_init_(sdiffs, reciever_ecef):
  num_sdiffs = len(sdiffs)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(reciever_ecef, dtype=np.double)
  cdef sdiff_t sdiffs_[32]
  mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
  dgnss_init(num_sdiffs, &sdiffs_[0], &ref_ecef_[0])

def dgnss_update_(sdiffs, reciever_ecef, disable_raim=False):
  num_sdiffs = len(sdiffs)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(reciever_ecef, dtype=np.double)
  cdef sdiff_t sdiffs_[32]
  cdef sdiff_t s_
  for (i,sdiff) in enumerate(sdiffs):
    s_ = (<SingleDiff ?> sdiff)._thisptr
    memcpy(&sdiffs_[i], &s_, sizeof(sdiff_t))
  dgnss_update(num_sdiffs, &sdiffs_[0], &ref_ecef_[0], disable_raim, DEFAULT_RAIM_THRESHOLD)

# def dgnss_rebase_ref_(sdiffs, reciever_ecef, old_prns):
#   num_sdiffs = len(sdiffs)
#   cdef sdiff_t sdiffs_[32], corrected_sdiffs_[32]
#   mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
#   mk_sdiff_array(sdiffs, 32, &corrected_sdiffs_[0])
#   cdef np.ndarray[np.double_t, ndim=1, mode="c"] receiver_ecef_ = np.array(reciever_ecef, dtype=np.double)
#   cdef np.ndarray[np.uint8_t, ndim=1, mode="c"] old_prns_ = np.array(old_prns, dtype=np.uint8)
#   dgnss_rebase_ref(num_sdiffs, &sdiffs_[0], &receiver_ecef_[0], &old_prns_[0], &corrected_sdiffs_[0])
#   raise NotImplementedError()

def get_dgnss_nkf_():
  cdef nkf_t* nkf = get_dgnss_nkf()
  return

# def get_sats_management_():
#   cdef sats_management_t * sats_man = get_sats_management()
#   cdef np.ndarray[np.uint8_t, ndim=1, mode="c"] prns = np.empty(sats_man.num_sats, dtype=np.uint8)
#   memcpy(&prns[0], sats_man.prns, sats_man.num_sats * sizeof(u8))
#   return (sats_man.num_sats, prns)

def get_ambiguity_test_():
  cdef ambiguity_test_t* test = get_ambiguity_test()
  amb_test = AmbiguityTest()
  memcpy(test, &amb_test._thisptr, sizeof(ambiguity_test_t))
  return amb_test

def dgnss_iar_resolved_():
  return dgnss_iar_resolved() > 0

def dgnss_iar_num_hyps_():
  return dgnss_iar_num_hyps()

def dgnss_iar_num_sats_():
  return dgnss_iar_num_sats_()

def dgnss_iar_get_single_hyp_(num_dds):
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] hyp = np.empty(num_dds, dtype=np.double)
  dgnss_iar_get_single_hyp(&hyp[0])
  return hyp

def dgnss_init_known_baseline_(sdiffs, receiver_ecef, b):
  num_sats = len(sdiffs)
  cdef sdiff_t sdiffs_[32]
  mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] receiver_ecef_ = np.array(receiver_ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b_ = np.array(b, dtype=np.double)
  dgnss_init_known_baseline(num_sats, &sdiffs_[0], &receiver_ecef_[0], &b_[0])

def dgnss_update_ambiguity_state_(AmbiguityState s):
  dgnss_update_ambiguity_state(&s._thisptr)

def dgnss_baseline_(sdiffs, ref_ecef, AmbiguityState s, disable_raim=False):
  cdef u8 num_sdiffs = len(sdiffs)
  cdef u8 num_used
  cdef sdiff_t sdiffs_[32]
  mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(ref_ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = np.empty(3, dtype=np.double)
  cdef s8 flag = dgnss_baseline(num_sdiffs, &sdiffs_[0], &ref_ecef_[0], &s._thisptr,
                                &num_used, &b[0], disable_raim, DEFAULT_RAIM_THRESHOLD)
  return (flag, num_used, b)


def dgnss_fixed_baseline(sdiffs, ref_ecef, AmbiguityState s, disable_raim=False):
  cdef num_sdiffs = len(sdiffs)
  cdef u8 num_used
  cdef sdiff_t sdiffs_[32]
  mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(ref_ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = np.empty(3, dtype=np.double)
  cdef s8 flag = baseline(num_sdiffs, &sdiffs_[0], &ref_ecef_[0],
                          &s._thisptr.fixed_ambs, &num_used, &b[0],
                          disable_raim, DEFAULT_RAIM_THRESHOLD)
  return (flag, num_used, b)

def dgnss_float_baseline(sdiffs, ref_ecef, AmbiguityState s):
  cdef num_sdiffs = len(sdiffs)
  cdef u8 num_used
  cdef sdiff_t sdiffs_[32]
  mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(ref_ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = np.empty(3, dtype=np.double)
  cdef s8 flag = baseline(num_sdiffs, &sdiffs_[0], &ref_ecef_[0],
                          &s._thisptr.float_ambs, &num_used, &b[0],
                          False, DEFAULT_RAIM_THRESHOLD)
  return flag, num_used, b

def measure_amb_kf_b_(sdiffs, receiver_ecef, AmbiguityState s, disable_raim=False):
  cdef u8 num_sdiffs = len(sdiffs)
  cdef sdiff_t sdiffs_[32]
  mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] receiver_ecef_ = np.array(receiver_ecef, dtype=np.double)
  cdef double b
  measure_amb_kf_b(num_sdiffs, &sdiffs_[0], &receiver_ecef_[0], &b)
  return b

# TODO (Buro): Fix later
#TODO eventually, want to get reciever_ecef from data
# def measure_float_b(sdiffs, reciever_ecef):
#   num_sdiffs = len(sdiffs)
#   cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(reciever_ecef, dtype=np.double)
#   cdef sdiff_t sdiffs_[32]
#   mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
#   cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = np.empty(3, dtype=np.double)
#   measure_amb_kf_b(&ref_ecef_[0], num_sdiffs, &sdiffs_[0], &b[0])
#   return b

def measure_b_with_external_ambs_(sdiffs, ambs, reciever_ecef): #TODO eventually, want to get reciever_ecef from data
  num_sdiffs = len(sdiffs)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(reciever_ecef, dtype=np.double)
  cdef sdiff_t sdiffs_[32]
  mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ambs_ = np.array(ambs, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = np.empty(3, dtype=np.double)
  measure_b_with_external_ambs(len(ambs), &ambs_[0], num_sdiffs, &sdiffs_[0], &ref_ecef_[0], &b[0])
  return b

def measure_iar_b_with_external_ambs_(sdiffs, ambs, reciever_ecef): #TODO eventually, want to get reciever_ecef from data
  num_sdiffs = len(sdiffs)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(reciever_ecef, dtype=np.double)
  cdef sdiff_t sdiffs_[32]
  mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ambs_ = np.array(ambs, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = np.empty(3, dtype=np.double)
  measure_iar_b_with_external_ambs(&ref_ecef_[0], num_sdiffs, &sdiffs_[0], &ambs_[0], &b[0])
  return b

def get_amb_kf_de_and_phase_(sdiffs, ref_ecef):
  num_sdiffs = len(sdiffs)
  cdef sdiff_t sdiffs_[32]
  mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
  cdef double de, phase
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(ref_ecef, dtype=np.double)
  cdef u8 num_sats = get_amb_kf_de_and_phase(num_sdiffs, &sdiffs_[0], &ref_ecef_[0], &de, &phase)
  return (num_sats, de, phase)

def get_iar_de_and_phase_(sdiffs, ref_ecef):
  num_sdiffs = len(sdiffs)
  cdef sdiff_t sdiffs_[32]
  mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
  cdef double de, phase
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(ref_ecef, dtype=np.double)
  cdef u8 num_sats = get_iar_de_and_phase(num_sdiffs, &sdiffs_[0], &ref_ecef_[0], &de, &phase)
  return (num_sats, de, phase)

def dgnss_iar_pool_contains_(ambs):
  num_ambs = len(ambs)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ambs_ = np.array(ambs, dtype=np.double)
  return dgnss_iar_pool_contains(&ambs_[0])

def dgnss_iar_pool_ll_(ambs):
  num_ambs = len(ambs)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ambs_ = np.array(ambs, dtype=np.double)
  return dgnss_iar_pool_ll(num_ambs, &ambs_[0])

def dgnss_iar_pool_prob_(ambs):
  num_ambs = len(ambs)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ambs_ = np.array(ambs, dtype=np.double)
  return dgnss_iar_pool_prob(num_ambs, &ambs_[0])

def get_float_de_and_phase(sdiffs, ref_ecef):
  num_sdiffs = len(sdiffs)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(ref_ecef, dtype=np.double)
  cdef sdiff_t sdiffs_[32]
  mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] de = np.empty((32,3), dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] phase = np.empty(32, dtype=np.double)
  m = get_amb_kf_de_and_phase(num_sdiffs, &sdiffs_[0], &ref_ecef_[0], &de[0,0], &phase[0])
  return (de[:m-1], phase[:m-1])

def get_iar_de_and_phase_(sdiffs, ref_ecef):
  num_sdiffs = len(sdiffs)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(ref_ecef, dtype=np.double)
  cdef sdiff_t sdiffs_[32]
  mk_sdiff_array(sdiffs, 32, &sdiffs_[0])
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] de = np.empty((32,3), dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] phase = np.empty(32, dtype=np.double)
  m = get_iar_de_and_phase(num_sdiffs, &sdiffs_[0], &ref_ecef_[0], &de[0,0], &phase[0])
  return (de[:m-1], phase[:m-1])

def dgnss_iar_pool_contains_(ambs):
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ambs_ = np.array(ambs, dtype=np.double)
  return dgnss_iar_pool_contains(&ambs_[0]) == 1

def get_amb_kf_mean_():
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ambs = np.empty(32, dtype=np.double)
  num_dds = get_amb_kf_mean(&ambs[0])
  return ambs[:num_dds]

def get_amb_kf_cov_(num_dds):
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] cov = np.empty((num_dds, num_dds), dtype=np.double)
  num_dds2 = get_amb_kf_cov(&cov[0,0])
  if num_dds == num_dds2:
    return cov
  else:
    raise ValueError("Was given the wrong num_dds.")

def get_amb_kf_cov2():
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] cov = np.empty((32*32, 1), dtype=np.double)
  num_dds = get_amb_kf_cov(&cov[0,0])
  return cov[:num_dds*num_dds,0].reshape((num_dds,num_dds))

# def get_amb_kf_prns_():
#   cdef np.ndarray[np.uint8_t, ndim=1, mode="c"] prns = np.empty(32, dtype=np.uint8)
#   num_sats = get_amb_kf_prns(&prns[0])
#   return prns[:num_sats]

# def get_amb_test_prns_():
#   cdef np.ndarray[np.uint8_t, ndim=1, mode="c"] prns = np.empty(32, dtype=np.uint8)
#   num_sats = get_amb_test_prns(&prns[0])
#   return prns[:num_sats]

def dgnss_iar_MLE_ambs_():
  cdef np.ndarray[np.int32_t, ndim=1, mode="c"] ambs = np.empty(32, dtype=np.int32)
  num_dds = dgnss_iar_MLE_ambs(&ambs[0])
  return ambs[:num_dds]
