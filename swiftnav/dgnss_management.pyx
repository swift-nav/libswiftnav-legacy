# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

import numpy as np
cimport numpy as np
cimport dgnss_management_c
from amb_kf import KalmanFilter
from amb_kf cimport *
from amb_kf_c cimport *
from single_diff_c cimport *
from single_diff import SingleDiff
from single_diff cimport SingleDiff
from almanac cimport *
from almanac_c cimport *
from gpstime cimport *
from gpstime_c cimport *
from libc.string cimport memcpy
from libc.stdio cimport printf
from sats_management_c cimport *


def set_settings(phase_var_test, code_var_test,
                 phase_var_kf, code_var_kf,
                 amb_drift_var,
                 amb_init_var,
                 new_int_var):
  dgnss_management_c.dgnss_set_settings(phase_var_test, code_var_test,
                                        phase_var_kf, code_var_kf,
                                        amb_drift_var,
                                        amb_init_var,
                                        new_int_var)

def dgnss_init(sdiffs,
               reciever_ecef):
  num_sdiffs = len(sdiffs)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(reciever_ecef, dtype=np.double)
  cdef sdiff_t sdiffs_[32]
  cdef sdiff_t s_
  for (i,sdiff) in enumerate(sdiffs):
    s_ = (<SingleDiff> sdiff).sdiff
    memcpy(&sdiffs_[i], &s_, sizeof(sdiff_t))

  dgnss_management_c.dgnss_init(num_sdiffs, &sdiffs_[0], &ref_ecef_[0])


def dgnss_update(sdiffs,
                 reciever_ecef):
  num_sdiffs = len(sdiffs)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(reciever_ecef, dtype=np.double)
  cdef sdiff_t sdiffs_[32]
  cdef sdiff_t s_
  for (i,sdiff) in enumerate(sdiffs):
    s_ = (<SingleDiff> sdiff).sdiff
    memcpy(&sdiffs_[i], &s_, sizeof(sdiff_t))

  dgnss_management_c.dgnss_update(num_sdiffs, &sdiffs_[0], &ref_ecef_[0])


def alm_dgnss_init(alms, GpsTime timestamp,
               numpy_measurements,
               reciever_ecef, b=None):
  n = len(alms)
  state_dim = n + 5
  obs_dim = 2 * (n-1)

  cdef almanac_t al[32]
  cdef almanac_t a_
  cdef np.ndarray[np.uint8_t, ndim=1, mode="c"] prns = \
        np.empty(n, dtype=np.uint8)
  for i, a in enumerate(alms):
    a_ = (<Almanac> a).almanac
    memcpy(&al[i], &a_, sizeof(almanac_t))
    prns[i] = (<Almanac> a).prn

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(reciever_ecef, dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b_
  if (b == None):
    b_ = np.empty(3, dtype=np.double)
  else:
    b_ = np.array(b, dtype=np.double)

  cdef gps_time_t timestamp_ = timestamp.gps_time

  cdef sdiff_t sdiffs[32]
  almanacs_to_single_diffs(len(alms), &al[0], timestamp_, sdiffs)

  for i, (l,c) in enumerate(numpy_measurements):
    sdiffs[i].pseudorange = c
    sdiffs[i].carrier_phase = l

  # dgnss_management_c.dgnss_init(n, &sdiffs[0], &ref_ecef_[0], &b_[0])
  dgnss_management_c.dgnss_init(n, &sdiffs[0], &ref_ecef_[0])

def alm_dgnss_update(alms, GpsTime timestamp,
               numpy_measurements,
               reciever_ecef):
  n = len(alms)
  state_dim = n + 5
  obs_dim = 2 * (n-1)

  cdef almanac_t al[32]
  cdef almanac_t a_
  cdef np.ndarray[np.uint8_t, ndim=1, mode="c"] prns = \
        np.empty(n, dtype=np.uint8)
  for i, a in enumerate(alms):
    a_ = (<Almanac> a).almanac
    memcpy(&al[i], &a_, sizeof(almanac_t))
    prns[i] = (<Almanac> a).prn

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(reciever_ecef, dtype=np.double)

  cdef gps_time_t timestamp_ = timestamp.gps_time

  cdef sdiff_t sdiffs[32]
  almanacs_to_single_diffs(len(alms), &al[0], timestamp_, sdiffs)

  for i, (l,c) in enumerate(numpy_measurements):
    sdiffs[i].pseudorange = c
    sdiffs[i].carrier_phase = l

  cdef double b_[3]
  # dgnss_management_c.dgnss_update(n, &sdiffs[0], &ref_ecef_[0], dt, 1, b_)
  # print n
  dgnss_management_c.dgnss_update(n, &sdiffs[0], &ref_ecef_[0])
  
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = \
    np.empty(3, dtype=np.double)
  memcpy(&b[0],b_, 3*sizeof(double))
  return b

def dgnss_iar_resolved():
  return dgnss_management_c.dgnss_iar_resolved() > 0

def dgnss_iar_num_hyps():
  return dgnss_management_c.dgnss_iar_num_hyps()
  
def dgnss_iar_num_sats():
  return dgnss_management_c.dgnss_iar_num_sats()

def dgnss_iar_get_single_hyp(num_dds):
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] hyp = \
    np.empty(num_dds, dtype=np.double)
  dgnss_management_c.dgnss_iar_get_single_hyp(&hyp[0])
  return hyp

def get_sats_management():
  cdef sats_management_t * sats_man = dgnss_management_c.get_sats_management()
  cdef np.ndarray[np.uint8_t, ndim=1, mode="c"] prns = \
    np.empty(sats_man.num_sats, dtype=np.uint8)
  memcpy(&prns[0], sats_man.prns, sats_man.num_sats * sizeof(u8))
  return sats_man.num_sats, prns

def dgnss_new_float_baseline(sdiffs, ref_ecef):
  cdef u8 num_sats = len(sdiffs)
  cdef u8 num_used
  cdef sdiff_t sdiffs_[32]
  cdef sdiff_t s_
  for (i,sdiff) in enumerate(sdiffs):
    s_ = (<SingleDiff> sdiff).sdiff
    memcpy(&sdiffs_[i], &s_, sizeof(sdiff_t))
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(ref_ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = \
    np.empty(3, dtype=np.double)
  dgnss_management_c.dgnss_new_float_baseline(num_sats,
                                              &sdiffs_[0], &ref_ecef_[0],
                                              &num_used, &b[0])
  return num_used, b

def dgnss_fixed_baseline(sdiffs, ref_ecef):
  cdef u8 num_sats = len(sdiffs)
  cdef u8 num_used
  cdef sdiff_t sdiffs_[32]
  cdef sdiff_t s_
  for (i,sdiff) in enumerate(sdiffs):
    s_ = (<SingleDiff> sdiff).sdiff
    memcpy(&sdiffs_[i], &s_, sizeof(sdiff_t))
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(ref_ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = \
    np.empty(3, dtype=np.double)
  cdef s8 flag = dgnss_management_c.dgnss_fixed_baseline(num_sats,
                                                         &sdiffs_[0], &ref_ecef_[0],
                                                         &num_used, &b[0])
  if flag == 1:
    return num_used, b
  else:
    return 0, np.array([np.nan]*3)

def measure_float_b(sdiffs, reciever_ecef): #TODO eventually, want to get reciever_ecef from data
  num_sdiffs = len(sdiffs)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(reciever_ecef, dtype=np.double)
  cdef sdiff_t sdiffs_[32]
  cdef sdiff_t s_
  for (i,sdiff) in enumerate(sdiffs):
    s_ = (<SingleDiff> sdiff).sdiff
    memcpy(&sdiffs_[i], &s_, sizeof(sdiff_t))
  
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = \
    np.empty(3, dtype=np.double)

  dgnss_management_c.measure_amb_kf_b(&ref_ecef_[0],
                                      num_sdiffs, &sdiffs_[0],
                                      &b[0])
  return b

def alm_measure_float_b(alms, GpsTime timestamp,
                numpy_measurements,
                reciever_ecef): #TODO eventually, want to get reciever_ecef from data
  n = len(alms)
  cdef almanac_t al[32]
  cdef almanac_t a_
  cdef np.ndarray[np.uint8_t, ndim=1, mode="c"] prns = \
        np.empty(n, dtype=np.uint8)
  for i, a in enumerate(alms):
    a_ = (<Almanac> a).almanac
    memcpy(&al[i], &a_, sizeof(almanac_t))
    prns[i] = (<Almanac> a).prn

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(reciever_ecef, dtype=np.double)

  cdef gps_time_t timestamp_ = timestamp.gps_time

  cdef sdiff_t sdiffs[32]
  almanacs_to_single_diffs(len(alms), &al[0], timestamp_, sdiffs)

  for i, (l,c) in enumerate(numpy_measurements):
    sdiffs[i].pseudorange = c
    sdiffs[i].carrier_phase = l
  
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = \
    np.empty(3, dtype=np.double)

  dgnss_management_c.measure_amb_kf_b(&ref_ecef_[0],
                                      n, &sdiffs[0],
                                      &b[0])
  return b

def measure_b_with_external_ambs(sdiffs, ambs, reciever_ecef): #TODO eventually, want to get reciever_ecef from data
  num_sdiffs = len(sdiffs)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(reciever_ecef, dtype=np.double)
  cdef sdiff_t sdiffs_[32]
  cdef sdiff_t s_
  for (i,sdiff) in enumerate(sdiffs):
    s_ = (<SingleDiff> sdiff).sdiff
    memcpy(&sdiffs_[i], &s_, sizeof(sdiff_t))
  
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ambs_ = \
    np.array(ambs, dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = \
    np.empty(3, dtype=np.double)

  dgnss_management_c.measure_b_with_external_ambs(len(ambs), &ambs_[0],
                                num_sdiffs, &sdiffs_[0],
                                &ref_ecef_[0], &b[0])
  return b

def alm_measure_b_with_external_ambs(alms, GpsTime timestamp,
                                   numpy_measurements,
                                   ambs,
                                   reciever_ecef): #TODO eventually, want to get reciever_ecef from data
  n = len(alms)
  cdef almanac_t al[32]
  cdef almanac_t a_
  cdef np.ndarray[np.uint8_t, ndim=1, mode="c"] prns = \
        np.empty(n, dtype=np.uint8)
  for i, a in enumerate(alms):
    a_ = (<Almanac> a).almanac
    memcpy(&al[i], &a_, sizeof(almanac_t))
    prns[i] = (<Almanac> a).prn

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(reciever_ecef, dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ambs_ = \
    np.array(ambs, dtype=np.double)

  cdef gps_time_t timestamp_ = timestamp.gps_time

  cdef sdiff_t sdiffs[32]
  almanacs_to_single_diffs(len(alms), &al[0], timestamp_, sdiffs)

  for i, (l,c) in enumerate(numpy_measurements):
    sdiffs[i].pseudorange = c
    sdiffs[i].carrier_phase = l
  
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = \
    np.empty(3, dtype=np.double)
  dgnss_management_c.measure_b_with_external_ambs(len(ambs), &ambs_[0],
                                  n, &sdiffs[0],
                                  &ref_ecef_[0], &b[0])
  return b

def measure_iar_b_with_external_ambs(sdiffs, ambs, reciever_ecef): #TODO eventually, want to get reciever_ecef from data
  num_sdiffs = len(sdiffs)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(reciever_ecef, dtype=np.double)
  cdef sdiff_t sdiffs_[32]
  cdef sdiff_t s_
  for (i,sdiff) in enumerate(sdiffs):
    s_ = (<SingleDiff> sdiff).sdiff
    memcpy(&sdiffs_[i], &s_, sizeof(sdiff_t))
  
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ambs_ = \
    np.array(ambs, dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = \
    np.empty(3, dtype=np.double)

  dgnss_management_c.measure_iar_b_with_external_ambs(&ref_ecef_[0],
                                                      num_sdiffs, &sdiffs_[0],
                                                      &ambs_[0],
                                                      &b[0])
  return b

def alm_measure_iar_b_with_external_ambs(alms, GpsTime timestamp,
                                     numpy_measurements,
                                     ambs,
                                     reciever_ecef): #TODO eventually, want to get reciever_ecef from data
  n = len(alms)
  cdef almanac_t al[32]
  cdef almanac_t a_
  cdef np.ndarray[np.uint8_t, ndim=1, mode="c"] prns = \
        np.empty(n, dtype=np.uint8)
  for i, a in enumerate(alms):
    a_ = (<Almanac> a).almanac
    memcpy(&al[i], &a_, sizeof(almanac_t))
    prns[i] = (<Almanac> a).prn

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(reciever_ecef, dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ambs_ = \
    np.array(ambs, dtype=np.double)

  cdef gps_time_t timestamp_ = timestamp.gps_time

  cdef sdiff_t sdiffs[32]
  almanacs_to_single_diffs(len(alms), &al[0], timestamp_, sdiffs)

  for i, (l,c) in enumerate(numpy_measurements):
    sdiffs[i].pseudorange = c
    sdiffs[i].carrier_phase = l
  
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = \
    np.empty(3, dtype=np.double)

  dgnss_management_c.measure_iar_b_with_external_ambs(&ref_ecef_[0],
                                                      n, &sdiffs[0],
                                                      &ambs_[0],
                                                      &b[0])
  return b


def get_float_de_and_phase(sdiffs, ref_ecef):
  num_sdiffs = len(sdiffs)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(ref_ecef, dtype=np.double)
  cdef sdiff_t sdiffs_[32]
  cdef sdiff_t s_
  for (i,sdiff) in enumerate(sdiffs):
    s_ = (<SingleDiff> sdiff).sdiff
    memcpy(&sdiffs_[i], &s_, sizeof(sdiff_t))

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] de = \
    np.empty((32,3), dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] phase = \
    np.empty(32, dtype=np.double)

  m = dgnss_management_c.get_amb_kf_de_and_phase(num_sdiffs, &sdiffs_[0],
                                                 &ref_ecef_[0],
                                                 &de[0,0], &phase[0])

  return de[:m-1], phase[:m-1]

def alms_get_float_de_and_phase(alms, GpsTime timestamp,
                           numpy_measurements,
                           ref_ecef):
  n = len(alms)
  cdef almanac_t al[32]
  cdef almanac_t a_
  cdef np.ndarray[np.uint8_t, ndim=1, mode="c"] prns = \
        np.empty(n, dtype=np.uint8)
  for i, a in enumerate(alms):
    a_ = (<Almanac> a).almanac
    memcpy(&al[i], &a_, sizeof(almanac_t))
    prns[i] = (<Almanac> a).prn

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(ref_ecef, dtype=np.double)

  cdef gps_time_t timestamp_ = timestamp.gps_time

  cdef sdiff_t sdiffs[32]
  almanacs_to_single_diffs(len(alms), &al[0], timestamp_, sdiffs)
  
  for i, (l,c) in enumerate(numpy_measurements):
    sdiffs[i].pseudorange = c
    sdiffs[i].carrier_phase = l

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] de = \
    np.empty((32,3), dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] phase = \
    np.empty(32, dtype=np.double)

  m = dgnss_management_c.get_amb_kf_de_and_phase(n, &sdiffs[0],
                                                 &ref_ecef_[0],
                                                 &de[0,0], &phase[0])

  return de[:m-1], phase[:m-1]

def get_iar_de_and_phase(sdiffs, ref_ecef):
  num_sdiffs = len(sdiffs)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(ref_ecef, dtype=np.double)
  cdef sdiff_t sdiffs_[32]
  cdef sdiff_t s_
  for (i,sdiff) in enumerate(sdiffs):
    s_ = (<SingleDiff> sdiff).sdiff
    memcpy(&sdiffs_[i], &s_, sizeof(sdiff_t))
  
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] de = \
    np.empty((32,3), dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] phase = \
    np.empty(32, dtype=np.double)

  m = dgnss_management_c.get_iar_de_and_phase(num_sdiffs, &sdiffs_[0],
                                              &ref_ecef_[0],
                                              &de[0,0], &phase[0])

  return de[:m-1], phase[:m-1]

def alm_get_iar_de_and_phase(alms, GpsTime timestamp,
                         numpy_measurements,
                         ref_ecef):
  n = len(alms)
  cdef almanac_t al[32]
  cdef almanac_t a_
  cdef np.ndarray[np.uint8_t, ndim=1, mode="c"] prns = \
        np.empty(n, dtype=np.uint8)
  for i, a in enumerate(alms):
    a_ = (<Almanac> a).almanac
    memcpy(&al[i], &a_, sizeof(almanac_t))
    prns[i] = (<Almanac> a).prn

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(ref_ecef, dtype=np.double)

  cdef gps_time_t timestamp_ = timestamp.gps_time

  cdef sdiff_t sdiffs[32]
  almanacs_to_single_diffs(len(alms), &al[0], timestamp_, sdiffs)
  
  for i, (l,c) in enumerate(numpy_measurements):
    sdiffs[i].pseudorange = c
    sdiffs[i].carrier_phase = l

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] de = \
    np.empty((32,3), dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] phase = \
    np.empty(32, dtype=np.double)

  m = dgnss_management_c.get_iar_de_and_phase(n, &sdiffs[0],
                                              &ref_ecef_[0],
                                              &de[0,0], &phase[0])

  return de[:m-1], phase[:m-1]

def dgnss_iar_pool_contains(ambs):
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ambs_ = \
    np.array(ambs, dtype=np.double)
  return dgnss_management_c.dgnss_iar_pool_contains(&ambs_[0]) == 1

def get_amb_kf_mean():
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ambs = \
    np.empty(32, dtype=np.double)
  num_dds = dgnss_management_c.get_amb_kf_mean(&ambs[0])
  return ambs[:num_dds]

def get_amb_kf_cov(num_dds):
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] cov = \
    np.empty((num_dds, num_dds), dtype=np.double)
  num_dds2 = dgnss_management_c.get_amb_kf_cov(&cov[0,0])
  if num_dds == num_dds2:
    return cov
  else:
    raise ValueError("Was given the wrong num_dds.")

def get_amb_kf_cov2():
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] cov = \
    np.empty((32*32, 1), dtype=np.double)
  num_dds = dgnss_management_c.get_amb_kf_cov(&cov[0,0])
  return cov[:num_dds*num_dds,0].reshape((num_dds,num_dds))

def get_amb_kf_prns():
  cdef np.ndarray[np.uint8_t, ndim=1, mode="c"] prns = \
    np.empty(32, dtype=np.uint8)
  num_sats = dgnss_management_c.get_amb_kf_prns(&prns[0])
  return prns[:num_sats]

def get_amb_test_prns():
  cdef np.ndarray[np.uint8_t, ndim=1, mode="c"] prns = \
    np.empty(32, dtype=np.uint8)
  num_sats = dgnss_management_c.get_amb_test_prns(&prns[0])
  return prns[:num_sats]

def dgnss_iar_MLE_ambs():
  cdef np.ndarray[np.int32_t, ndim=1, mode="c"] ambs = \
    np.empty(32, dtype=np.int32)
  num_dds = dgnss_management_c.dgnss_iar_MLE_ambs(&ambs[0])
  return ambs[:num_dds]