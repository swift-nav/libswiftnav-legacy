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
from almanac cimport *
from almanac_c cimport *
from gpstime cimport *
from gpstime_c cimport *
from libc.string cimport memcpy
from libc.stdio cimport printf
from sats_management_c cimport *

def dgnss_init(alms, GpsTime timestamp,
               numpy_measurements,
               reciever_ecef, dt, b=None):
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

  dgnss_management_c.dgnss_init(n, &sdiffs[0], &ref_ecef_[0], &b_[0])

def dgnss_update(alms, GpsTime timestamp,
               numpy_measurements,
               reciever_ecef, dt):
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
  dgnss_management_c.dgnss_update(n, &sdiffs[0], &ref_ecef_[0], dt, 1, b_)
  
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = \
    np.empty(3, dtype=np.double)
  memcpy(&b[0],b_, 3*sizeof(double))
  return b
  

def get_dgnss_kf():
  cdef nkf_t *kf = dgnss_management_c.get_dgnss_kf()
  cdef KalmanFilter pykf = KalmanFilter()
  memcpy(&(pykf.kf), kf, sizeof(amb_kf_c.nkf_t))
  return pykf


def get_stupid_state(num):
  cdef np.ndarray[np.int32_t, ndim=1, mode="c"] N = \
    np.empty(num, dtype=np.int32)
  cdef s32 * c_N = dgnss_management_c.get_stupid_filter_ints()
  memcpy(&N[0], c_N, num * sizeof(s32))
  return N

def get_sats_management():
  cdef sats_management_t * sats_man = dgnss_management_c.get_sats_management()
  cdef np.ndarray[np.uint8_t, ndim=1, mode="c"] prns = \
    np.empty(sats_man.num_sats, dtype=np.uint8)
  memcpy(&prns[0], sats_man.prns, sats_man.num_sats * sizeof(u8))
  return sats_man.num_sats, prns








