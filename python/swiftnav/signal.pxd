# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from libcpp cimport bool

cdef extern from "libswiftnav/signal.h":
  enum:
    GPS_L1_SATS
    WAAS_SATS
    EGNOS_SATS
    GAGAN_SATS
    MSAS_SATS
    SDCM_SATS
    SBAS_SATS
    ALL_SATS
    CONSTELLATION_GPS
    CONSTELLATION_SBAS
    BAND_L1

  ctypedef struct gnss_signal_t:
    u16 sat
    u8 band
    u8 constellation

  int sid_compare(const gnss_signal_t a, const gnss_signal_t b)
  int cmp_sid_sid(const void *a, const void *b)
  bool sid_is_equal(const gnss_signal_t a, const gnss_signal_t b)

  int sid_to_string(char *s, int n, gnss_signal_t sid)
  bool sid_valid(gnss_signal_t sid)
  gnss_signal_t sid_from_index(u32 i)
  u32 sid_to_index(gnss_signal_t sid)

cdef class GNSSSignal:
  cdef gnss_signal_t _thisptr

cdef mk_signal_array(py_signals, u8 n_c_signals, gnss_signal_t *c_signals)
cdef read_signal_array(u8 n_c_signals, gnss_signal_t *c_signals)
