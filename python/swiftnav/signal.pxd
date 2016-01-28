# Copyright (C) 2016 Swift Navigation Inc.
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
    CONSTELLATION_GPS
    CONSTELLATION_SBAS
    CONSTELLATION_COUNT
    CODE_GPS_L1CA
    CODE_GPS_L2CL
    CODE_GPS_L2CM
    CODE_SBAS_L1CA
    CODE_COUNT
    INDEX_ALL_CODES
    INDEX_SATS_IN_CONS

  ctypedef struct gnss_signal_t:
    u16 sat
    u8 code

  int sid_compare(const gnss_signal_t a, const gnss_signal_t b)
  int cmp_sid_sid(const void *a, const void *b)
  bool sid_is_equal(const gnss_signal_t a, const gnss_signal_t b)

  int sid_to_string(char *s, int n, gnss_signal_t sid)
  bool sid_valid(gnss_signal_t sid)

  gnss_signal_t construct_sid(u8 code, u16 sat);
  gnss_signal_t sid_from_index(u16 i)
  u16 sid_to_index(gnss_signal_t sid, u8 it)

cdef class GNSSSignal:
  cdef gnss_signal_t _thisptr

cdef mk_signal_array(py_signals, u8 n_c_signals, gnss_signal_t *c_signals)
cdef read_signal_array(u8 n_c_signals, gnss_signal_t *c_signals)
