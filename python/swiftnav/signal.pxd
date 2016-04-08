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
  cdef enum constellation:
    CONSTELLATION_INVALID = -1
    CONSTELLATION_GPS
    CONSTELLATION_SBAS
    CONSTELLATION_GLO
    CONSTELLATION_COUNT

  cdef enum code:
    CODE_INVALID = -1
    CODE_GPS_L1CA
    CODE_GPS_L2CM
    CODE_SBAS_L1CA
    CODE_GLO_L1CA
    CODE_GLO_L2CA
    CODE_COUNT

  enum:
    NUM_SATS_GPS
    NUM_SATS_SBAS
    NUM_SATS_GLO
    NUM_SATS
    NUM_CODES_GPS
    NUM_CODES_SBAS
    NUM_CODES_GLO
    NUM_SIGNALS_GPS_L1CA
    NUM_SIGNALS_GPS_L2CM
    NUM_SIGNALS_SBAS_L1CA
    NUM_SIGNALS_GLO_L1CA
    NUM_SIGNALS_GLO_L2CA
    NUM_SIGNALS_GPS
    NUM_SIGNALS_SBAS
    NUM_SIGNALS_GLO
    NUM_SIGNALS
    GPS_FIRST_PRN
    SBAS_FIRST_PRN
    GLO_FIRST_PRN
    SID_STR_LEN_MAX

  ctypedef struct gnss_signal_t:
    u16 sat
    code code

  int sid_compare(const gnss_signal_t a, const gnss_signal_t b)
  int cmp_sid_sid(const void *a, const void *b)
  bool sid_is_equal(const gnss_signal_t a, const gnss_signal_t b)

  gnss_signal_t construct_sid(code code, u16 sat);
  int sid_to_string(char *s, int n, gnss_signal_t sid)
  bool sid_valid(gnss_signal_t sid)
  bool code_valid(code code)
  bool constellation_valid(constellation constellation)

  gnss_signal_t sid_from_code_index(code code, u16 code_index);
  u16 sid_to_code_index(gnss_signal_t sid);
  constellation sid_to_constellation(gnss_signal_t sid);
  constellation code_to_constellation(code code);

cdef class GNSSSignal:
  cdef gnss_signal_t _thisptr

cdef mk_signal_array(py_signals, u8 n_c_signals, gnss_signal_t *c_signals)
cdef read_signal_array(u8 n_c_signals, gnss_signal_t *c_signals)
