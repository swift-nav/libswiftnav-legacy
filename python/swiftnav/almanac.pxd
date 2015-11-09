# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *

cdef extern from "libswiftnav/almanac.h":
  ctypedef struct almanac_t:
    double ecc
    double toa
    double inc
    double rora
    double a
    double raaw
    double argp
    double ma
    double af0
    double af1
    u16 week
    u8 prn
    u8 healthy

  void calc_sat_state_almanac(almanac_t* alm, double t, s16 week, double pos[3], double vel[3])
  void calc_sat_az_el_almanac(almanac_t* alm, double t, s16 week, double ref[3], double* az, double* el)
  double calc_sat_doppler_almanac(almanac_t* alm, double t, s16 week, double ref[3])

cdef class Almanac:
  cdef almanac_t _thisptr
  cdef public double ecc, toa, inc, rora, a, raaw, argp, ma, af0, af1
  cdef public u16 week
  cdef public u8 prn, healthy
