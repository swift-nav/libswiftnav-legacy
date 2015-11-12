# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from gpstime cimport gps_time_t

cdef extern from "libswiftnav/ephemeris.h":
  ctypedef struct ephemeris_t:
   double tgd
   double crs, crc, cuc, cus, cic, cis
   double dn, m0, ecc, sqrta, omega0, omegadot, w, inc, inc_dot
   double af0, af1, af2
   gps_time_t toe, toc
   u8 valid
   u8 healthy
   u8 prn

  s8 calc_sat_state(const ephemeris_t *ephemeris, gps_time_t t,
                    double pos[3], double vel[3],
                    double *clock_err, double *clock_rate_err)
  u8 ephemeris_good(ephemeris_t *eph, gps_time_t t)
  void decode_ephemeris(u32 frame_words[3][8], ephemeris_t *e)
  u8 ephemeris_equal(ephemeris_t *a, ephemeris_t *b)

cdef class Ephemeris:
  cdef ephemeris_t _thisptr
  cdef public double tgd
  cdef public double crs, crc, cuc, cus, cic, cis
  cdef public double dn, m0, ecc, sqrta, omega0, omegadot, w, inc, inc_dot
  cdef public double af0, af1, af2
  cdef public gps_time_t toe, toc
  cdef public u8 valid
  cdef public u8 healthy
  cdef public u8 prn
