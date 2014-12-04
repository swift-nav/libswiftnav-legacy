# Copyright (C) 2012 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from gpstime_c cimport gps_time_t

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

  int calc_sat_pos(double pos[3], double vel[3],
                   double *clock_err, double *clock_rate_err,
                   ephemeris_t *ephemeris,
                   gps_time_t tot)