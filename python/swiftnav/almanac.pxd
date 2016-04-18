# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from time cimport gps_time_t
from signal cimport gnss_signal_t

cdef extern from "libswiftnav/almanac.h":
  ctypedef struct almanac_kepler_t:
    double m0, ecc, sqrta, omega0, omegadot, w, inc
    double af0, af1

  ctypedef struct almanac_xyz_t:
    double pos[3]
    double vel[3]
    double acc[3]

  ctypedef struct almanac_t:
    gnss_signal_t sid
    gps_time_t toa
    float ura
    u32 fit_interval
    u8 valid
    u8 healthy
    # HACK: Actually an anonymous union in libswiftnat!
    almanac_kepler_t kepler
    almanac_xyz_t xyz

  void calc_sat_state_almanac(const almanac_t *a, const gps_time_t *t,
                              double pos[3], double vel[3],
                              double *clock_err, double *clock_rate_err)
  void calc_sat_az_el_almanac(const almanac_t *a, const gps_time_t *t,
                              double ref[3], double *az, double *el)
  double calc_sat_doppler_almanac(const almanac_t *a, const gps_time_t *t,
                                  double ref[3], double *doppler)
  u8 almanac_valid(const almanac_t *a, const gps_time_t *t)
  u8 satellite_healthy_almanac(const almanac_t *a)
  u8 almanac_equal(const almanac_t *a, const almanac_t *b)

cdef class Almanac:
  cdef almanac_t _thisptr
