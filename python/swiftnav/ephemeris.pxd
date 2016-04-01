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

cdef extern from "libswiftnav/ephemeris.h":

  ctypedef struct ephemeris_kepler_t:
    double tgd
    double crs, crc, cuc, cus, cic, cis
    double dn, m0, ecc, sqrta, omega0, omegadot, w, inc, inc_dot
    double af0, af1, af2
    gps_time_t toc
    u16 iodc
    u8 iode

  ctypedef struct ephemeris_xyz_t:
    double pos[3]
    double vel[3]
    double acc[3]
    double a_gf0
    double a_gf1

  ctypedef struct ephemeris_t:
    gnss_signal_t sid
    gps_time_t toe
    float ura
    u32 fit_interval
    u8 valid
    u8 healthy
    # HACK: Actually an anonymous union in libswiftnat!
    ephemeris_kepler_t kepler
    ephemeris_xyz_t xyz

  s8 calc_sat_state(const ephemeris_t *e, const gps_time_t *t,
                    double pos[3], double vel[3],
                    double *clock_err, double *clock_rate_err)
  s8 calc_sat_az_el(const ephemeris_t *e, const gps_time_t *t,
                    double ref[3], double *az, double *el)
  s8 calc_sat_doppler(const ephemeris_t *e, const gps_time_t *t,
                    double ref[3], double *doppler)
  u8 ephemeris_valid(const ephemeris_t *e, const gps_time_t *t)
  u8 satellite_healthy(const ephemeris_t *e)
  void decode_ephemeris(u32 frame_words[3][8], ephemeris_t *e)
  u8 ephemeris_equal(const ephemeris_t *a, const ephemeris_t *b)
  u8 ephemeris_params_valid(const u8 valid, const u32 fit_interval,
                      const gps_time_t* toe, const gps_time_t *t)


cdef class Ephemeris:
  cdef ephemeris_t _thisptr

cdef mk_ephemeris_array(py_ephs, u8 n_c_ephs, const ephemeris_t **c_ephs)
