# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

cimport numpy as np
from cpython.object cimport Py_EQ
from gpstime cimport *
from gpstime import GpsTime
from libc.string cimport memcpy
import numpy as np

cdef class Ephemeris:

  def __init__(self,
               tgd,
               crs, crc, cuc, cus, cic, cis,
               dn, m0, ecc, sqrta, omega0, omegadot, w, inc, inc_dot,
               af0, af1, af2,
               toe, toc,
               valid,
               healthy,
               prn):
    self._thisptr.tgd = tgd
    self._thisptr.crs = crs
    self._thisptr.crc = crc
    self._thisptr.cuc = cuc
    self._thisptr.cus = cus
    self._thisptr.cic = cic
    self._thisptr.cis = cis
    self._thisptr.dn = dn
    self._thisptr.m0 = m0
    self._thisptr.ecc = ecc
    self._thisptr.sqrta = sqrta
    self._thisptr.omega0 = omega0
    self._thisptr.omegadot = omegadot
    self._thisptr.w   = w
    self._thisptr.inc = inc
    self._thisptr.inc_dot = inc_dot
    self._thisptr.af0 = af0
    self._thisptr.af1 = af1
    self._thisptr.af2 = af2
    self._thisptr.toe = toe
    self._thisptr.toc = toc
    self._thisptr.valid   = int(valid)
    self._thisptr.healthy = healthy
    self._thisptr.prn = int(prn)

  def __repr__(self):
    return "<Ephemeris prn=%d, toe=%d>" % (self._thisptr.prn, self._thisptr.toe)

  def calc_sat_state(self, GpsTime time):
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] pos = np.array([0,0,0], dtype=np.double)
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] vel = np.array([0,0,0], dtype=np.double)
    cdef double clock_err, clock_rate_err
    calc_sat_state(&(self._thisptr), time._thisptr, &pos[0], &vel[0], &clock_err, &clock_rate_err)
    return (pos, vel, clock_err, clock_rate_err)

  def ephemeris_good(self, GpsTime time):
    return ephemeris_good(&(self._thisptr), time._thisptr)

  def __rich_cmp__(self, eph2, op):
    if op != Py_EQ:
      return False
    cdef ephemeris_t eph = eph2._thisptr
    return ephemeris_equal(&self._thisptr, &eph)

# TODO (Buro): Fix this!
# def decode_eph(np.ndarray[np.uint32_t, ndim=2] frame_words):
#   cdef ephemeris_t tmp
#   decode_ephemeris(frame_words, &tmp)
#   return Ephemeris(tmp.tgd, tmp.crs, tmp.crc, tmp.cuc, tmp.cus, tmp.cic, tmp.cis,
#                    tmp.dn, tmp.m0, tmp.ecc, tmp.sqrta, tmp.omega0, tmp.omegadot,
#                    tmp.w, tmp.inc, tmp.inc_dot, tmp.af0, tmp.af1, tmp.af2,
#                    tmp.toe, tmp.toc,
#                    tmp.valid, tmp.healthy, tmp.prn)
