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
cimport ephemeris_c
from gpstime cimport *
from gpstime_c cimport *
from gpstime import GpsTime
from libc.string cimport memcpy

cdef class Ephemeris:
  cdef ephemeris_c.ephemeris_t ephemeris
  def __init__(self,
               tgd,
               crs, crc, cuc, cus, cic, cis,
               dn, m0, ecc, sqrta, omega0, omegadot, w, inc, inc_dot,
               af0, af1, af2,
               toe, toc,
               valid,
               healthy,
               prn):
    self.tgd = tgd
    self.crs = crs
    self.crc = crc
    self.cuc = cuc
    self.cus = cus
    self.cic = cic
    self.cis = cis
    self.dn  = dn
    self.m0  = m0
    self.ecc   = ecc
    self.sqrta = sqrta
    self.omega0   = omega0
    self.omegadot = omegadot
    self.w   = w
    self.inc = inc
    self.inc_dot = inc_dot
    self.af0 = af0
    self.af1 = af1
    self.af2 = af2
    self.toe = toe
    self.toc = toc
    self.valid   = int(valid)
    self.healthy = healthy
    self.prn = int(prn)

  def __repr__(self):
    return "<Ephemeris prn=" + str(self.prn) + ", toe=" + str(self.toe) + ">"

  property tgd:
    def __get__(self):
        return self.ephemeris.tgd
    def __set__(self, tgd):
        self.ephemeris.tgd = tgd

  property crs:
    def __get__(self):
        return self.ephemeris.crs
    def __set__(self, crs):
        self.ephemeris.crs = crs

  property crc:
    def __get__(self):
        return self.ephemeris.crc
    def __set__(self, crc):
        self.ephemeris.crc = crc

  property cuc:
    def __get__(self):
        return self.ephemeris.cuc
    def __set__(self, cuc):
        self.ephemeris.cuc = cuc

  property cus:
    def __get__(self):
        return self.ephemeris.cus
    def __set__(self, cus):
        self.ephemeris.cus = cus

  property cic:
    def __get__(self):
        return self.ephemeris.cic
    def __set__(self, cic):
        self.ephemeris.cic = cic

  property cis:
    def __get__(self):
        return self.ephemeris.cis
    def __set__(self, cis):
        self.ephemeris.cis = cis

  property dn:
    def __get__(self):
        return self.ephemeris.dn
    def __set__(self, dn):
        self.ephemeris.dn = dn

  property m0:
    def __get__(self):
        return self.ephemeris.m0
    def __set__(self, m0):
        self.ephemeris.m0 = m0

  property ecc:
    def __get__(self):
        return self.ephemeris.ecc
    def __set__(self, ecc):
        self.ephemeris.ecc = ecc

  property sqrta:
    def __get__(self):
        return self.ephemeris.sqrta
    def __set__(self, sqrta):
        self.ephemeris.sqrta = sqrta

  property omega0:
    def __get__(self):
        return self.ephemeris.omega0
    def __set__(self, omega0):
        self.ephemeris.omega0 = omega0

  property omegadot:
    def __get__(self):
        return self.ephemeris.omegadot
    def __set__(self, omegadot):
        self.ephemeris.omegadot = omegadot

  property w:
    def __get__(self):
        return self.ephemeris.w
    def __set__(self, w):
        self.ephemeris.w = w

  property inc:
    def __get__(self):
        return self.ephemeris.inc
    def __set__(self, inc):
        self.ephemeris.inc = inc

  property inc_dot:
    def __get__(self):
        return self.ephemeris.inc_dot
    def __set__(self, inc_dot):
        self.ephemeris.inc_dot = inc_dot

  property af0:
    def __get__(self):
        return self.ephemeris.af0
    def __set__(self, af0):
        self.ephemeris.af0 = af0

  property af1:
    def __get__(self):
        return self.ephemeris.af1
    def __set__(self, af1):
        self.ephemeris.af1 = af1

  property af2:
    def __get__(self):
        return self.ephemeris.af2
    def __set__(self, af2):
        self.ephemeris.af2 = af2

  property toe:
    def __get__(self):
        return GpsTime(self.ephemeris.toe.wn, self.ephemeris.toe.tow)
    def __set__(self, toe):
        self.ephemeris.toe.wn  = toe.wn
        self.ephemeris.toe.tow = toe.tow

  property toc:
    def __get__(self):
        return GpsTime(self.ephemeris.toc.wn, self.ephemeris.toc.tow)
    def __set__(self, toc):
        self.ephemeris.toc.wn  = toc.wn
        self.ephemeris.toc.tow = toc.tow

  property valid:
    def __get__(self):
        return self.ephemeris.valid
    def __set__(self, valid):
        self.ephemeris.valid = valid

  property healthy:
    def __get__(self):
        return self.ephemeris.healthy
    def __set__(self, healthy):
        self.ephemeris.healthy = healthy

  property prn:
    def __get__(self):
        return self.ephemeris.prn
    def __set__(self, prn):
        self.ephemeris.prn = prn

def calc_sat_pos(Ephemeris eph, GpsTime time):
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] pos_ = \
    np.array([0,0,0], dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] vel_ = \
    np.array([0,0,0], dtype=np.double)

  cdef double clock_err
  cdef double clock_rate_err

  cdef ephemeris_c.ephemeris_t eph_ = eph.ephemeris
  cdef gpstime_c.gps_time_t t_      = time.gps_time

  ephemeris_c.calc_sat_pos(&pos_[0], &vel_[0],
                 &clock_err, &clock_rate_err,
                 &eph_,
                 t_)

  return pos_, vel_, clock_err, clock_rate_err

# def measure_b_with_external_ambs(alms, GpsTime timestamp,
#                                    numpy_measurements,
#                                    ambs,
#                                    reciever_ecef): #TODO eventually, want to get reciever_ecef from data
#   n = len(alms)
#   cdef almanac_t al[32]
#   cdef almanac_t a_
#   cdef np.ndarray[np.uint8_t, ndim=1, mode="c"] prns = \
#         np.empty(n, dtype=np.uint8)
#   for i, a in enumerate(alms):
#     a_ = (<Almanac> a).almanac
#     memcpy(&al[i], &a_, sizeof(almanac_t))
#     prns[i] = (<Almanac> a).prn

#   cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
#     np.array(reciever_ecef, dtype=np.double)

#   cdef np.ndarray[np.double_t, ndim=1, mode="c"] ambs_ = \
#     np.array(ambs, dtype=np.double)

#   cdef gps_time_t timestamp_ = timestamp.gps_time

#   cdef sdiff_t sdiffs[32]
#   almanacs_to_single_diffs(len(alms), &al[0], timestamp_, sdiffs)

#   for i, (l,c) in enumerate(numpy_measurements):
#     sdiffs[i].pseudorange = c
#     sdiffs[i].carrier_phase = l
  
#   cdef np.ndarray[np.double_t, ndim=1, mode="c"] b = \
#     np.empty(3, dtype=np.double)

#   dgnss_management_c.measure_b_with_external_ambs(&ref_ecef_[0],
#                                                   n, &sdiffs[0],
#                                                   &ambs_[0],
#                                                   &b[0])
#   return b
