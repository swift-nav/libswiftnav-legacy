# Copyright (C) 2013 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

# cython: embedsignature=True

import numpy as np
cimport numpy as np
cimport almanac_c

cdef class Almanac:
  cdef almanac_c.almanac_t almanac

  def __init__(self, ecc, toa, inc, rora, a, raaw,
               argp, ma, af0, af1, week, prn, healthy):
    self.ecc = ecc
    self.toa = toa
    self.inc = inc
    self.rora = rora
    self.a = a
    self.raaw = raaw
    self.argp = argp
    self.ma = ma
    self.af0 = af0
    self.af1 = af1
    self.week = week
    self.prn = prn
    self.healthy = healthy

  def calc_state(self, t, week=None):
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] pos = \
      np.empty(3, dtype=np.double)
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] vel = \
      np.empty(3, dtype=np.double)

    if week is None:
      week = -1

    almanac_c.calc_sat_state_almanac(&self.almanac, t, week, &pos[0], &vel[0])

    return (pos, vel)

  def calc_az_el(self, t, ref, week=None):
    if len(ref) != 3:
      raise ValueError("ECEF coordinates must have dimension 3.")

    cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ = \
      np.array(ref, dtype=np.double)

    cdef double az, el

    if week is None:
      week = -1

    almanac_c.calc_sat_az_el_almanac(&self.almanac, t, week,
                                     &ref_[0], &az, &el)

    return (az, el)

  def calc_doppler(self, t, ref, week=None):
    if len(ref) != 3:
      raise ValueError("ECEF coordinates must have dimension 3.")

    cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ = \
      np.array(ref, dtype=np.double)

    cdef double az, el

    if week is None:
      week = -1

    return almanac_c.calc_sat_doppler_almanac(&self.almanac, t, week,
                                              &ref_[0])

  def __repr__(self):
    return "<Almanac PRN %02d, Week %d, ToA %.1f>" % (self.prn, self.week,
                                                      self.toa)

  property ecc:
    def __get__(self):
      return self.almanac.ecc
    def __set__(self, ecc):
      self.almanac.ecc = ecc

  property toa:
    def __get__(self):
      return self.almanac.toa
    def __set__(self, toa):
      self.almanac.toa = toa

  property inc:
    def __get__(self):
      return self.almanac.inc
    def __set__(self, inc):
      self.almanac.inc = inc

  property rora:
    def __get__(self):
      return self.almanac.rora
    def __set__(self, rora):
      self.almanac.rora = rora

  property a:
    def __get__(self):
      return self.almanac.a
    def __set__(self, a):
      self.almanac.a = a

  property raaw:
    def __get__(self):
      return self.almanac.raaw
    def __set__(self, raaw):
      self.almanac.raaw = raaw

  property argp:
    def __get__(self):
      return self.almanac.argp
    def __set__(self, argp):
      self.almanac.argp = argp

  property ma:
    def __get__(self):
      return self.almanac.ma
    def __set__(self, ma):
      self.almanac.ma = ma

  property af0:
    def __get__(self):
      return self.almanac.af0
    def __set__(self, af0):
      self.almanac.af0 = af0

  property af1:
    def __get__(self):
      return self.almanac.af1
    def __set__(self, af1):
      self.almanac.af1 = af1

  property week:
    def __get__(self):
      return self.almanac.week
    def __set__(self, week):
      self.almanac.week = week

  property prn:
    def __get__(self):
      return self.almanac.prn
    def __set__(self, prn):
      self.almanac.prn = prn

  property healthy:
    def __get__(self):
      return (self.almanac.healthy > 0)
    def __set__(self, healthy):
      self.almanac.healthy = 1 if healthy else 0

