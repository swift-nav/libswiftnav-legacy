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
  """
  Wraps the :libswiftnav:`almanac_t` structure and associated functions for
  performing calculations with almanacs.

  Parameters
  ----------
  ecc : float
    Eccentricity in radians.
  toa : float
    Time of Applicability in seconds since Sunday.
  inc : float
    Inclination in radians.
  rora : float
    Rate of Right Ascension in radians/sec.
  a : float
    Semi-major axis in meters.
  raaw : float
    Right Ascension at Week in radians.
  argp : float
    Argument of Perigee in radians.
  ma : float
    Mean Anomaly at Time of Applicability in radians.
  af0 : float
    0-order clock correction in seconds.
  af1 : float
    1-order clock correction in seconds/second.
  week : int
    GPS week number, modulo 1024.
  prn : int
    PRN number of the satellite.
  healthy : bool
    Satellite health status.

  """
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
    """
    Wraps the function :libswiftnav:`calc_sat_state_almanac`.

    Parameters
    ----------
    t : float
      The GPS time at which to calculate the azimuth and elevation.
    week : int, optional
      The GPS week number modulo 1024. If `None` then it is assumed that `t` is
      within one half week of the almanac time of applicability.

    Returns
    -------
    out : (:class:`numpy.ndarray`, :class:`numpy.ndarray`)
      The tuple (position, velocity) in meters and meters/sec in the ECEF
      coordinate system.

    """
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] pos = \
      np.empty(3, dtype=np.double)
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] vel = \
      np.empty(3, dtype=np.double)

    if week is None:
      week = -1

    almanac_c.calc_sat_state_almanac(&self.almanac, t, week, &pos[0], &vel[0])

    return (pos, vel)

  def calc_az_el(self, t, ref, week=None):
    """
    Wraps the function :libswiftnav:`calc_sat_az_el_almanac`.

    Parameters
    ----------
    t : float
      The GPS time at which to calculate the azimuth and elevation.
    ref : (float, float, float)
      The tuple of coordinates of the reference position, `(x, y, z)`
    week : int, optional
      The GPS week number modulo 1024. If `None` then it is assumed that `t` is
      within one half week of the almanac time of applicability.

    Returns
    -------
    out : (float, float)
      The tuple (azimuth, elevation) in radians.

    """
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
    """
    Wraps the function :libswiftnav:`calc_sat_doppler_almanac`.

    Parameters
    ----------
    t : float
      The GPS time at which to calculate the Doppler shift.
    ref : (float, float, float)
      The tuple of coordinates of the reference position, `(x, y, z)`
    week : int, optional
      The GPS week number modulo 1024. If `None` then it is assumed that `t` is
      within one half week of the almanac time of applicability.

    Returns
    -------
    out : float
      The Doppler shift in Hz.

    """
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
    return "<Almanac PRN %02d, Week %d, ToA %.1f>" % \
        (self.prn, self.week, self.toa)

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

