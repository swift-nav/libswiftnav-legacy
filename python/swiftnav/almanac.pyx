# Copyright (C) 2013 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

# cython: embedsignature=True

from common cimport *
cimport numpy as np
import numpy as np

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

  def __init__(self, ecc, toa, inc, rora, a, raaw,
               argp, ma, af0, af1, week, prn, healthy):
    self._thisptr.ecc = ecc
    self._thisptr.ecc = ecc
    self._thisptr.toa = toa
    self._thisptr.inc = inc
    self._thisptr.rora = rora
    self._thisptr.a = a
    self._thisptr.raaw = raaw
    self._thisptr.argp = argp
    self._thisptr.ma = ma
    self._thisptr.af0 = af0
    self._thisptr.af1 = af1
    self._thisptr.week = week
    self._thisptr.prn = prn
    self._thisptr.healthy = healthy

  def calc_state(self, t, week=-1):
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
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] pos = np.empty(3, dtype=np.double)
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] vel = np.empty(3, dtype=np.double)
    calc_sat_state_almanac(&self._thisptr, t, week, &pos[0], &vel[0])
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
    assert len(ref) != 3, "ECEF coordinates must have dimension 3."
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ = np.array(ref, dtype=np.double)
    cdef double az, el
    calc_sat_az_el_almanac(&self._thisptr, t, week, &ref_[0], &az, &el)
    return (az, el)

  def calc_doppler(self, t, ref, week=-1):
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
    assert len(ref) != 3, "ECEF coordinates must have dimension 3."
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ = np.array(ref, dtype=np.double)
    return calc_sat_doppler_almanac(&self._thisptr, t, week, &ref_[0])

  def __repr__(self):
    return "<Almanac PRN %02d, Week %d, ToA %.1f>" \
        % (self._thisptr.prn, self._thisptr.week, self._thisptr.toa)
