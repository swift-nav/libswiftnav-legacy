# Copyright (C) 2013 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

# cython: embedsignature=True

"""Almanac

Functions and calculations related to the GPS almanac.

Note: All positions are referenced to the WGS84 coordinate system.
"""

from common cimport *
from time cimport *
from time import GpsTime
from fmt_utils import fmt_repr
from libc.string cimport memset
cimport numpy as np
import numpy as np

cdef class Almanac:

  def __init__(self, **kwargs):
    memset(&self._thisptr, 0, sizeof(almanac_t))
    self._thisptr.sid = kwargs.pop('sid')
    self._thisptr.toa = kwargs.pop('toa')
    self._thisptr.ura = kwargs.pop('ura')
    self._thisptr.fit_interval = kwargs.pop('fit_interval')
    self._thisptr.valid = kwargs.pop('valid')
    self._thisptr.healthy = kwargs.pop('healthy')
    if 'kepler' in kwargs:
      self._thisptr.kepler = kwargs.pop('kepler')
    elif 'xyz' in kwargs:
      self._thisptr.xyz = kwargs.pop('xyz')

  def __getattr__(self, k):
    return self._thisptr.get(k)

  def __repr__(self):
    return fmt_repr(self)

  def to_dict(self):
    return self._thisptr

  def from_dict(self, d):
    self._thisptr = d

  def calc_state(self, GpsTime t):
    """
    Wraps the function :libswiftnav:`calc_sat_state_almanac`.

    Parameters
    ----------
    t : GpsTime
      The GPS time at which to calculate the azimuth and elevation.

    Returns
    -------
    out : (:class:`numpy.ndarray`, :class:`numpy.ndarray`)
      The tuple (position, velocity) in meters and meters/sec in the ECEF
      coordinate system.

    """
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] pos = np.array([0,0,0], dtype=np.double)
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] vel = np.array([0,0,0], dtype=np.double)
    cdef double clock_err, clock_rate_err
    calc_sat_state_almanac(&self._thisptr, &t._thisptr, &pos[0], &vel[0], &clock_err, &clock_rate_err)
    return (pos, vel, clock_err, clock_rate_err)

  def calc_az_el(self, GpsTime t, ref):
    """
    Wraps the function :libswiftnav:`calc_sat_az_el_almanac`.

    Parameters
    ----------
    t : GpsTime
      The GPS time at which to calculate the azimuth and elevation.
    ref : (float, float, float)
      The tuple of coordinates of the reference position, `(x, y, z)`

    Returns
    -------
    out : (float, float)
      The tuple (azimuth, elevation) in radians.

    """
    assert len(ref) == 3, "ECEF coordinates must have dimension 3."
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ = np.array(ref, dtype=np.double)
    cdef double az, el
    calc_sat_az_el_almanac(&self._thisptr, &t._thisptr, &ref_[0], &az, &el)
    return (az, el)

  def calc_doppler(self, GpsTime t, ref):
    """
    Wraps the function :libswiftnav:`calc_sat_doppler_almanac`.

    Parameters
    ----------
    t : GpsTime
      The GPS time at which to calculate the Doppler shift.
    ref : (float, float, float)
      The tuple of coordinates of the reference position, `(x, y, z)`

    Returns
    -------
    out : float
      The Doppler shift in Hz.

    """
    assert len(ref) == 3, "ECEF coordinates must have dimension 3."
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ = np.array(ref, dtype=np.double)
    cdef double doppler
    calc_sat_doppler_almanac(&self._thisptr, &t._thisptr, &ref_[0], &doppler)
    return doppler
