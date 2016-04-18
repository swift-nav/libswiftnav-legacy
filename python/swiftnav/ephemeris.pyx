# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

""" Ephemeris utilities

"""

cimport numpy as np
from cpython.object cimport Py_EQ
from fmt_utils import fmt_repr
from time cimport *
from time import GpsTime
from libc.string cimport memcpy, memset
import numpy as np


cdef class Ephemeris:

  def __init__(self, **kwargs):
    memset(&self._thisptr, 0, sizeof(ephemeris_t))
    self._thisptr.sid = kwargs.pop('sid')
    self._thisptr.toe = kwargs.pop('toe')
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

  def calc_sat_state(self, GpsTime t):
    """
    Wraps the function :libswiftnav:`calc_sat_state`.

    Parameters
    ----------
    t : GpsTime
      The GPS time at which to calculate the satellite position.

    Returns
    -------
    out : (:class:`numpy.ndarray`, :class:`numpy.ndarray`)
      The tuple (position, velocity) in meters and meters/sec in the ECEF
      coordinate system.

    """
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] pos = np.array([0,0,0], dtype=np.double)
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] vel = np.array([0,0,0], dtype=np.double)
    cdef double clock_err, clock_rate_err
    calc_sat_state(&self._thisptr, &t._thisptr, &pos[0], &vel[0], &clock_err, &clock_rate_err)
    return (pos, vel, clock_err, clock_rate_err)

  def calc_az_el(self, GpsTime t, ref):
    """
    Wraps the function :libswiftnav:`calc_sat_az_el`.

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
    calc_sat_az_el(&self._thisptr, &t._thisptr, &ref_[0], &az, &el)
    return (az, el)

  def calc_doppler(self, GpsTime t, ref):
    """
    Wraps the function :libswiftnav:`calc_sat_doppler`.

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
    calc_sat_doppler(&self._thisptr, &t._thisptr, &ref_[0], &doppler)
    return doppler

  def is_valid(self, GpsTime t):
    return ephemeris_valid(&self._thisptr, &t._thisptr)

  @staticmethod
  def is_params_valid(u8 valid, u32 fit_interval, GpsTime toe, GpsTime t):
    return ephemeris_params_valid(valid, fit_interval, &toe._thisptr, &t._thisptr)

  def is_healthy(self):
    return satellite_healthy(&self._thisptr)

  def __rich_cmp__(self, eph2, op):
    if op != Py_EQ:
      return False
    cdef ephemeris_t eph = eph2._thisptr
    return ephemeris_equal(&self._thisptr, &eph)

# TODO (Buro):  Fix this
# def decode_eph(np.ndarray[np.uint32_t, ndim=2] frame_words):
#   cdef ephemeris_t tmp
#   decode_ephemeris(&frame_words[0], &tmp)
#   return Ephemeris(tmp.tgd, tmp.crs, tmp.crc, tmp.cuc, tmp.cus, tmp.cic, tmp.cis,
#                    tmp.dn, tmp.m0, tmp.ecc, tmp.sqrta, tmp.omega0, tmp.omegadot,
#                    tmp.w, tmp.inc, tmp.inc_dot, tmp.af0, tmp.af1, tmp.af2,
#                    tmp.toe, tmp.toc,
#                    tmp.valid, tmp.healthy, tmp.prn)

cdef mk_ephemeris_array(py_ephemerides, u8 n_c_ephemerides,
                        const ephemeris_t **c_ephemerides):
  """Given an array of python Ephemerides, copies their contents to an
  array of ephemeris_t's.

  """
  for (i, ephemeris) in enumerate(py_ephemerides):
    c_ephemerides[i] = &((<Ephemeris ?>ephemeris)._thisptr)
