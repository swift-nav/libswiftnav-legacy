# Copyright (C) 2012 Swift Navigation Inc.
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
cimport coord_system_c

def wgsllh2ecef(lat, lon, hgt):
  """
  Wraps function :libswiftnav:`wgsllh2ecef`.

  Parameters
  ----------
  lat : float
    Latitude
  lon : float
    Longitude
  hgt : float
    Height

  Returns
  -------
  out : :class:`numpy.ndarray`, shape(3,)
    The array `[x, y, z]`.

  """

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] llh = \
    np.array([lat, lon, hgt], dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ecef = \
    np.empty(3, dtype=np.double)

  coord_system_c.wgsllh2ecef(&llh[0], &ecef[0])

  return ecef

def wgsecef2llh(x, y, z):
  """
  Wraps function :libswiftnav:`wgsecef2llh`.

  Parameters
  ----------
  x : float
    ECEF frame X coordinate
  y : float
    ECEF frame Y coordinate
  z : float
    ECEF frame Z coordinate

  Returns
  -------
  out : :class:`numpy.ndarray`, shape(3,)
    The array `[Latitude, Longitude, Height]`.

  """

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ecef = \
    np.array([x, y, z], dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] llh = \
    np.empty(3, dtype=np.double)
  
  coord_system_c.wgsecef2llh(&ecef[0], &llh[0])
  
  return llh

def wgsecef2ned(ecef, ref_ecef):
  """
  Wraps function :libswiftnav:`wgsecef2ned`.

  Parameters
  ----------
  ecef : (float, float, float)
    The tuple of coordinates, `(x, y, z)`
  ref_ecef : (float, float, float)
    The tuple of coordinates of the reference position, `(x, y, z)`

  Returns
  -------
  out : :class:`numpy.ndarray`, shape(3,)
    The array `[North, East, Down]`.

  """
  
  if len(ecef) != 3 or len(ref_ecef) != 3:
    raise ValueError("ECEF coordinates must have dimension 3.")
  
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ecef_ = \
    np.array(ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(ref_ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ned = \
    np.empty(3, dtype=np.double)
  
  coord_system_c.wgsecef2ned(&ecef_[0], &ref_ecef_[0], &ned[0])
  
  return ned

def wgsecef2ned_d(ecef, ref_ecef):
  """
  Wraps function :libswiftnav:`wgsecef2ned_d`.

  Parameters
  ----------
  ecef : (float, float, float)
    The tuple of coordinates, `(x, y, z)`
  ref_ecef : (float, float, float)
    The tuple of coordinates of the reference position, `(x, y, z)`

  Returns
  -------
  out : :class:`numpy.ndarray`, shape(3,)
    The array `[North, East, Down]`.

  """

  if len(ecef) != 3 or len(ref_ecef) != 3:
    raise ValueError("ECEF coordinates must have dimension 3.")

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ecef_ = \
    np.array(ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(ref_ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ned = \
    np.empty(3, dtype=np.double)

  coord_system_c.wgsecef2ned_d(&ecef_[0], &ref_ecef_[0], &ned[0])

  return ned

def wgsecef2azel(ecef, ref_ecef):
  """
  Wraps function :libswiftnav:`wgsecef2azel`.

  Parameters
  ----------
  ecef : (float, float, float)
    The tuple of coordinates, `(x, y, z)`
  ref_ecef : (float, float, float)
    The tuple of coordinates of the reference position, `(x, y, z)`

  Returns
  -------
  out : (float, float)
    The tuple `(azimuth, elevation)`.

  """

  if len(ecef) != 3 or len(ref_ecef) != 3:
    raise ValueError("ECEF coordinates must have dimension 3.")

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ecef_ = \
    np.array(ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(ref_ecef, dtype=np.double)

  cdef double az
  cdef double el

  coord_system_c.wgsecef2azel(&ecef_[0], &ref_ecef_[0], &az, &el)

  return (az, el)

