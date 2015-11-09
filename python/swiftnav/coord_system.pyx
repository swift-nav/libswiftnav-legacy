# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

# cython: embedsignature=True

cimport numpy as np
import numpy as np

cdef extern from "libswiftnav/coord_system.h":
  float WGS84_A
  float WGS84_IF
  float WGS84_F
  float WGS84_B
  float WGS84_E

  void llhrad2deg(const double llh_rad[3], double llh_deg[3])
  void llhdeg2rad(const double llh_deg[3], double llh_rad[3])
  void wgsllh2ecef(const double llh[3], double ecef[3])
  void wgsecef2llh(const double ecef[3], double llh[3])
  void wgsecef2ned(const double ecef[3], const double ref_ecef[3], double ned[3])
  void wgsecef2ned_d(const double ecef[3], const double ref_ecef[3], double ned[3])
  void wgsned2ecef(const double ned[3], const double ref_ecef[3], double ecef[3])
  void wgsned2ecef_d(const double ned[3], const double ref_ecef[3], double ecef[3])
  void wgsecef2azel(const double ecef[3], const double ref_ecef[3], double* azimuth, double* elevation)
  void ecef2ned_matrix(const double ref_ecef[3], double M[3][3])

def llhrad2deg(llh_rad):
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] llh_rad_ = np.array(llh_rad, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] llh_deg_ = np.empty(3, dtype=np.double)
  llhrad2deg(llh_rad_[0], llh_deg_)
  return llh_deg_

def llhdeg2rad(llh_deg):
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] llh_deg_ = np.array(llh_deg, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] llh_rad_ = np.empty(3, dtype=np.double)
  llhdeg2rad(llh_deg_[0], llh_rad_)
  return llh_rad_

def wgsllh2ecef_(lat, lon, hgt):
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

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] llh = np.array([lat, lon, hgt], dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ecef = np.empty(3, dtype=np.double)
  wgsllh2ecef(&llh[0], &ecef[0])
  return ecef

def wgsecef2llh_(x, y, z):
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
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ecef = np.array([x, y, z], dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] llh = np.empty(3, dtype=np.double)
  wgsecef2llh(&ecef[0], &llh[0])
  return llh

def wgsecef2ned_(ecef, ref_ecef):
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

  assert len(ecef) == 3 and len(ref_ecef) == 3, "ECEF coordinates must have dimension 3."
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ecef_ = np.array(ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(ref_ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ned = np.empty(3, dtype=np.double)
  wgsecef2ned(&ecef_[0], &ref_ecef_[0], &ned[0])
  return ned

def wgsecef2ned_d_(ecef, ref_ecef):
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
  assert len(ecef) == 3 and len(ref_ecef) == 3, "ECEF coordinates must have dimension 3."
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ecef_ = np.array(ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(ref_ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ned = np.empty(3, dtype=np.double)
  wgsecef2ned_d(&ecef_[0], &ref_ecef_[0], &ned[0])
  return ned

def wgsned2ecef_(ned, ref_ecef):
  """
  Wraps function :libswiftnav:`wgsned2ecef`.

  Parameters
  ----------
  ned : (float, float, float)
    The tuple of coordinates, `(n, e, d)`
  ref_ecef : (float, float, float)
    The tuple of coordinates of the reference position, `(x, y, z)`

  Returns
  -------
  out : :class:`numpy.ndarray`, shape(3,)
    The array `[x, y, z]`.

  """

  assert len(ned) == 3 and len(ref_ecef) == 3, "ECEF coordinates must have dimension 3."
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ned_ = np.array(ned, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(ref_ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ecef = np.empty(3, dtype=np.double)
  wgsned2ecef(&ned_[0], &ref_ecef_[0], &ecef[0])
  return ecef

def wgsned2ecef_d_(ned, ref_ecef):
  """
  Wraps function :libswiftnav:`wgsned2ecef_d`.

  Parameters
  ----------
  ned : (float, float, float)
    The tuple of coordinates, `(n, e, d)`
  ref_ecef : (float, float, float)
    The tuple of coordinates of the reference position, `(x, y, z)`

  Returns
  -------
  out : :class:`numpy.ndarray`, shape(3,)
    The array `[x, y, z]`.

  """
  assert len(ned) == 3 and len(ref_ecef) == 3, "ECEF coordinates must have dimension 3."
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ned_ = np.array(ned, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(ref_ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ecef = np.empty(3, dtype=np.double)
  wgsned2ecef_d(&ned_[0], &ref_ecef_[0], &ecef[0])
  return ecef

def wgsecef2azel_(ecef, ref_ecef):
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
  assert len(ecef) == 3 and len(ref_ecef) == 3, "ECEF coordinates must have dimension 3."
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ecef_ = np.array(ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(ref_ecef, dtype=np.double)
  cdef double az, el
  wgsecef2azel(&ecef_[0], &ref_ecef_[0], &az, &el)
  return (az, el)

def ecef2ned_matrix_(ref_ecef, M):
  assert len(ref_ecef) == 3, "ECEF coordinates must have dimension 3."
  assert M.shape == (3, 3), "M must be a 3x3 matrix."
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = np.array(ref_ecef, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] M_ = np.array(M, dtype=np.double)
  ecef2ned_matrix_(ref_ecef_[0], M[0][0])
  return M
