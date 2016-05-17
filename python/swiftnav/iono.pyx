# Copyright (C) 2016 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from time cimport *
from time import GpsTime

def calc_iono_correction(GpsTime t, lat, lon, azimuth, elevation, alpha, beta):
  """
  Wraps function :libswiftnav:`calc_ionosphere`.

  Parameters
  ----------
  t : GpsTime
    GPS time for iono computation (only time of day used)
  lat, lon : float
    User latitude and longitude (rad)
  elevation, azimuth : float
    Satellite az-el (rad)
  alpha, beta : (float, float, float, float)
    ionospheric parameters

  Returns
  -------
  out : float
    Ionospheric delay in meters

  """
  assert len(alpha) == 4 and len(beta) == 4, "Ionospheric parameter vectors alpha and beta must have 4 elements each."
  
  cdef ionosphere_t i
  
  i.a0, i.a1, i.a2, i.a3 = alpha[0], alpha[1], alpha[2], alpha[3]
  i.b0, i.b1, i.b2, i.b3 = beta[0], beta[1], beta[2], beta[3]
  return calc_ionosphere(&t._thisptr, 
                         lat, lon,
                         azimuth, elevation,
                         &i);
