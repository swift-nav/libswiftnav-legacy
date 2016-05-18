# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from time cimport *
from time import GpsTime

def calc_tropo_correction(GpsTime t, lat, h, el):
  """
  Wraps function :libswiftnav:`calc_troposphere`.

  Parameters
  ----------
  t : GpsTime 
    Time at which to calculate tropospheric delay (only doy matters)
  lat : float
    Latitude of the receiver [rad]
  h : float
    Orthometric height of the receiver [m]
  el : float
    Elevation of the satellite [rad]

  Returns
  -------
  out : float
    Tropospheric delay in meters

  """
  
  return calc_troposphere(&t._thisptr, lat, h, el)
