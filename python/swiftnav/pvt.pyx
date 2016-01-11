# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from libc.stdlib cimport malloc, free
from libc.string cimport memset
from track cimport NavigationMeasurement
from track cimport navigation_measurement_t
import warnings

cdef class GNSSSolution:

  def __init__(self, **kwargs):
    memset(&self._thisptr, 0, sizeof(gnss_solution))
    if kwargs:
      self._thisptr = kwargs

  def __getattr__(self, k):
    return self._thisptr.get(k)


cdef class DOPS:

  def __init__(self, **kwargs):
    memset(&self._thisptr, 0, sizeof(dops_t))
    if kwargs:
      self._thisptr = kwargs

  def __getattr__(self, k):
    return self._thisptr.get(k)

_calc_pvt_codes = {2: "Solution converged but RAIM unavailable or disabled",
                   1: "Solution converged, failed RAIM but was successfully repaired",
                   -1: "PDOP is too high to yield a good solution.",
                   -2: "Altitude is unreasonable.",
                   -3: "Velocity is greater than or equal to 1000 kts.",
                   -4: "RAIM check failed and repair was unsuccessful",
                   -5: "RAIM check failed and repair was impossible (not enough measurements)",
                   -6: "pvt_iter didn't converge",
                   -7: "< 4 measurements"}

def calc_PVT_(nav_meas, disable_raim=False):
  """Wraps the function :libswiftnav:`calc_PVT`.

  Parameters
  ----------
  nav_meas : [NavigationMeasurement]
    Navigation measurements
  disable_raim : bool
    Disable RAIM

  Returns
  -------
  (GNSSSolution, DOPS)

  """
  n_used = len(nav_meas)
  cdef navigation_measurement_t* nav_meas_ = <navigation_measurement_t *>malloc(n_used*sizeof(navigation_measurement_t))
  for n in range(n_used):
    nav_meas_[n] = (<NavigationMeasurement?>nav_meas[n])._thisptr
  dops_ = DOPS()
  soln = GNSSSolution()
  cdef s8 ret = calc_PVT(n_used, nav_meas_, disable_raim, &soln._thisptr, &dops_._thisptr)
  if not ret == 0:
    warnings.warn(_calc_pvt_codes[ret])
  free(nav_meas_)
  return (ret, soln, dops_)
