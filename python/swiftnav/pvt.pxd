# Copyright (C) 2012 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from gpstime cimport gps_time_t
from track cimport navigation_measurement_t

cdef extern from "libswiftnav/pvt.h":
  enum:
    PVT_MAX_ITERATIONS
    PVT_CONVERGED_NO_RAIM
    PVT_CONVERGED_RAIM_REPAIR
    PVT_CONVERGED_RAIM_OK
    PVT_PDOP_TOO_HIGH
    PVT_BAD_ALTITUDE
    PVT_VELOCITY_LOCKOUT
    PVT_RAIM_REPAIR_FAILED
    PVT_RAIM_REPAIR_IMPOSSIBLE
    PVT_UNCONVERGED
    PVT_INSUFFICENT_MEAS

  ctypedef struct dops_t:
    double pdop
    double gdop
    double tdop
    double hdop
    double vdop

  ctypedef struct gnss_solution:
    double pos_llh[3]
    double pos_ecef[3]
    double vel_ned[3]
    double vel_ecef[3]
    double err_cov[7]
    double clock_offset
    double clock_bias
    gps_time_t time
    u8 valid
    u8 n_used

  u8 calc_PVT(u8 n_used,
              navigation_measurement_t nav_meas[],
              u8 disable_raim,
              gnss_solution *soln,
              dops_t *dops)

cdef class GNSSSolution:
  cdef gnss_solution _thisptr

cdef class DOPS:
  cdef dops_t _thisptr
