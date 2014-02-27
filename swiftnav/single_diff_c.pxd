# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from almanac_c cimport *
from gpstime_c cimport *

cdef extern from "libswiftnav/single_diff.h":
  ctypedef struct sdiff_t:
    double pseudorange
    double carrier_phase
    double doppler
    double *sat_pos
    double *sat_vel
    double snr
    u8 prn

  void almanacs_to_single_diffs(u8 n, almanac_t *alms, gps_time_t timestamp, sdiff_t *sdiffs)
