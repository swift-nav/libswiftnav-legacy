# Copyright (C) 2016 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from time cimport gps_time_t

cdef extern from "libswiftnav/ionosphere.h":
  double calc_ionosphere(const gps_time_t *t_gps,
                       double lat_u, double lon_u,
                       double a, double e,
                       const ionosphere_t *i);

  ctypedef struct ionosphere_t:
    double a0, a1, a2, a3;
    double b0, b1, b2, b3;

