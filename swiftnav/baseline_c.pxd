# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from observation_c cimport *

cdef extern from "libswiftnav/amb_kf.h":
  ctypedef struct ambiguities_t:
    double *ambs
    u8 *prns
    u8 n
  s8 baseline(u8 num_sdiffs, const sdiff_t *sdiffs, const double ref_ecef[3],
              const ambiguities_t *ambs, u8 *num_used, double b[3])

