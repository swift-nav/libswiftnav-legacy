# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from observation cimport *

cdef extern from "libswiftnav/filter_utils.h":
  double simple_amb_measurement(double carrier, double code)
  s8 assign_de_mtx(u8 num_sats, const sdiff_t *sats_with_ref_first,
                   const double ref_ecef[3], double *DE)
