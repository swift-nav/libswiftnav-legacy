# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *

cdef extern from "libswiftnav/correlate.h":
  void track_correlate(s8* samples, s8* code,
                       double* init_code_phase, double code_step,
                       double* init_carr_phase, double carr_step,
                       double* I_E, double* Q_E,
                       double* I_P, double* Q_P,
                       double* I_L, double* Q_L,
                       u32* num_samples)
