# Copyright (C) 2015,2016 Swift Navigation Inc.
# Contact: Adel Mamin <adelm@exafore.com>
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *

cdef extern from "libswiftnav/correlate.h":
  void l1_ca_track_correlate(s8* samples, size_t samples_len, s8* code,
                       u32 chips_to_correlate,
                       double* init_code_phase, double code_step,
                       double* init_carr_phase, double carr_step,
                       double* I_E, double* Q_E,
                       double* I_P, double* Q_P,
                       double* I_L, double* Q_L,
                       u32* num_samples)

  void l2c_cm_track_correlate(s8* samples, size_t samples_len, s8* code,
                       u32 chips_to_correlate,
                       double* init_code_phase, double code_step,
                       double* init_carr_phase, double carr_step,
                       double* I_E, double* Q_E,
                       double* I_P, double* Q_P,
                       double* I_L, double* Q_L,
                       u32* num_samples)

  void glo_ca_track_correlate(s8* samples, size_t samples_len, s8* code,
                       double* init_code_phase, double code_step,
                       double* init_carr_phase, double carr_step,
                       double* I_E, double* Q_E,
                       double* I_P, double* Q_P,
                       double* I_L, double* Q_L,
                       u32* num_samples);

