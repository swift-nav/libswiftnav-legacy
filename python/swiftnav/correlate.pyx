# Copyright (C) 2012,2016 Swift Navigation Inc.
# Contact: Adel Mamin <adelm@exafore.com>
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

# cython: profile=True

"""Correlation

Correlators used for tracking.
"""

cimport cython
cimport libc.math
cimport numpy as np
from common cimport *
from libc.math cimport M_PI
import numpy as np

cdef extern from "complexobject.h":
  struct Py_complex:
    double real
    double imag

  ctypedef class __builtin__.complex [object PyComplexObject]:
    cdef Py_complex cval

def track_correlate(np.ndarray[char, ndim=1, mode="c"] samples,
                    code_freq, code_phase, carr_freq, carr_phase,
                    np.ndarray[char, ndim=1, mode="c"] code,
                    sampling_freq, signal):
  cdef double init_code_phase = code_phase
  cdef double init_carr_phase = carr_phase
  cdef double I_E, Q_E, I_P, Q_P, I_L, Q_L
  cdef unsigned int blksize
  if signal == "l1ca":
    l1_ca_track_correlate(<s8*>&samples[0], len(samples), <s8*>&code[0],
                  &init_code_phase, code_freq / sampling_freq,
                  &init_carr_phase, carr_freq * 2.0 * M_PI / sampling_freq,
                  &I_E, &Q_E, &I_P, &Q_P, &I_L, &Q_L, &blksize)
  elif signal == "gloca":
    glo_ca_track_correlate(<s8*>&samples[0], len(samples), <s8*>&code[0],
                  &init_code_phase, code_freq / sampling_freq,
                  &init_carr_phase, carr_freq * 2.0 * M_PI / sampling_freq,
                  &I_E, &Q_E, &I_P, &Q_P, &I_L, &Q_L, &blksize)
  else:
    l2c_cm_track_correlate(<s8*>&samples[0], len(samples), <s8*>&code[0],
                  &init_code_phase, code_freq / sampling_freq,
                  &init_carr_phase, carr_freq * 2.0 * M_PI / sampling_freq,
                  &I_E, &Q_E, &I_P, &Q_P, &I_L, &Q_L, &blksize)
  # TODO: Pass pointers to Python complex number real and imag parts directly
  # to C function
  E = I_E + Q_E*1.j
  P = I_P + Q_P*1.j
  L = I_L + Q_L*1.j
  return (E, P, L, blksize, init_code_phase, init_carr_phase)
