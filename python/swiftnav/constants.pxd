# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from constants cimport *

cdef extern from "libswiftnav/constants.h":
  # HACK: This is declared as an enum so that Cython knows its a constant
  # https://groups.google.com/d/msg/cython-users/-fLG08E5lYM/xC93UHSvLF0J
  enum:
    MAX_CHANNELS
    MAX_SATS

  float R2D
  float D2R

  float GPS_PI
  float GPS_L1_HZ
  float GPS_OMEGAE_DOT
  float GPS_GM
  float GPS_C
  float GPS_F
  float GPS_C_NO_VAC
  float GPS_L1_LAMBDA
  float GPS_L1_LAMBDA_NO_VAC
  float GPS_NOMINAL_RANGE
  float GPS_CA_CHIPPING_RATE

  float DEFAULT_PHASE_VAR_TEST
  float DEFAULT_CODE_VAR_TEST
  float DEFAULT_PHASE_VAR_KF
  float DEFAULT_CODE_VAR_KF
  float DEFAULT_AMB_DRIFT_VAR
  float DEFAULT_AMB_INIT_VAR
  float DEFAULT_NEW_INT_VAR
