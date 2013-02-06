# Copyright (C) 2012 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from ephemeris_c cimport ephemeris_t

cdef extern from "libswiftnav/track.h":
  ctypedef struct channel_measurement_t:
    u8 prn
    double code_phase_chips
    double code_phase_rate
    double carrier_phase
    double carrier_freq
    u32 time_of_week_ms
    double receiver_time
    double snr

  ctypedef struct navigation_measurement_t:
    double pseudorange
    double pseudorange_rate
    double TOT
    double sat_pos[3]
    double sat_vel[3]

  ctypedef struct simple_lf_state_t:
    double y

  ctypedef struct correlation_t:
    double I
    double Q

  void calc_navigation_measurement_(u8 n_channels, channel_measurement_t* meas[], navigation_measurement_t* nav_meas[], double nav_time, ephemeris_t* ephemerides[])

  void calc_loop_gains(double bw, double zeta, double k, double loop_freq,
                       double *pgain, double *igain)
  double costas_discriminator(double I, double Q)
  double dll_discriminator(correlation_t cs[3])

  void simple_lf_init(simple_lf_state_t *s, double y0,
                      double pgain, double igain)
  double simple_lf_update(simple_lf_state_t *s, double error)


