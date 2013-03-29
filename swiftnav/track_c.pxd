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
    float y

  ctypedef struct simple_tl_state_t:
    float code_freq
    float carr_freq

  ctypedef struct comp_tl_state_t:
    float code_freq
    float carr_freq

  ctypedef struct cn0_est_state_t:
    pass

  ctypedef struct correlation_t:
    float I
    float Q

  void calc_navigation_measurement_(u8 n_channels, channel_measurement_t* meas[], navigation_measurement_t* nav_meas[], double nav_time, ephemeris_t* ephemerides[])

  void calc_loop_gains(float bw, float zeta, float k, float loop_freq,
                       float *pgain, float *igain)
  float costas_discriminator(float I, float Q)
  float dll_discriminator(correlation_t cs[3])

  void simple_lf_init(simple_lf_state_t *s, float y0,
                      float pgain, float igain)
  float simple_lf_update(simple_lf_state_t *s, float error)

  void simple_tl_init(simple_tl_state_t *s, float loop_freq,
                      float code_freq, float code_bw,
                      float code_zeta, float code_k,
                      float carr_freq, float carr_bw,
                      float carr_zeta, float carr_k)
  void simple_tl_update(simple_tl_state_t *s, correlation_t cs[3])

  void comp_tl_init(comp_tl_state_t *s, float loop_freq,
                      float code_freq, float code_bw,
                      float code_zeta, float code_k,
                      float carr_freq, float carr_bw,
                      float carr_zeta, float carr_k,
                      float tau, float cpc, u32 sched)
  void comp_tl_update(comp_tl_state_t *s, correlation_t cs[3])

  void cn0_est_init(cn0_est_state_t* s, float bw, float cn0_0,
                    float cutoff_freq, float loop_freq)
  float cn0_est(cn0_est_state_t* s, float I)



