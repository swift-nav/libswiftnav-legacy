# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from ephemeris cimport ephemeris_t
from time cimport gps_time_t
from libcpp cimport bool
from signal cimport gnss_signal_t

cdef extern from "libswiftnav/track.h":

  ctypedef struct correlation_t:
    float I
    float Q

  # Discriminators
  void calc_loop_gains(float bw, float zeta, float k, float loop_freq, float *b0, float *b1)
  float costas_discriminator(float I, float Q)
  float frequency_discriminator(float I, float Q, float prev_I, float prev_Q)
  float dll_discriminator(correlation_t cs[3])

  # Tracking loop: Aided loop filter
  ctypedef struct aided_lf_state_t:
    float b0
    float b1
    float aiding_igain
    float prev_error
    float y

  void aided_lf_init(aided_lf_state_t *s, float y0, float pgain, float igain, float aiding_igain)
  float aided_lf_update(aided_lf_state_t *s, float p_i_error, float aiding_error)

  # Tracking loop: Simple loop (tracking) filter
  ctypedef struct simple_lf_state_t:
    float b0
    float b1
    float prev_error
    float y

  void simple_lf_init(simple_lf_state_t *s, float y0, float pgain, float igain)
  float simple_lf_update(simple_lf_state_t *s, float error)

  ctypedef struct simple_tl_state_t:
    float code_freq
    float carr_freq
    simple_lf_state_t code_filt
    simple_lf_state_t carr_filt

  void simple_tl_init(simple_tl_state_t *s,
                      float loop_freq,
                      float code_freq,
                      float code_bw,
                      float code_zeta,
                      float code_k,
                      float carr_freq,
                      float carr_bw,
                      float carr_zeta,
                      float carr_k)
  void simple_tl_update(simple_tl_state_t *s, correlation_t cs[3])

  # Tracking Loop: Aided tracking loop
  ctypedef struct aided_tl_state_t:
    float carr_freq
    aided_lf_state_t carr_filt
    float code_freq
    simple_lf_state_t code_filt
    float prev_I, prev_Q
    float carr_to_code

  void aided_tl_init(aided_tl_state_t *s,
                     float loop_freq,
                     float code_freq,
                     float code_bw,
                     float code_zeta,
                     float code_k,
                     float carr_to_code,
                     float carr_freq,
                     float carr_bw,
                     float carr_zeta,
                     float carr_k,
                     float carr_freq_b1)
  void aided_tl_retune(aided_tl_state_t *s, float loop_freq,
                       float code_bw, float code_zeta, float code_k,
                       float carr_to_code,
                       float carr_bw, float carr_zeta, float carr_k,
                       float carr_freq_b1)
  void aided_tl_update(aided_tl_state_t *s, correlation_t cs[3])

  # Tracking loop: Comp
  ctypedef struct comp_tl_state_t:
    float code_freq
    float carr_freq
    simple_lf_state_t code_filt
    simple_lf_state_t carr_filt
    u32 sched
    u32 n
    float A
    float carr_to_code

  void comp_tl_init(comp_tl_state_t *s, float loop_freq,
                    float code_freq, float code_bw,
                    float code_zeta, float code_k,
                    float carr_freq, float carr_bw,
                    float carr_zeta, float carr_k,
                    float tau, float cpc, u32 sched)
  void comp_tl_update(comp_tl_state_t *s, correlation_t cs[3])

  # Tracking loop: Alias detect
  ctypedef struct alias_detect_t:
    u8 acc_len
    float dt
    float dot
    float cross
    u8 fl_count
    float first_I
    float first_Q

  void alias_detect_init(alias_detect_t *a, u32 acc_len, float time_diff)
  void alias_detect_first(alias_detect_t *a, float I, float Q)
  float alias_detect_second(alias_detect_t *a, float I, float Q)

  # Tracking loop: Lock detect
  struct loop_detect_lpf:
    float k1
    float y

  ctypedef struct lock_detect_t:
    loop_detect_lpf lpfi
    loop_detect_lpf lpfq
    float k2
    u16 lo, lp
    u16 pcount1, pcount2
    bool outo, outp

  void lock_detect_init(lock_detect_t *l, float k1, float k2, u16 lp, u16 lo)
  void lock_detect_reinit(lock_detect_t *l, float k1, float k2, u16 lp, u16 lo)
  void lock_detect_update(lock_detect_t *l, float I, float Q, float DT)

  # Tracking loop: CN0 est
  ctypedef struct cn0_est_state_t:
    float log_bw
    float b
    float a
    float I_prev_abs
    float Q_prev_abs
    float nsr
    float xn

  void cn0_est_init(cn0_est_state_t *s, float bw, float cn0_0, float cutoff_freq, float loop_freq)
  float cn0_est(cn0_est_state_t *s, float I, float Q)

  # Tracking loop: Navigation measurement
  ctypedef struct channel_measurement_t:
    gnss_signal_t sid
    double code_phase_chips
    double code_phase_rate
    double carrier_phase
    double carrier_freq
    u32 time_of_week_ms
    double receiver_time
    double snr
    u16 lock_counter

  ctypedef struct navigation_measurement_t:
    double raw_pseudorange
    double pseudorange
    double raw_carrier_phase
    double carrier_phase
    double raw_doppler
    double doppler
    double sat_pos[3]
    double sat_vel[3]
    double snr
    double lock_time
    gps_time_t tot
    gnss_signal_t sid
    u16 lock_counter

  void calc_navigation_measurement(u8 n_channels, const channel_measurement_t *meas[],
                                   navigation_measurement_t *nav_meas[],
                                   double nav_time_tc, gps_time_t *nav_time,
                                   const e* e[])
  int nav_meas_cmp(const void *a, const void *b)
  u8 tdcp_doppler(u8 n_new, navigation_measurement_t *m_new,
                  u8 n_old, navigation_measurement_t *m_old,
                  navigation_measurement_t *m_corrected)

cdef class Correlation:
  cdef correlation_t _thisptr

cdef class SimpleLoopFilter:
  cdef simple_lf_state_t _thisptr

cdef class SimpleTrackingLoop:
  cdef simple_tl_state_t _thisptr

cdef class AidedLoopFilter:
   cdef aided_lf_state_t _thisptr

cdef class AidedTrackingLoop:
  cdef aided_tl_state_t _thisptr

cdef class CompTrackingLoop:
  cdef comp_tl_state_t _thisptr

cdef class AliasDetector:
  cdef alias_detect_t _thisptr

cdef class LockDetectLoopFilter:
  cdef loop_detect_lpf _thisptr

cdef class LockDetector:
  cdef lock_detect_t _thisptr

cdef class CN0Estimator:
  cdef cn0_est_state_t _thisptr

cdef class ChannelMeasurement:
  cdef channel_measurement_t _thisptr

cdef class NavigationMeasurement:
 cdef navigation_measurement_t _thisptr

cdef mk_nav_meas_array(py_nav_meas, u8 n_c_nav_meas, navigation_measurement_t *c_nav_meas)
