/*
 * Copyright (C) 2012 Swift Navigation Inc.
 * Contact: Fergus Noble <fergus@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef LIBSWIFTNAV_TRACK_H
#define LIBSWIFTNAV_TRACK_H

#include "common.h"
#include "ephemeris.h"

/** \addtogroup track
 * \{ */

/** \addtogroup track_loop
 * \{ */

/** State structure for the simple loop filter.
 * Should be initialised with simple_lf_init().
 */
typedef struct {
  float pgain;      /**< Proportional gain. */
  float igain;      /**< Integral gain. */
  float prev_error; /**< Previous error. */
  float y;          /**< Output variable. */
} simple_lf_state_t;

/** State structure for a simple tracking loop.
 * Should be initialised with simple_tl_init().
 */
typedef struct {
  float code_freq;             /**< Code phase rate (i.e. frequency). */
  float carr_freq;             /**< Carrier frequency. */
  simple_lf_state_t code_filt; /**< Code loop filter state. */
  simple_lf_state_t carr_filt; /**< Carrier loop filter state. */
} simple_tl_state_t;

/** State structure for a code/carrier phase complimentary filter tracking
 * loop.
 * Should be initialised with comp_tl_init().
 */
typedef struct {
  float code_freq;             /**< Code phase rate (i.e. frequency). */
  float carr_freq;             /**< Carrier frequency. */
  simple_lf_state_t code_filt; /**< Code loop filter state. */
  simple_lf_state_t carr_filt; /**< Carrier loop filter state. */
  u32 sched;                   /**< Gain scheduling count. */
  u32 n;                       /**< Iteration counter. */
  float A;                     /**< Complementary filter crossover gain. */
  float carr_to_code;          /**< Scale factor from carrier to code. */
} comp_tl_state_t;


/** \} */

/** Structure representing a complex valued correlation. */
typedef struct {
  float I; /**< In-phase correlation. */
  float Q; /**< Quadrature correlation. */
} correlation_t;

/** State structure for the \f$ C / N_0 \f$ estimator.
 * Should be initialised with cn0_est_init().
 */
typedef struct {
  float log_bw;     /**< Noise bandwidth in dBHz. */
  float A;          /**< IIR filter coeff. */
  float I_prev_abs; /**< Abs. value of the previous in-phase correlation. */
  float nsr;        /**< Noise-to-signal ratio (1 / SNR). */
} cn0_est_state_t;

/** \} */

typedef struct {
  u8 prn;
  double code_phase_chips;
  double code_phase_rate;
  double carrier_phase;
  double carrier_freq;
  u32 time_of_week_ms;
  double receiver_time;
  double snr;
} channel_measurement_t;

typedef struct __attribute__((packed)) {
  double raw_pseudorange;
  double raw_pseudorange_rate;
  double pseudorange;
  double pseudorange_rate;
  double carrier_phase;
  double doppler;
  double sat_pos[3];
  double sat_vel[3];
  double snr;
  gps_time_t tot;
  u8 prn;
  u8 padding;
} navigation_measurement_t;

void calc_loop_gains(float bw, float zeta, float k, float loop_freq,
                     float *pgain, float *igain);
float costas_discriminator(float I, float Q);
float dll_discriminator(correlation_t cs[3]);

void simple_lf_init(simple_lf_state_t *s, float y0,
                    float pgain, float igain);
float simple_lf_update(simple_lf_state_t *s, float error);

void simple_tl_init(simple_tl_state_t *s, float loop_freq,
                    float code_freq, float code_bw,
                    float code_zeta, float code_k,
                    float carr_freq, float carr_bw,
                    float carr_zeta, float carr_k);
void simple_tl_update(simple_tl_state_t *s, correlation_t cs[3]);

void comp_tl_init(comp_tl_state_t *s, float loop_freq,
                    float code_freq, float code_bw,
                    float code_zeta, float code_k,
                    float carr_freq, float carr_bw,
                    float carr_zeta, float carr_k,
                    float tau, float cpc, u32 sched);
void comp_tl_update(comp_tl_state_t *s, correlation_t cs[3]);

void cn0_est_init(cn0_est_state_t *s, float bw, float cn0_0,
                  float cutoff_freq, float loop_freq);
float cn0_est(cn0_est_state_t *s, float I);

void calc_navigation_measurement(u8 n_channels, channel_measurement_t meas[],
                                 navigation_measurement_t nav_meas[],
                                 double nav_time, ephemeris_t ephemerides[]);
void calc_navigation_measurement_(u8 n_channels, channel_measurement_t* meas[],
                                  navigation_measurement_t* nav_meas[],
                                  double nav_time, ephemeris_t* ephemerides[]);

#endif /* LIBSWIFTNAV_TRACK_H */

