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
  double pgain;      /**< Proportional gain. */
  double igain;      /**< Integral gain. */
  double prev_error; /**< Previous error. */
  double y;          /**< Output variable. */
} simple_lf_state_t;

/** State structure for a simple tracking loop.
 * Should be initialised with simple_tl_init().
 */
typedef struct {
  double code_freq;            /**< Code phase rate (i.e. frequency). */
  double carr_freq;            /**< Carrier frequency. */
  simple_lf_state_t code_filt; /**< Code loop filter state. */
  simple_lf_state_t carr_filt; /**< Carrier loop filter state. */
} simple_tl_state_t;

/** \} */

/** Structure representing a complex valued correlation. */
typedef struct {
  double I; /**< In-phase correlation. */
  double Q; /**< Quadrature correlation. */
} correlation_t;

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

typedef struct {
  double pseudorange;
  double pseudorange_rate;
  double TOT;
  double sat_pos[3];
  double sat_vel[3];
} navigation_measurement_t;

void calc_loop_gains(double bw, double zeta, double k, double loop_freq,
                     double *pgain, double *igain);
double costas_discriminator(double I, double Q);
double dll_discriminator(correlation_t cs[3]);

void simple_lf_init(simple_lf_state_t *s, double y0,
                    double pgain, double igain);
double simple_lf_update(simple_lf_state_t *s, double error);

void simple_tl_init(simple_tl_state_t *s, double loop_freq,
                    double code_freq, double code_bw,
                    double code_zeta, double code_k,
                    double carr_freq, double carr_bw,
                    double carr_zeta, double carr_k);
void simple_tl_update(simple_tl_state_t *s, correlation_t cs[3]);

void calc_navigation_measurement(u8 n_channels, channel_measurement_t meas[],
                                 navigation_measurement_t nav_meas[],
                                 double nav_time, ephemeris_t ephemerides[]);
void calc_navigation_measurement_(u8 n_channels, channel_measurement_t* meas[],
                                  navigation_measurement_t* nav_meas[],
                                  double nav_time, ephemeris_t* ephemerides[]);

#endif /* LIBSWIFTNAV_TRACK_H */

