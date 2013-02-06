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

typedef struct {
  double I, Q;
} correlation_t;

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


void calc_navigation_measurement(u8 n_channels, channel_measurement_t meas[],
                                 navigation_measurement_t nav_meas[],
                                 double nav_time, ephemeris_t ephemerides[]);
void calc_navigation_measurement_(u8 n_channels, channel_measurement_t* meas[],
                                  navigation_measurement_t* nav_meas[],
                                  double nav_time, ephemeris_t* ephemerides[]);

#endif /* LIBSWIFTNAV_TRACK_H */

