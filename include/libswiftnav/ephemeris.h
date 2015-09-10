/*
 * Copyright (C) 2010 Swift Navigation Inc.
 * Contact: Henry Hallam <henry@swift-nav.com>
 *          Fergus Noble <fergus@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef LIBSWIFTNAV_EPHEMERIS_H
#define LIBSWIFTNAV_EPHEMERIS_H

#include "signal.h"
#include "gpstime.h"
#include "common.h"

typedef struct {
  double tgd;
  double crs, crc, cuc, cus, cic, cis;
  double dn, m0, ecc, sqrta, omega0, omegadot, w, inc, inc_dot;
  double af0, af1, af2;
  gps_time_t toe, toc;
  u8 valid;
  u8 healthy;
  signal_t sid;
  u8 iode;
} ephemeris_kepler_t;

typedef struct {
  gps_time_t toe;
  double pos[3];
  double rate[3];
  double acc[3];
  u8 valid;
  u8 healthy;
  signal_t sid;
  u8 iod;
  u16 toa;
  u8 ura;
  double a_gf0;
  double a_gf1;
} ephemeris_xyz_t;

typedef struct {
  ephemeris_kepler_t *ephemeris_kep;
  ephemeris_xyz_t *ephemeris_xyz;
} ephemeris_t;

s8 calc_sat_state(const ephemeris_kepler_t *ephemeris, gps_time_t t,
                  double pos[3], double vel[3],
                  double *clock_err, double *clock_rate_err);

u8 ephemeris_good(ephemeris_t *eph, signal_t sid, gps_time_t t);

void decode_ephemeris(u32 frame_words[3][8], ephemeris_kepler_t *e);
bool ephemeris_xyz_equal(ephemeris_xyz_t *a, ephemeris_xyz_t *b);
bool ephemeris_kepler_equal(ephemeris_kepler_t *a, ephemeris_kepler_t *b);

#endif /* LIBSWIFTNAV_EPHEMERIS_H */

