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

#include <libswiftnav/signal.h>
#include <libswiftnav/time.h>
#include <libswiftnav/common.h>

typedef struct {
  double tgd;
  double crs, crc, cuc, cus, cic, cis;
  double dn, m0, ecc, sqrta, omega0, omegadot, w, inc, inc_dot;
  double af0, af1, af2;
  gps_time_t toc;
  u16 iodc;
  u8 iode;
} ephemeris_kepler_t;

typedef struct {
  double pos[3];
  double rate[3];
  double acc[3];
  u8 iod;
  u16 toa;
  double a_gf0;
  double a_gf1;
} ephemeris_xyz_t;

typedef struct {
  gnss_signal_t sid;
  gps_time_t toe;
  float ura;
  u8 fit_interval;
  u8 valid;
  u8 healthy;
  union {
    ephemeris_kepler_t kepler;
    ephemeris_xyz_t xyz;
  };
} ephemeris_t;

s8 calc_sat_state(const ephemeris_t *e, const gps_time_t *t,
                  double pos[3], double vel[3],
                  double *clock_err, double *clock_rate_err);

u8 ephemeris_valid(const ephemeris_t *eph);
u8 ephemeris_valid_at_time(const ephemeris_t *eph, const gps_time_t *t);
u8 satellite_healthy(const ephemeris_t *eph);

float decode_ura_index(const u8 index);
u8 decode_fit_interval(u8 fit_interval_flag, u16 iodc);
void decode_ephemeris(u32 frame_words[3][8], ephemeris_t *e);
bool ephemeris_equal(const ephemeris_t *a, const ephemeris_t *b);

#endif /* LIBSWIFTNAV_EPHEMERIS_H */
