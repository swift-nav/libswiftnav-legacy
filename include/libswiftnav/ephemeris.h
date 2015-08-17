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
} ephemeris_t;

s8 calc_sat_state(const ephemeris_t *ephemeris, gps_time_t t,
                  double pos[3], double vel[3],
                  double *clock_err, double *clock_rate_err);

u8 ephemeris_good(ephemeris_t *eph, gps_time_t t);

void decode_ephemeris(u32 frame_words[3][8], ephemeris_t *e);
bool ephemeris_equal(ephemeris_t *a, ephemeris_t *b);

#endif /* LIBSWIFTNAV_EPHEMERIS_H */

