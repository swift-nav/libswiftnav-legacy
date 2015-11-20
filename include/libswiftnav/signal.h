/*
 * Copyright (c) 2015 Swift Navigation Inc.
 * Contact: Vlad Ungureanu <vvu@vdev.ro>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef LIBSWIFTNAV_SIGNAL_H
#define LIBSWIFTNAV_SIGNAL_H

#include "common.h"

#define NUM_SATS_GPS 32
#define NUM_SATS_SBAS 22
#define NUM_SATS (NUM_SATS_GPS + NUM_SATS_SBAS)

enum constellation {
  CONSTELLATION_GPS,
  CONSTELLATION_SBAS,
  CONSTELLATION_COUNT,
};

enum band {
  BAND_L1,
  BAND_COUNT,
};

#define GPS_FIRST_PRN 1
#define SBAS_FIRST_PRN 120

typedef struct __attribute__((packed)) {
  u16 sat;
  u8 band;
  u8 constellation;
} gnss_signal_t;

static inline int sid_compare(const gnss_signal_t a, const gnss_signal_t b)
{
  s32 x;
  if ((x = a.constellation - b.constellation))
    return x;
  if ((x = a.band - b.band))
    return x;
  return a.sat - b.sat;
}

static inline int cmp_sid_sid(const void *a, const void *b)
{
  return sid_compare(*(gnss_signal_t *)a, *(gnss_signal_t*)b);
}

static inline bool sid_is_equal(const gnss_signal_t a, const gnss_signal_t b)
{
  return sid_compare(a, b) == 0;
}

#endif /* LIBSWIFTNAV_SIGNAL_H */

