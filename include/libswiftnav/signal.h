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

#define GPS_L1_SATS 32
#define WAAS_SATS 3
#define EGNOS_SATS 3
#define GAGAN_SATS 2
#define MSAS_SATS 2
#define SDCM_SATS 2
#define SBAS_SATS (WAAS_SATS + EGNOS_SATS + GAGAN_SATS + MSAS_SATS + SDCM_SATS)
#define ALL_SATS (GPS_L1_SATS + SBAS_SATS)

#define GPS_CONSTELLATION 0
#define SBAS_CONSTELLATION 1

#define L1_BAND 0

typedef struct __attribute__((packed)) {
  u16 sat;
  u8 band;
  u8 constellation;
} signal_t;

static inline int cmp_signal_signal(const void *a, const void *b)
{
  const signal_t *a_ = a, *b_ = b;
  s32 x;
  if ((x = a_->constellation - b_->constellation))
    return x;
  if ((x = a_->band - b_->band))
    return x;
  return a_->sat - b_->sat;
}

static inline bool signal_is_equal(const signal_t a, const signal_t b)
{
  return cmp_signal_signal(&a, &b) == 0;
}

#endif /* LIBSWIFTNAV_SIGNAL_H */

