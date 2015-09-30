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
  u8 constellation;
  u8 band;
  u16 prn;
} signal_t;

signal_t sbas_index_to_sid(u8 index);
u8 sbas_sid_to_index(signal_t sid);

#endif /* LIBSWIFTNAV_SIGNAL_H */

