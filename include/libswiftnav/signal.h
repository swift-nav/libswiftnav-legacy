/*
 * Copyright (c) 2016 Swift Navigation Inc.
 * Contact: Vlad Ungureanu <vvu@vdev.ro>
 *          Pasi Miettinen <pasi.miettinen@exafore.com>
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

#include <libswiftnav/common.h>

/** \addtogroup signal
 * \{ */

/* Number of satellites in each constellation. */
#define NUM_SATS_GPS 32
#define NUM_SATS_SBAS 19
#define NUM_SATS (NUM_SATS_GPS + NUM_SATS_SBAS)

/* Number of codes in each constellation. */
#define NUM_CODES_GPS 2
#define NUM_CODES_SBAS 1

/* Number of signals in each code. */
#define NUM_SIGNALS_GPS_L1CA NUM_SATS_GPS
#define NUM_SIGNALS_GPS_L2CM NUM_SATS_GPS
#define NUM_SIGNALS_SBAS_L1CA NUM_SATS_SBAS

/* Number of signals in each constellation. */
#define NUM_SIGNALS_GPS (NUM_SIGNALS_GPS_L1CA + NUM_SIGNALS_GPS_L2CM)
#define NUM_SIGNALS_SBAS (NUM_SIGNALS_SBAS_L1CA)
#define NUM_SIGNALS (NUM_SIGNALS_GPS + NUM_SIGNALS_SBAS)

#define GPS_FIRST_PRN 1
#define SBAS_FIRST_PRN 120

#define SID_STR_LEN_MAX 16

/** Constellation identifier. */
typedef enum constellation {
  CONSTELLATION_INVALID = -1,
  CONSTELLATION_GPS,
  CONSTELLATION_SBAS,
  CONSTELLATION_COUNT,
} constellation_t;

/** Code identifier. */
typedef enum code {
  CODE_INVALID = -1,
  CODE_GPS_L1CA,
  CODE_GPS_L2CM,
  CODE_SBAS_L1CA,
  CODE_COUNT,
} code_t;

/** GNSS signal identifier. */
typedef struct {
  u16 sat;
  code_t code;
} gnss_signal_t;

/** Signal comparison function. */
static inline int sid_compare(const gnss_signal_t a, const gnss_signal_t b)
{
  s32 x;
  if ((x = a.code - b.code))
    return x;
  return a.sat - b.sat;
}

/** Untyped signal comparison function. */
static inline int cmp_sid_sid(const void *a, const void *b)
{
  return sid_compare(*(gnss_signal_t *)a, *(gnss_signal_t*)b);
}

/** Signal equality function. */
static inline bool sid_is_equal(const gnss_signal_t a, const gnss_signal_t b)
{
  return sid_compare(a, b) == 0;
}

/* \} */

gnss_signal_t construct_sid(code_t code, u16 sat);
int sid_to_string(char *s, int n, gnss_signal_t sid);
bool sid_valid(gnss_signal_t sid);
bool code_valid(code_t code);
bool constellation_valid(constellation_t constellation);
gnss_signal_t sid_from_code_index(code_t code, u16 code_index);
u16 sid_to_code_index(gnss_signal_t sid);
enum constellation sid_to_constellation(gnss_signal_t sid);
enum constellation code_to_constellation(code_t code);
double code_to_carr_freq(code_t code);

#endif /* LIBSWIFTNAV_SIGNAL_H */
