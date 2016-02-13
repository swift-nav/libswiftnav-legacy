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

#define NUM_SATS_GPS 32
#define NUM_SATS_SBAS 19
#define NUM_SATS (NUM_SATS_GPS + NUM_SATS_SBAS)

#define NUM_CODES_GPS 2
#define NUM_CODES_SBAS 1

#define NUM_SIGNALS_GPS_L1CA NUM_SATS_GPS
#define NUM_SIGNALS_GPS_L2CM NUM_SATS_GPS
#define NUM_SIGNALS_SBAS_L1CA NUM_SATS_SBAS

#define NUM_SIGNALS_GPS (NUM_SIGNALS_GPS_L1CA + NUM_SIGNALS_GPS_L2CM)
#define NUM_SIGNALS_SBAS (NUM_SIGNALS_SBAS_L1CA)
#define NUM_SIGNALS (NUM_SIGNALS_GPS + NUM_SIGNALS_SBAS)

#define GPS_FIRST_PRN 1
#define SBAS_FIRST_PRN 120

#define SID_STR_LEN_MAX 16

enum constellation {
  CONSTELLATION_INVALID = -1,
  CONSTELLATION_GPS,
  CONSTELLATION_SBAS,
  CONSTELLATION_COUNT,
};

enum code {
  CODE_INVALID = -1,
  CODE_GPS_L1CA,
  CODE_GPS_L2CM,
  CODE_SBAS_L1CA,
  CODE_COUNT,
};

typedef struct {
  u16 sat;
  enum code code;
} gnss_signal_t;

static inline int sid_compare(const gnss_signal_t a, const gnss_signal_t b)
{
  s32 x;
  if ((x = a.code - b.code))
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

/** Constructs and returns a gnss_signal_t structure using given parameters.
 */
gnss_signal_t construct_sid(enum code code, u16 sat);

/** Print a string representation of sid to the buffer s of capacity n.
 *  Returns the number of characters printed, excluding the terminating null.
 */
int sid_to_string(char *s, int n, gnss_signal_t sid);

/** Returns true if sid corresponds to a known constellation, code, and
 *  satellite identifier.
 */
bool sid_valid(gnss_signal_t sid);

/** Converts a code-specific index in [0, NUM_SIGNALS_<code>) to
 *  the corresponding sid.
 */
gnss_signal_t sid_from_code_index(enum code code, u16 code_index);

/** Returns the code index in [0, NUM_SIGNALS_<code>) corresponding to the sid.
 */
u16 sid_to_code_index(gnss_signal_t sid);

/** Returns the constellation which the sid belongs to.
 */
enum constellation sid_to_constellation(gnss_signal_t sid);

/** Returns the constellation which the code belongs to.
 */
enum constellation code_to_constellation(enum code code);

#endif /* LIBSWIFTNAV_SIGNAL_H */
