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

#define NUM_SIGNALS_GPS (NUM_SATS_GPS * NUM_CODES_GPS)
#define NUM_SIGNALS_SBAS (NUM_SATS_SBAS * NUM_CODES_SBAS)
#define NUM_SIGNALS (NUM_SIGNALS_GPS + NUM_SIGNALS_SBAS)

#define GPS_FIRST_PRN 1
#define SBAS_FIRST_PRN 120

#define INVALID_SID_INDEX -1

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

enum indexing_type {
  INDEX_GLOBAL, /* [0, NUM_SIGNALS) index of the signal inside the full virtual
                 * list of all satellites  and codes */
  INDEX_CONSTELLATION , /* [0, NUM_SIGNALS_CONSTELLATION) index of the signal
                         * inside the constellation*/
  INDEX_SAT_IN_CONS, /* [0, NUM_SATS_CONSTELLATION) index of the satellite
                      * transmitting this sid inside the constellation it
                      * belongs to*/
};

#define GPS_FIRST_PRN 1
#define SBAS_FIRST_PRN 120

#define SID_STR_LEN_MAX 16

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

/** Returns true if sid corresponds to a known constellation, band, and
 *  satellite identifier.
 */
bool sid_valid(gnss_signal_t sid);

/** Converts the global index i in [0, NUM_SIGNALS) to the corresponding sid
 */
gnss_signal_t sid_from_index(u16 i);

/** Converts sid to a index according to given indexing type.
 */
u16 sid_to_index(gnss_signal_t sid, enum indexing_type it);

/** Returns the constellation in which the sid belongs to.
 */
enum constellation sid_to_constellation(gnss_signal_t sid);

#endif /* LIBSWIFTNAV_SIGNAL_H */
