/*
 * Copyright (c) 2015 Swift Navigation Inc.
 * Contact: Jacob McNamee <jacob@swiftnav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#include "signal.h"

#include <stdio.h>
#include <assert.h>

typedef struct {
  u8 constellation;
  u8 band;
  u16 sat_count;
  u16 (*sat_get)(u32 local_index);
  u32 (*local_index_get)(u16 sat);
} sat_group_t;

static const char * constellation_strs[CONSTELLATION_COUNT] = {
  [CONSTELLATION_GPS] = "GPS",
  [CONSTELLATION_SBAS] = "SBAS",
  [CONSTELLATION_QZSS] = "QZSS"
};

static const char * band_strs[BAND_COUNT] = {
  [BAND_L1] = "L1"
};

static const char * unknown_str = "?";

static u16 sat_get_gps(u32 local_index);
static u16 sat_get_sbas(u32 local_index);
static u16 sat_get_qzss(u32 local_index);
static u32 local_index_get_gps(u16 sat);
static u32 local_index_get_sbas(u16 sat);
static u32 local_index_get_qzss(u16 sat);

static const sat_group_t sat_groups[CONSTELLATION_COUNT] = {
  [CONSTELLATION_GPS] = {
    CONSTELLATION_GPS, BAND_L1, NUM_SATS_GPS,
     sat_get_gps, local_index_get_gps
  },
  [CONSTELLATION_SBAS] = {
    CONSTELLATION_SBAS, BAND_L1, NUM_SATS_SBAS,
    sat_get_sbas, local_index_get_sbas
  },
  [CONSTELLATION_QZSS] = {
    CONSTELLATION_QZSS, BAND_L1, NUM_SATS_QZSS,
    sat_get_qzss, local_index_get_qzss
  }
};

static u16 sat_get_gps(u32 local_index)
{
  return local_index + GPS_FIRST_PRN;
}

static u16 sat_get_sbas(u32 local_index)
{
  return local_index + SBAS_FIRST_PRN;
}

static u16 sat_get_qzss(u32 local_index)
{
  return local_index + QZSS_FIRST_PRN;
}

static u32 local_index_get_gps(u16 sat)
{
  return sat - GPS_FIRST_PRN;
}

static u32 local_index_get_sbas(u16 sat)
{
  return sat - SBAS_FIRST_PRN;
}

static u32 local_index_get_qzss(u16 sat)
{
  return sat - QZSS_FIRST_PRN;
}

int sid_to_string(char *s, int n, gnss_signal_t sid)
{
  const char *constellation_str = (sid.constellation < CONSTELLATION_COUNT) ?
      constellation_strs[sid.constellation] : unknown_str;
  const char *band_str = (sid.band < BAND_COUNT) ?
      band_strs[sid.band] : unknown_str;

  int nchars = snprintf(s, n, "%s %s %u", constellation_str, band_str, sid.sat);
  s[n-1] = 0;
  return nchars;
}

bool sid_valid(gnss_signal_t sid)
{
  if (sid.constellation >= CONSTELLATION_COUNT)
    return false;
  if (sid.band >= BAND_COUNT)
    return false;
  if (sat_groups[sid.constellation].local_index_get(sid.sat)
        >= sat_groups[sid.constellation].sat_count)
    return false;
  return true;
}

gnss_signal_t sid_from_index(u32 i)
{
  assert(i < NUM_SATS);
  gnss_signal_t sid = {0, 0, 0};
  u32 offset = i;
  u32 group_index = 0;
  while (group_index < sizeof(sat_groups) / sizeof(sat_groups[0])) {
    if (offset >= sat_groups[group_index].sat_count) {
      offset -= sat_groups[group_index].sat_count;
      group_index++;
    } else {
      sid = (gnss_signal_t) {
        .constellation = sat_groups[group_index].constellation,
        .band = sat_groups[group_index].band,
        .sat = sat_groups[group_index].sat_get(offset)
      };
      return sid;
    }
  }
  assert("sid_from_index() failed");
  return sid;
}

u32 sid_to_index(gnss_signal_t sid)
{
  assert(sid_valid(sid));
  u32 offset = 0;
  for (u32 i=0; i<sid.constellation; i++)
    offset += sat_groups[i].sat_count;
  return offset + sat_groups[sid.constellation].local_index_get(sid.sat);
}
