/*
 * Copyright (c) 2016 Swift Navigation Inc.
 * Contact: Jacob McNamee <jacob@swiftnav.com>
 *          Pasi Miettinen <pasi.miettinen@exafore.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#include <stdio.h>
#include <assert.h>

#include <libswiftnav/signal.h>

/*
The implementation assumes that we have the following virtual list of codes:

     prn         code
 0: {1, CODE_GPS_L1CA},
 1: {1, CODE_GPS_L2CM},
 2: {2, CODE_GPS_L1CA},
 3: {2, CODE_GPS_L2CM},
...
81: {137, CODE_SBAS_L1CA},
82: {138, CODE_SBAS_L1CA},
*/

static const u16 number_of_sats[] = {
  NUM_SATS_GPS,
  NUM_SATS_SBAS,
};

static const u16 first_sat_index[] = {
  0,              /* GPS  */
  NUM_SIGNALS_GPS,  /* SBAS */
};

static const u16 first_sat_id[] = {
  GPS_FIRST_PRN,     /* GPS  */
  SBAS_FIRST_PRN,    /* SBAS */
};

static const u8 codes_on_cons[] = {
  NUM_CODES_GPS,   /* GPS  */
  NUM_CODES_SBAS, /* SBAS */
};

static const u8 codes_on_gps[NUM_CODES_GPS] = {
  CODE_GPS_L1CA,
  CODE_GPS_L2CM,
};

static const u8 codes_on_sbas[NUM_CODES_SBAS] = {
  CODE_SBAS_L1CA,
};

static const u8 local_code_indexes[CODE_COUNT] = {
  [CODE_GPS_L1CA] = 0,
  [CODE_GPS_L2CM] = 1,
  [CODE_SBAS_L1CA] = 0,
};

static const char * code_strs[CODE_COUNT] = {
  [CODE_GPS_L1CA] = "GPS L1CA",
  [CODE_GPS_L2CM] = "GPS L2CM",
  [CODE_SBAS_L1CA] = "SBAS L1CA"
};

static const char * unknown_str = "?";

gnss_signal_t construct_sid(enum code code, u16 sat)
{
  gnss_signal_t sid = { .code = code, .sat = sat};
  return sid;
}

int sid_to_string(char *s, int n, gnss_signal_t sid)
{
  const char *code_str = (sid.code >= CODE_COUNT || sid.code == CODE_INVALID) ?
      unknown_str : code_strs[sid.code];

  int nchars = snprintf(s, n, "%s %u", code_str, sid.sat);
  s[n-1] = 0;
  return nchars;
}

bool sid_valid(gnss_signal_t sid)
{
  if (sid_to_constellation(sid) == CONSTELLATION_INVALID) {
    return false;
  }

  if (sid.sat < first_sat_id[sid_to_constellation(sid)] ||
      sid.sat >= first_sat_id[sid_to_constellation(sid)] +
      number_of_sats[sid_to_constellation(sid)])
    return false;

  return true;
}

gnss_signal_t sid_from_index(u16 idx)
{
  gnss_signal_t sid;
  u8 code_idx = 0;
  enum constellation cons = CONSTELLATION_INVALID;

  if (idx < first_sat_index[CONSTELLATION_GPS] +
      NUM_SATS_GPS * codes_on_cons[CONSTELLATION_GPS]) {
    cons = CONSTELLATION_GPS;
  } else if (idx < first_sat_index[CONSTELLATION_SBAS] +
      NUM_SATS_SBAS * codes_on_cons[CONSTELLATION_SBAS]) {
    cons = CONSTELLATION_SBAS;
  } else {
    assert("sid_from_index() invalid index");
  }

  code_idx = (idx - first_sat_index[cons]) % codes_on_cons[cons];

  switch (cons) {
  case CONSTELLATION_GPS:
    sid.code = codes_on_gps[code_idx];
    break;
  case CONSTELLATION_SBAS:
    if (sid.code != CODE_SBAS_L1CA)
    sid.code = codes_on_sbas[code_idx];
    break;
  default:
    break;
  }

  sid.sat = (idx - first_sat_index[cons]) /
      codes_on_cons[cons] + first_sat_id[cons];

  return sid;
}

/******************************************************************************
  sid    INDEX_GLOBAL    INDEX_CONSTELLATION I  NDEX_SAT_IN_CONS
sat code
 1  L1CA      0                   0                     0
 1  L2CM      1                   1                     0
 2  L1CA      2                   2                     1
 2  L2CM      3                   3                     1
...  ...    ...                 ...                   ...
 31 L1CA     60                  60                    30
 31 L2CM     61                  61                    30
 32 L1CA     62                  62                    31
 32 L2CM     63                  63                    31
120 L1CA     64                   0                     0
121 L1CA     65                   1                     1
...  ...    ...                 ...                   ...
137 L1CA     81                  17                    17
138 L1CA     82                  18                    18

'sat': Satellite ID in its constellation
'code':  A combination of constellation, band and signal encoding
'sid': Signal ID
*******************************************************************************/

u16 sid_to_index(gnss_signal_t sid, enum indexing_type it)
{
  enum constellation cons = sid_to_constellation(sid);

  switch (it) {
  case INDEX_GLOBAL:
    return first_sat_index[cons] +
        (sid.sat - first_sat_id[cons]) * codes_on_cons[cons] +
        local_code_indexes[sid.code];
    break;
  case INDEX_CONSTELLATION:
    return sid_to_index(sid, INDEX_GLOBAL) - first_sat_index[cons];
    break;
  case INDEX_SAT_IN_CONS:
    return sid.sat - first_sat_id[cons];
    break;
  default:
    assert("unknown indexing type");
    break;
  }

  return INVALID_SID_INDEX;
}

enum constellation sid_to_constellation(gnss_signal_t sid)
{
  switch (sid.code) {
  case CODE_GPS_L1CA:
  case CODE_GPS_L2CM:
    return CONSTELLATION_GPS;
    break;
  case CODE_SBAS_L1CA:
    return CONSTELLATION_SBAS;
    break;
  default:
    break;
  }
  return CONSTELLATION_INVALID;
}
