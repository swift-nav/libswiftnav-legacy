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
#include <libswiftnav/constants.h>

/** \defgroup signal GNSS signal identifiers (SID)
 * \{ */

/** Element in the code data table. */
typedef struct {
  constellation_t constellation;
  u16 sat_count;
  u16 sat_start;
  const char *str;
} code_table_element_t;

/** Table of useful data for each code. */
static const code_table_element_t code_table[CODE_COUNT] = {
  [CODE_GPS_L1CA] =
    {CONSTELLATION_GPS, NUM_SIGNALS_GPS_L1CA, GPS_FIRST_PRN, "GPS L1CA"},
  [CODE_GPS_L2CM] =
    {CONSTELLATION_GPS, NUM_SIGNALS_GPS_L2CM, GPS_FIRST_PRN, "GPS L2CM"},
  [CODE_SBAS_L1CA] =
    {CONSTELLATION_SBAS, NUM_SIGNALS_SBAS_L1CA, SBAS_FIRST_PRN, "SBAS L1CA"},
};

/** String representation used for unknown code values. */
static const char * unknown_str = "?";

/** Construct a gnss_signal_t.
 *
 * \note This function does not check the validity of the resulting signal.
 *
 * \param code  Code to use.
 * \param sat   Satellite identifier to use.
 *
 * \return gnss_signal_t corresponding to the specified arguments.
 */
gnss_signal_t construct_sid(code_t code, u16 sat)
{
  gnss_signal_t sid = {.code = code, .sat = sat};
  return sid;
}

/** Print a string representation of a gnss_signal_t.
 *
 * \param s     Buffer of capacity n to which the string will be written.
 * \param n     Capacity of buffer s.
 * \param sid   gnss_signal_t to use.
 *
 * \return Number of characters written to s, excluding the terminating null.
 */
int sid_to_string(char *s, int n, gnss_signal_t sid)
{
  const char *code_str = ((sid.code < 0) || (sid.code >= CODE_COUNT)) ?
      unknown_str : code_table[sid.code].str;

  int nchars = snprintf(s, n, "%s %u", code_str, sid.sat);
  s[n-1] = 0;
  return nchars;
}

/** Determine if a gnss_signal_t corresponds to a known code and
 * satellite identifier.
 *
 * \param sid   gnss_signal_t to use.
 *
 * \return true if sid exists, false otherwise.
 */
bool sid_valid(gnss_signal_t sid)
{
  if (!code_valid(sid.code))
    return false;

  const code_table_element_t *e = &code_table[sid.code];
  if ((sid.sat < e->sat_start) || (sid.sat >= e->sat_start + e->sat_count))
    return false;

  return true;
}

/** Determine if a code is valid.
 *
 * \param code    Code to use.
 *
 * \return true if code is valid, false otherwise
 */
bool code_valid(code_t code)
{
  return ((code >= 0) && (code < CODE_COUNT));
}

/** Determine if a constellation is valid.
 *
 * \param constellation   Constellation to use.
 *
 * \return true if constellation is valid, false otherwise
 */
bool constellation_valid(constellation_t constellation)
{
  return ((constellation >= 0) && (constellation < CONSTELLATION_COUNT));
}

/** Convert a code-specific signal index to a gnss_signal_t.
 *
 * \param code          Code to use.
 * \param code_index    Code-specific signal index in
 *                      [0, SIGNAL_COUNT_\<code\>).
 *
 * \return gnss_signal_t corresponding to code and code_index.
 */
gnss_signal_t sid_from_code_index(code_t code, u16 code_index)
{
  assert(code_valid(code));
  assert(code_index < code_table[code].sat_count);
  return construct_sid(code, code_table[code].sat_start + code_index);
}

/** Return the code-specific signal index for a gnss_signal_t.
 *
 * \param sid   gnss_signal_t to use.
 *
 * \return Code-specific signal index in [0, SIGNAL_COUNT_\<code\>).
 */
u16 sid_to_code_index(gnss_signal_t sid)
{
  assert(sid_valid(sid));
  return sid.sat - code_table[sid.code].sat_start;
}

/** Get the constellation to which a gnss_signal_t belongs.
 *
 * \param sid   gnss_signal_t to use.
 *
 * \return Constellation to which sid belongs.
 */
constellation_t sid_to_constellation(gnss_signal_t sid)
{
  return code_to_constellation(sid.code);
}

/** Get the constellation to which a code belongs.
 *
 * \param code  Code to use.
 *
 * \return Constellation to which code belongs.
 */
constellation_t code_to_constellation(code_t code)
{
  assert(code_valid(code));
  return code_table[code].constellation;
}

/** Return the center carrier frequency for a gnss_signal_t.
 *
 * \param sid  gnss_signal_t to use.
 *
 * \return center carrier frequency
 */
double sid_to_carr_freq(gnss_signal_t code)
{
  double f;
  assert(sid_valid(code));
  switch (code.code) {
    case CODE_GPS_L1CA:
    case CODE_SBAS_L1CA:
      f = GPS_L1_HZ;
      break;
    case CODE_GPS_L2CM:
      f = GPS_L2_HZ;
      break;
    default:
      f = 0.0;
  }
  return f;
}

/* \} */
