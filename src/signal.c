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
    {CONSTELLATION_GPS, NUM_SIGNALS_GPS_L1CA, GPS_FIRST_PRN,
     "GPS L1CA"},
  [CODE_GPS_L2CM] =
    {CONSTELLATION_GPS, NUM_SIGNALS_GPS_L2CM, GPS_FIRST_PRN,
     "GPS L2CM"},
  [CODE_SBAS_L1CA] =
    {CONSTELLATION_SBAS, NUM_SIGNALS_SBAS_L1CA, SBAS_FIRST_PRN,
     "SBAS L1CA"},
  [CODE_GLO_L1CA] =
    {CONSTELLATION_GLO, NUM_SIGNALS_GLO_L1CA, GLO_FIRST_PRN,
     "GLO L1CA"},
  [CODE_GLO_L2CA] =
    {CONSTELLATION_GLO, NUM_SIGNALS_GLO_L2CA, GLO_FIRST_PRN,
     "GLO L2CA"},
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
  assert(n >= SID_STR_LEN_MAX);
  const char *code_str = ((sid.code < 0) || (sid.code >= CODE_COUNT)) ?
      unknown_str : code_table[sid.code].str;
  int nchars = 0;

  /* Copy code string */
  for (u32 i=0; code_str[i] != 0; i++) {
    s[nchars++] = code_str[i];
  }

  s[nchars++] = ' ';

  /* Print sat value */
  bool started = false;
  u16 div = 10000;
  while (div > 0) {
    u8 digit = (sid.sat / div) % 10;
    div /= 10;
    if (started || (digit != 0)) {
      s[nchars++] = '0' + digit;
      started = true;
    }
  }

  s[nchars] = 0;
  assert(nchars < SID_STR_LEN_MAX);
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

/** Return the center carrier frequency for a code_t.
 *
 * \param code  code_t to use.
 * \return center carrier frequency
 */
double code_to_carr_freq(code_t code)
{
  double f;
  assert(code_valid(code));
  switch (code) {
  case CODE_GPS_L1CA:
  case CODE_SBAS_L1CA:
    f = GPS_L1_HZ;
    break;
  case CODE_GPS_L2CM:
    f = GPS_L2_HZ;
    break;
  default:
    assert(!"Unsupported code_t ID");
    f = 0.0;
    break;
  }
  return f;
}

/** Return the chips count for a code_t.
 *
 * \param code  code_t to use.
 * \return chips count
 */
u16 code_to_chip_count(code_t code)
{
  u16 cn;
  assert(code_valid(code));
  switch (code) {
  case CODE_GPS_L1CA:
  case CODE_SBAS_L1CA:
    cn = GPS_L1CA_CHIPS_NUM;
    break;
  case CODE_GPS_L2CM:
    cn = 2 * GPS_L2CM_CHIPS_NUM;
    break;
  default:
    assert(!"Unsupported code_t ID");
    cn = 0;
    break;
  }
  return cn;
}

/** Return the chip rate for a code_t.
 *
 * \param code  code_t to use.
 * \return chip rate
 */
double code_to_chip_rate(code_t code)
{
  double cr;
  assert(code_valid(code));
  switch (code) {
  case CODE_GPS_L1CA:
  case CODE_SBAS_L1CA:
  case CODE_GPS_L2CM:
    cr = GPS_CA_CHIPPING_RATE;
    break;
  default:
    assert(!"Unsupported code_t ID");
    cr = 0;
    break;
  }
  return cr;
}

/* \} */
