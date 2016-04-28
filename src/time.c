/*
 * Copyright (C) 2013 Swift Navigation Inc.
 * Contact: Fergus Noble <fergus@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <libswiftnav/time.h>
#include <libswiftnav/constants.h>

/** \defgroup time Time functions
 * Functions to handle GPS and UTC time values.
 * \{ */

/** Normalize a `gps_time_t` GPS time struct.
 * Ensures that the time of week is greater than zero and less than one week by
 * wrapping and adjusting the week number accordingly.
 *
 * \param t GPS time struct.
 * \return Normalized GPS time struct.
 */
void normalize_gps_time(gps_time_t *t)
{
  while(t->tow < 0) {
    t->tow += WEEK_SECS;
    t->wn -= 1;
  }

  while(t->tow > WEEK_SECS) {
    t->tow -= WEEK_SECS;
    t->wn += 1;
  }
}

/** Convert a `gps_time_t` GPS time to a Unix `time_t`.
 * \note Adjusts for leap seconds using the current constant offset.
 *
 * \param t_gps GPS time struct.
 * \return Unix time.
 */
time_t gps2time(const gps_time_t *t_gps)
{
  time_t t = GPS_EPOCH - GPS_MINUS_UTC_SECS;

  t += WEEK_SECS * t_gps->wn;
  t += (s32)t_gps->tow;

  return t;
}

/** Time difference in seconds between two GPS times.
 * If the week number field of either time is unspecified, the result
 * will be as if the week numbers had been chosen for the times to be
 * as close as possible.
 * \param end Higher bound of the time interval whose length is calculated.
 * \param beginning Lower bound of the time interval whose length is
 *                  calculated. If this describes a time point later than end,
 *                  the result is negative.
 * \return The time difference in seconds between `beginning` and `end`.
 */
double gpsdifftime(const gps_time_t *end, const gps_time_t *beginning)
{
  double dt = end->tow - beginning->tow;
  if (end->wn == WN_UNKNOWN || beginning->wn == WN_UNKNOWN) {
    /* One or both of the week numbers is unspecified.  Assume times
       are within +/- 0.5 weeks of each other. */
    if (dt > WEEK_SECS / 2)
      dt -= WEEK_SECS;
    if (dt < -WEEK_SECS / 2)
      dt += WEEK_SECS;
  } else {
    /* Week numbers were provided - use them. */
    dt += (end->wn - beginning->wn) * WEEK_SECS;
  }
  return dt;
}

/** Given an uknown week number in t, fill in the week number from ref
 * under the assumption the two times are seperated by less than a week.
 *
 * \param t Pointer to GPS time whose week number will be set
 * \param ref Reference GPS time
 */
void gps_time_match_weeks(gps_time_t *t, const gps_time_t *ref)
{
  t->wn = ref->wn;
  double dt = t->tow - ref->tow;
  if (dt > WEEK_SECS / 2)
    t->wn--;
  else if (dt < -WEEK_SECS / 2)
    t->wn++;
}

/** Adjust the week number of t to correctly refelct the current week cycle.
 *
 * Assumes the current week number cannot be earlier than the reference WN. So
 * will return the correct WN for at most 20 years after the reference WN.
 *
 * \param wn_raw Raw week number from NAV data stream that is modulo 1024
 * \param wn_ref Reference week number that is from some point in the past
 * |return The abosolute week number counted from 1980
 */
u16 gps_adjust_week_cycle(u16 wn_raw, u16 wn_ref)
{
  if (wn_raw >= wn_ref)
  return wn_raw;

  return wn_raw + 1024 * ((wn_ref + 1023 - wn_raw) / 1024);
}

#define GLO_NT_0_FLOOR      1
#define GLO_NT_0_CEILING    LEAP_YEAR_DAYS                  /* 366 */
#define GLO_NT_1_CEILING    (GLO_NT_0_CEILING + YEAR_DAYS)  /* 731 */
#define GLO_NT_2_CEILING    (GLO_NT_1_CEILING + YEAR_DAYS)  /* 1096 */
#define GLO_NT_3_CEILING    (GLO_NT_2_CEILING + YEAR_DAYS)  /* 1461 */

#define N4_MIN              1
#define N4_MAX              31

/** Transformation of GLONASS-M current data information into gps_time_t.
 *
 *  Reference: GLONASS ICD Edition 5.1 2008
 *
 * \param nt current date, calendar number of day within four-year interval
 *           starting from the 1-st of January in a leap year
 * \param n4 four-year interval number starting from 1996
 * \param h  hours part of tk field
 * \param m  minutes part of tk field
 * \param s  seconds part of tk field
 * |return   converted gps time
 */

gps_time_t glo_time2gps_time(u16 nt, u8 n4, u8 h, u8 m, u8 s)
{
  u8 j = 0;
  u16 day_of_year = 0;
  u32 glo_year = 0;
  gps_time_t gps_t = GPS_TIME_UNKNOWN;

  if (n4 < N4_MIN || n4 > N4_MAX)
    return GPS_TIME_UNKNOWN;

  if (nt < GLO_NT_0_FLOOR)
    return GPS_TIME_UNKNOWN;
  else if (nt <= GLO_NT_0_CEILING) {
    j = 1;
    day_of_year = nt;
  }
  else if (nt <= GLO_NT_1_CEILING) {
    j = 2;
    day_of_year = nt - LEAP_YEAR_DAYS;
  }
  else if (nt <= GLO_NT_2_CEILING) {
    j = 3;
    day_of_year = nt - LEAP_YEAR_DAYS - YEAR_DAYS;
  }
  else if (nt <= GLO_NT_3_CEILING) {
    j = 4;
    day_of_year = nt - LEAP_YEAR_DAYS - YEAR_DAYS * 2;
  }
  else
    return GPS_TIME_UNKNOWN;

  glo_year = 1996 + 4 * (n4 - 1) + (j - 1);

  /* Calculate days since GPS epoch */
  s64 days_gps_epoch = YEAR_1980_GPS_DAYS + day_of_year - 1;
  u32 y = GPS_EPOCH_YEAR;
  while (++y < glo_year) {
    days_gps_epoch += is_leap_year(y) ? LEAP_YEAR_DAYS : YEAR_DAYS;
  }

  gps_t.wn = days_gps_epoch / WEEK_DAYS;
  gps_t.tow = GPS_MINUS_UTC_SECS + (days_gps_epoch % WEEK_DAYS) * DAY_SECS +
              (h - UTC_SU_OFFSET) * HOUR_SECS + m * MINUTE_SECS + s;

  normalize_gps_time(&gps_t);

  return gps_t;
}

/** \} */
