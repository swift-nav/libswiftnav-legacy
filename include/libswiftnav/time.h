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

#ifndef LIBSWIFTNAV_TIME_H
#define LIBSWIFTNAV_TIME_H

#include <time.h>

#include <libswiftnav/common.h>

/** Number of days in a common (non-leap) year. */
#define COMMON_YEAR_DAYS (365)

/** Number of days in a leap year. */
#define LEAP_YEAR_DAYS (366)

/** Number of days in four years. */
#define FOUR_YEARS_DAYS (3 * COMMON_YEAR_DAYS + LEAP_YEAR_DAYS)

/** Number of days in 100 years. */
#define HUNDRED_YEARS_DAYS (24 * FOUR_YEARS_DAYS + 4 * COMMON_YEAR_DAYS)

/** Number of days in 400 years. */
#define  FOUR_HUNDRED_YEARS_DAYS (3 * HUNDRED_YEARS_DAYS + 25 * FOUR_YEARS_DAYS)

/** Number of seconds in a minute. */
#define MINUTE_SECS (60)

/** Number of seconds in a hour. */
#define HOUR_SECS (60 * MINUTE_SECS)

/** Number of seconds in a day. */
#define DAY_SECS (24 * HOUR_SECS)

/** Number of seconds in a week. */
#define WEEK_SECS (7 * DAY_SECS)

/** Number of rollovers in the 10-bit broadcast GPS week number.
 * Update on next rollover on April 7, 2019.
 * \todo Detect and handle rollover more gracefully. */
#define GPS_WEEK_CYCLE 1

/** Offset between GPS and UTC times in seconds.
 * Update when a new leap second is inserted and be careful about times in the
 * past when this offset was different.
 * TODO handle leap seconds properly!
 */
#define GPS_MINUS_UTC_SECS 17

/** Unix timestamp of the GPS epoch 1980-01-06 00:00:00 UTC */
#define GPS_EPOCH 315964800

/** Modified Julian days of the GPS epoch 1980-01-06 00:00:00 UTC */
#define MJD_JAN_6_1980 44244

/** Modified Julian days of 1601-01-01 */
#define MJD_JAN_1_1601 -94187
//#define MJD_JAN_1_1601 15385

#define WN_UNKNOWN -1

/** Structure representing a GPS time. */
typedef struct __attribute__((packed)) {
  double tow; /**< Seconds since the GPS start of week. */
  s16 wn;     /**< GPS week number. */
} gps_time_t;

/** Structure containing GPS UTC correction parameters. */
typedef struct {
  double a0;
  double a1;
  gps_time_t tot;
  gps_time_t t_lse;
  u8 dt_ls;
  u8 dt_lsf;
} utc_params_t;

/** Structure representing UTC time. */
typedef struct {
  u16 year;           /**< Number of years AD. In four digit format. */
  u16 year_day;       /**< Day of the year (1 - 366). */
  u8 month;           /**< Month of the year (1 - 12). 0 = January, 12 = December. */
  u8 month_day;       /**< Day of the month (1 - 31). */
  u8 week_day;        /**< Day of the week (1 - 7). 1 = Monday, 7 = Sunday. */
  u8 hour;            /**< Minutes of the hour (0 - 59). */
  u8 minute;          /**< Minutes of the hour (0 - 59). */
  u8 second_int;      /**< Integer part of seconds of the minute (0 - 60). */
  double second_frac; /**< Fractional part of seconds (0 - .99...). */
} utc_time_t;

void normalize_gps_time(gps_time_t *t_gps);

void gps2utc(const utc_params_t *p, const gps_time_t *t, utc_time_t *u);
time_t gps2time(const utc_params_t *p, const gps_time_t *t_gps);

double gpsdifftime(const gps_time_t *end, const gps_time_t *beginning);
void gps_time_match_weeks(gps_time_t *t, const gps_time_t *ref);

#endif /* LIBSWIFTNAV_TIME_H */
