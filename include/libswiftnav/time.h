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

/** Number of days in a year. */
#define YEAR_DAYS           365
#define LEAP_YEAR_DAYS      (YEAR_DAYS + 1)

/** Number of days in a week. */
#define WEEK_DAYS           7

/** Number of months in a year. */
#define YEAR_MONTHS         12

/** Days in (leap) year 1980 since GPS epoch Jan 6th */
#define YEAR_1980_GPS_DAYS  361

/** Year of GPS epoch */
#define GPS_EPOCH_YEAR      1980

/** UTC (SU) offset */
#define UTC_SU_OFFSET       3

/** Number of seconds in a minute. */
#define MINUTE_SECS         60

/** Number of minutes in an hour. */
#define HOUR_MINUTES        60

/** Number of seconds in an hour. */
#define HOUR_SECS           (MINUTE_SECS * HOUR_MINUTES)

/** Number of hours in a day. */
#define DAY_HOURS           24

/** Number of seconds in a day. */
#define DAY_SECS            (DAY_HOURS * HOUR_MINUTES * MINUTE_SECS)

/** Number of seconds in a week. */
#define WEEK_SECS           (WEEK_DAYS * DAY_SECS)

/** Number of rollovers in the 10-bit broadcast GPS week number.
 * Update on next rollover on April 7, 2019.
 * \todo Detect and handle rollover more gracefully. */
#define GPS_WEEK_CYCLE 1

/** The GPS week reference number. The current GPS WN is always assumed to be
 * later than this reference number. It will keep the WN calculated from the
 * truncated 10-bit broadcast WN valid for ~20 years after this week.
 *
 * Current week number is set to 20 December 2015. */
#define GPS_WEEK_REFERENCE 1876

/** Offset between GPS and UTC times in seconds.
 * Update when a new leap second is inserted and be careful about times in the
 * past when this offset was different.
 * TODO handle leap seconds properly!
 */
#define GPS_MINUS_UTC_SECS 17

/** Unix timestamp of the GPS epoch 1980-01-06 00:00:00 UTC */
#define GPS_EPOCH 315964800

#define WN_UNKNOWN -1
#define TOW_UNKNOWN -1

/** Structure representing a GPS time. */
typedef struct __attribute__((packed)) {
  double tow; /**< Seconds since the GPS start of week. */
  s16 wn;     /**< GPS week number. */
} gps_time_t;

#define GPS_TIME_UNKNOWN ((gps_time_t){TOW_UNKNOWN, WN_UNKNOWN})

void normalize_gps_time(gps_time_t *t_gps);

time_t gps2time(const gps_time_t *t);

double gpsdifftime(const gps_time_t *end, const gps_time_t *beginning);
void gps_time_match_weeks(gps_time_t *t, const gps_time_t *ref);
u16 gps_adjust_week_cycle(u16 wn_raw, u16 wn_ref);

static inline bool is_leap_year(s32 year)
{
  return ((year%4==0) && (year%100!=0)) || (year%400==0);
}

gps_time_t glo_time2gps_time(u16 nt, u8 n4, s8 h, s8 m, s8 s);

#endif /* LIBSWIFTNAV_TIME_H */
