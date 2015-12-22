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

/** Number of seconds in a week. */
#define WEEK_SECS (7*24*60*60)

/** Number of seconds in a day. */
#define DAY_SECS (24*60*60)

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

/** Structure representing a GPS time. */
typedef struct __attribute__((packed)) {
  double tow; /**< Seconds since the GPS start of week. */
  s16 wn;     /**< GPS week number. */
} gps_time_t;

void normalize_gps_time(gps_time_t *t_gps);

time_t gps2time(const gps_time_t *t);

double gpsdifftime(const gps_time_t *end, const gps_time_t *beginning);
void gps_time_match_weeks(gps_time_t *t, const gps_time_t *ref);
u16 gps_adjust_week_cycle(u16 wn_raw, u16 wn_ref);

#endif /* LIBSWIFTNAV_TIME_H */
