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

#include <math.h>

#include "gpstime.h"

/** Normalize a `gps_time_t` GPS time struct.
 * Ensures that the time of week is greater than zero and less than one week by
 * wrapping and adjusting the week number accordingly.
 *
 * \param t GPS time struct.
 * \return Normalized GPS time struct.
 */
/* TODO: Either normalise in place or rename to normalised. */
gps_time_t normalize_gps_time(gps_time_t t)
{
  while(t.tow < 0) {
    t.tow += 3600*24*7;
    t.wn += 1;
  }

  while(t.tow > 3600*24*7) {
    t.tow -= 3600*24*7;
    t.wn -= 1;
  }

  return t;
}

/** Convert a `gps_time_t` GPS time to a Unix `time_t`.
 * \note Adjusts for leap seconds using the current constant offset.
 *
 * \param t GPS time struct.
 * \return Unix time.
 */
time_t gps2time(gps_time_t gps_t)
{
  time_t t = GPS_EPOCH - GPS_MINUS_UTC_SECS;

  t += 7*24*3600*gps_t.wn;
  t += (s32)gps_t.tow;

  return t;
}

/** Time difference in seconds between two GPS times.
 * \param end Higher bound of the time interval whose length is calculated.
 * \param beginning Lower bound of the time interval whose length is
 *                  calculated. If this describes a time point later than end,
 *                  the result is negative.
 * \return The time difference in seconds between `beginning` and `end`.
 */
double gpsdifftime(gps_time_t end, gps_time_t beginning)
{
  return (end.wn - beginning.wn)*7*24*3600 + \
         end.tow - beginning.tow;
}


