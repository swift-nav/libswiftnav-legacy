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

#include <libswiftnav/time.h>
#include <libswiftnav/constants.h>

// TODO add a doc group

/* TODO: does it make sense to be passing structs by value in all
   these functions? */

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
    t.tow += WEEK_SECS;
    t.wn -= 1;
  }

  while(t.tow > WEEK_SECS) {
    t.tow -= WEEK_SECS;
    t.wn += 1;
  }

  return t;
}

/** Convert a `gps_time_t` GPS time to a Unix `time_t`.
 * \note Adjusts for leap seconds using the current constant offset.
 *
 * \param gps_t GPS time struct.
 * \return Unix time.
 */
time_t gps2time(gps_time_t gps_t)
{
  time_t t = GPS_EPOCH - GPS_MINUS_UTC_SECS;

  t += WEEK_SECS*gps_t.wn;
  t += (s32)gps_t.tow;

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
double gpsdifftime(gps_time_t end, gps_time_t beginning)
{
  double dt = end.tow - beginning.tow;
  if (end.wn == WN_UNKNOWN || beginning.wn == WN_UNKNOWN) {
    /* One or both of the week numbers is unspecified.  Assume times
       are within +/- 0.5 weeks of each other. */
    if (dt > WEEK_SECS / 2)
      dt -= WEEK_SECS;
    if (dt < -WEEK_SECS / 2)
      dt += WEEK_SECS;
  } else {
    /* Week numbers were provided - use them. */
    dt += (end.wn - beginning.wn) * WEEK_SECS;
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
