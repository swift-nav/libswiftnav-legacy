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

/** \defgroup time Time functions
 * Functions to handle GPS and UTC time values.
 * \{ */

/** Normalize a `gps_time_t` GPS time struct in place.
 * Ensures that the time of week is greater than zero and less than one week by
 * wrapping and adjusting the week number accordingly.
 *
 * \param t GPS time struct.
 */
void normalize_gps_time(gps_time_t *t)
{
  while(t->tow < 0) {
    t->tow += WEEK_SECS;
    t->wn += 1;
  }

  while(t->tow > WEEK_SECS) {
    t->tow -= WEEK_SECS;
    t->wn -= 1;
  }
}

// http://www.ngs.noaa.gov/gps-toolbox/bwr-c.txt
// TODO note Beidou uses a different DN definition than GPS/Galileo
// TODO add day of week field to utc_time
void gps2utc(const utc_params_t *p, const gps_time_t *t, utc_time_t *u)
{
  /* Seconds from current time to leap second effctivity time */
  double dt_lse = gpsdifftime(t, &p->t_lse);

  /* Seconds from UTC data reference time (tot) */
  double dt = gpsdifftime(t, &p->tot); // TODO when making tot need to set WN correct

  // TODO move to header and python - split into hours?

  /* Check if we are near the LS effectivity time and pick correct LS offset */
  bool within_6hrs = false;
  u8 dt_ls;
  if ((dt_lse >= -6 * HOUR_SECS) && (dt_lse <= 6 * HOUR_SECS)) {
    /* We are within 6 hours before or after LS effectivity time */
    dt_ls = p->dt_ls;
    within_6hrs = true;

  } else if (dt_lse < 0.0) {
    /* We are before the LS effectivity time */
    dt_ls = p->dt_ls;

  } else {
    /* We are after the LS effectivity time */
    dt_ls = p->dt_lsf;
  }

  /* Calculate the UTC to GPS difference */
  double dt_utc = dt_ls + p->a0 + p->a1 * dt;

  /* Check if we need to correct for LS event */
  double t_utc;
  if (within_6hrs) {
    double w =
      fmod(t->tow - dt_utc - 12 * HOUR_SECS, DAY_SECS) + 12 * HOUR_SECS;
    t_utc = fmod(w, DAY_SECS + p->dt_lsf - p->dt_ls);
  } else {
    t_utc = fmod(t->tow - dt_utc, DAY_SECS);
  }

  /* t_utc now contains the seconds from the start of the day */
  /* It may be 1 second longer than a normal day due to LS */

  /* Convert this into hours, minutes and seconds */
  u32 second_int = t_utc;                /* The integer part of the seconds */
  u->second_frac = fmod(t_utc, 1.0);     /* The fractional part of the seconds */
  u->hour = second_int / HOUR_SECS;      /* The hours (1 - 23) */
  if (u->hour > 23) {
    u->hour = 23;                        /* Correct for a LS making it 24 */
  }
  second_int -= u->hour * HOUR_SECS;     /* Remove the hours from seconds */
  u->minute = second_int / MINUTE_SECS;  /* The minutes (1 - 59) */
  if (u->minute > 59) {
    u->minute = 59;                      /* Correct for a LS making it 60 */
  }
  second_int -= u->minute * MINUTE_SECS; /* Remove the minutes from seconds */
  u->second_int = second_int;            /* The seconds (1 - 60) */

  /* Calculate the years */

  /* Days from 1 Jan 1980. GPS epoch is 6 Jan 1980 */
  u8 day_of_week = t->tow / DAY_SECS;
  day_of_week += ((fmod(t->tow, DAY_SECS) > 23 * HOUR_SECS) &&
                 (t_utc < 1 * HOUR_SECS));
  u32 modified_julian_days = MJD_JAN_6_1980 + t->wn * 7 + day_of_week;
  u32 days_since_1601 = modified_julian_days - MJD_JAN_1_1601;

  /* Calculate the number of leap years */
  u16 num_400_years = days_since_1601 / FOUR_HUNDRED_YEARS_DAYS;
  u32 days_left = days_since_1601 -
                       num_400_years * FOUR_HUNDRED_YEARS_DAYS;
  u16 num_100_years = days_left / HUNDRED_YEARS_DAYS -
                                days_left / (FOUR_HUNDRED_YEARS_DAYS - 1);
  days_left -= num_100_years * HUNDRED_YEARS_DAYS;
  u16 num_4_years = days_left / FOUR_YEARS_DAYS;
  days_left -= num_4_years * FOUR_YEARS_DAYS;
  u16 num_non_leap_years = days_left / COMMON_YEAR_DAYS -
                                days_left / (FOUR_YEARS_DAYS - 1);

  /* Calculate the total number of years from 1980 */
  u->year = 1601 + num_400_years * 400 +
            num_100_years * 100 + num_4_years * 4 + num_non_leap_years;

  /* Calculate the month of the year */

  /* Calculate the day of the current year */
  u->year_day = days_left - num_non_leap_years * COMMON_YEAR_DAYS + 1;

  /* Work out if it is currently a leap year, 0 = no, 1 = yes` */
  u8 leap_year = ((u->year % 4) == 0) &&
                 (((u->year % 100) != 0) || ((u->year % 400) == 0));

  /* Roughly work out the month number */
  u8 month_guess = u->year_day * 0.032;

  /* Lookup table of the total number of days in the year after each month */
  /* First row is for a non-leap year, second row is for a leap year */
  static u16 days_after_month[2][13] = {
    {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365},
    {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366}
  };

  /* Check if our guess was out, and what the correction is, */
  /* 0 = correct, 1 = wrong */
  u8 month_correction = (u->year_day - days_after_month[leap_year][month_guess + 1]) > 0;

  /* Calculate the corrected number of months */
  u->month = month_guess + month_correction + 1;

  /* Calculate the day of the month */
  u->month_day = u->year_day - days_after_month[leap_year][month_guess + month_correction];

  /* Calculate the day of the week. 1 Jan 1601 was a Monday */
  u->week_day = days_since_1601 % 7 + 1;

  // TODO how use within_6hrs to handle roll over of leap second day past 60 seconds?
}

/** Convert a `gps_time_t` GPS time to a Unix `time_t`.
 * \note Adjusts for leap seconds using the current constant offset.
 *
 * \param p UTC correction parameters.
 * \param t_gps GPS time struct.
 * \return Unix time.
 */
time_t gps2time(const utc_params_t *p, const gps_time_t *t_gps)
{
  // TODO should go via gps2utc
  time_t t = GPS_EPOCH - p->dt_ls;

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

/** \} */
