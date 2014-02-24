/*
 * Copyright (C) 2010 Swift Navigation Inc.
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

#include "constants.h"
#include "linear_algebra.h"
#include "coord_system.h"
#include "almanac.h"

/** \defgroup almanac Almanac
 * Functions and calculations related to the GPS almanac.
 *
 * \note All positions are referenced to the WGS84 coordinate system.
 * \see coord_system
 * \{ */

/** Calculate the position / velocity state of a satellite from the almanac.
 *
 * \param alm  Pointer to an almanac structure for the satellite of interest.
 * \param t    GPS time of week at which to calculate the satellite state in
 *             seconds since Sunday.
 * \param week GPS week number modulo 1024 or pass -1 to assume within one
 *             half-week of the almanac time of applicability.
 * \param pos  The satellite position in ECEF coordinates is returned in this
 *             vector.
 * \param vel  The satellite velocity in ECEF coordinates is returned in this
 *             vector. Ignored if NULL.
 */
void calc_sat_state_almanac(almanac_t* alm, double t, s16 week,
                            double pos[3], double vel[3])
{
  /* Seconds since the almanac reference epoch. */
  double dt = t - alm->toa;

  if (week < 0) {
    /* Week number unknown, correct time for beginning or end of week
     * crossovers and limit to +/- 302400 (i.e. assume dt is within a
     * half-week). */
    if (dt > 302400)
      dt -= 604800;
    else if (dt < -302400)
      dt += 604800;
  } else {
    /* Week number specified, correct time using week number difference. */
    s32 dweeks = week - alm->week;
    dt += dweeks * 604800;
  }

  /* Calculate position and velocity per ICD-GPS-200D Table 20-IV. */

  /* Calculate mean motion in radians/sec. */
  double ma_dot = sqrt (GPS_GM / (alm->a * alm->a * alm->a));
  /* Calculate corrected mean anomaly in radians. */
  double ma = alm->ma + ma_dot * dt;

  /* Iteratively solve for the Eccentric Anomaly
   * (from Keith Alter and David Johnston). */
  double ea = ma;  /* Starting value for E. */
  double ea_old;
  double temp;
  double ecc = alm->ecc;
  u32 count = 0;

  /* TODO: Implement convergence test using integer difference of doubles,
   * http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm */
  do {
    ea_old = ea;
    temp = 1.0 - ecc * cos(ea_old);
    ea = ea + (ma - ea_old + ecc * sin(ea_old)) / temp;
    count++;
    if (count > 5)
      break;
  } while (fabs(ea - ea_old) > 1.0e-14);

  double ea_dot = ma_dot / temp;

  /* Begin calculation for True Anomaly and Argument of Latitude. */
  double temp2 = sqrt(1.0 - ecc * ecc);
  /* Argument of Latitude = True Anomaly + Argument of Perigee. */
  double al = atan2(temp2 * sin(ea), cos(ea) - ecc) + alm->argp;
  double al_dot = temp2 * ea_dot / temp;

  /* Calculate corrected radius based on argument of latitude. */
  double r = alm->a * temp;
  double r_dot = alm->a * ecc * sin(ea) * ea_dot;

  /* Calculate position and velocity in orbital plane. */
  double x = r * cos(al);
  double y = r * sin(al);
  double x_dot = r_dot * cos(al) - y * al_dot;
  double y_dot = r_dot * sin(al) + x * al_dot;

  /* Corrected longitude of ascending node. */
  double om_dot = alm->rora - GPS_OMEGAE_DOT;
  double om = alm->raaw + dt * om_dot - GPS_OMEGAE_DOT * alm->toa;

  /* Compute the satellite's position in Earth-Centered Earth-Fixed
   * coordinates. */
  pos[0] = x * cos(om) - y * cos(alm->inc) * sin(om);
  pos[1] = x * sin(om) + y * cos(alm->inc) * cos(om);
  pos[2] = y * sin(alm->inc);

  /* Compute the satellite's velocity in Earth-Centered Earth-Fixed
   * coordinates. */
  if (vel) {
    temp = y_dot * cos(alm->inc);
    vel[0] = -om_dot * pos[1] + x_dot * cos(om) - temp * sin(om);
    vel[1] =  om_dot * pos[0] + x_dot * sin(om) + temp * cos(om);
    vel[2] = y_dot * sin(alm->inc);
  }

}

/** Calculate the azimuth and elevation of a satellite from a reference
 * position given the satellite almanac.
 *
 * \param alm  Pointer to an almanac structure for the satellite of interest.
 * \param t    GPS time of week at which to calculate the az/el.
 * \param week GPS week number modulo 1024 or pass -1 to assume within one
 *             half-week of the almanac time of applicability.
 * \param ref  ECEF coordinates of the reference point from which the azimuth
 *             and elevation is to be determined, passed as [X, Y, Z], all in
 *             meters.
 * \param az   Pointer to where to store the calculated azimuth output.
 * \param el   Pointer to where to store the calculated elevation output.
 */
void calc_sat_az_el_almanac(almanac_t* alm, double t, s16 week,
                            double ref[3], double* az, double* el)
{
  double sat_pos[3];
  calc_sat_state_almanac(alm, t, week, sat_pos, 0);
  wgsecef2azel(sat_pos, ref, az, el);
}

/** Calculate the Doppler shift of a satellite as observed at a reference
 * position given the satellite almanac.
 *
 * \param alm  Pointer to an almanac structure for the satellite of interest.
 * \param t    GPS time of week at which to calculate the Doppler shift.
 * \param week GPS week number modulo 1024 or pass -1 to assume within one
 *             half-week of the almanac time of applicability.
 * \param ref  ECEF coordinates of the reference point from which the azimuth
 *             and elevation is to be determined, passed as [X, Y, Z], all in
 *             meters.
 * \return     The Doppler shift in Hz.
 */
double calc_sat_doppler_almanac(almanac_t* alm, double t, s16 week,
                                double ref[3])
{
  double sat_pos[3];
  double sat_vel[3];
  double vec_ref_sat[3];

  calc_sat_state_almanac(alm, t, week, sat_pos, sat_vel);

  /* Find the vector from the reference position to the satellite. */
  vector_subtract(3, sat_pos, ref, vec_ref_sat);

  /* Find the satellite velocity projected on the line of sight vector from the
   * reference position to the satellite. */
  double radial_velocity = vector_dot(3, vec_ref_sat, sat_vel) / \
                           vector_norm(3, vec_ref_sat);

  /* Return the Doppler shift. */
  return GPS_L1_HZ * radial_velocity / GPS_C;
}

/** \} */

