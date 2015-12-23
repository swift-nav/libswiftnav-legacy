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
#include <assert.h>
#include <memory.h>

#include <libswiftnav/constants.h>
#include <libswiftnav/linear_algebra.h>
#include <libswiftnav/coord_system.h>
#include <libswiftnav/time.h>
#include <libswiftnav/almanac.h>
#include <libswiftnav/ephemeris.h>

/** \defgroup almanac Almanac
 * Functions and calculations related to the GPS almanac.
 *
 * \note All positions are referenced to the WGS84 coordinate system.
 * \see coord_system
 * \{ */

/** Calculate satellite position, velocity and clock offset from SBAS ephemeris.
 *
 * \param a Pointer to an almanac structure for the satellite of interest
 * \param t GPS time at which to calculate the satellite state
 * \param pos Array into which to write calculated satellite position [m]
 * \param vel Array into which to write calculated satellite velocity [m/s]
 * \param clock_err Pointer to where to store the calculated satellite clock
 *                  error [s]
 * \param clock_rate_err Pointer to where to store the calculated satellite
 *                       clock error [s/s]
 */
static s8 calc_sat_state_xyz_almanac(const almanac_t* a, const gps_time_t *t,
                                     double pos[3], double vel[3],
                                     double *clock_err, double *clock_rate_err)
{
  ephemeris_t e;
  memset(&e, 0, sizeof(e));
  e.sid = a->sid;
  e.toe = a->toa;
  e.ura = a->ura;
  e.fit_interval = a->fit_interval;
  e.valid = a->valid;
  e.healthy = a->healthy;
  memcpy(e.xyz.pos, a->xyz.pos, sizeof(e.xyz.pos));
  memcpy(e.xyz.vel, a->xyz.vel, sizeof(e.xyz.vel));
  memcpy(e.xyz.acc, a->xyz.acc, sizeof(e.xyz.acc));

  return calc_sat_state(&e, t, pos, vel, clock_err, clock_rate_err);
}

/** Calculate satellite position, velocity and clock offset from GPS ephemeris.
 *
 * References:
 *   -# IS-GPS-200D, Section 20.3.3.5.2.1 and Table 20-VI
 *
 * \param a Pointer to an almanac structure for the satellite of interest
 * \param t GPS time at which to calculate the satellite state
 * \param pos Array into which to write calculated satellite position [m]
 * \param vel Array into which to write calculated satellite velocity [m/s]
 * \param clock_err Pointer to where to store the calculated satellite clock
 *                  error [s]
 * \param clock_rate_err Pointer to where to store the calculated satellite
 *                       clock error [s/s]
 *
 * \return  0 on success,
 *         -1 if almanac is not valid or too old
 */
static s8 calc_sat_state_kepler_almanac(const almanac_t* a, const gps_time_t *t,
                                          double pos[3], double vel[3],
                                          double *clock_err, double *clock_rate_err)
{
  ephemeris_t e;
  memset(&e, 0, sizeof(e));
  e.sid = a->sid;
  e.toe = a->toa;
  e.ura = a->ura;
  e.fit_interval = a->fit_interval;
  e.valid = a->valid;
  e.healthy = a->healthy;
  e.kepler.m0 = a->kepler.m0;
  e.kepler.ecc = a->kepler.ecc;
  e.kepler.sqrta = a->kepler.sqrta;
  e.kepler.omega0 = a->kepler.omega0;
  e.kepler.omegadot = a->kepler.omegadot;
  e.kepler.w = a->kepler.w;
  e.kepler.inc = a->kepler.inc;
  e.kepler.af0 = a->kepler.af0;
  e.kepler.af1 = a->kepler.af1;

  return calc_sat_state(&e, t, pos, vel, clock_err, clock_rate_err);
}

/** Calculate satellite position, velocity and clock offset from almanac.
 *
 * Dispatch to internal function for Kepler/XYZ almanac depending on
 * constellation.
 *
 * \param a Pointer to an almanac structure for the satellite of interest
 * \param t GPS time at which to calculate the satellite state
 * \param pos Array into which to write calculated satellite position [m]
 * \param vel Array into which to write calculated satellite velocity [m/s]
 * \param clock_err Pointer to where to store the calculated satellite clock
 *                  error [s]
 * \param clock_rate_err Pointer to where to store the calculated satellite
 *                       clock error [s/s]
 *
 * \return  0 on success,
 *         -1 if almanac is not valid or too old
 */
s8 calc_sat_state_almanac(const almanac_t* a, const gps_time_t *t,
                            double pos[3], double vel[3],
                            double *clock_err, double *clock_rate_err)
{
  switch(sid_to_constellation(a->sid)) {
  case CONSTELLATION_GPS:
    return calc_sat_state_kepler_almanac(a, t, pos, vel, clock_err, clock_rate_err);
  case CONSTELLATION_SBAS:
    return calc_sat_state_xyz_almanac(a, t, pos, vel, clock_err, clock_rate_err);
    break;
  default:
    assert(!"Unsupported constellation");
    return -1;
  }
}
// TODO why are there not equiv func in ephemeris
/** Calculate the azimuth and elevation of a satellite from a reference
 * position given the satellite almanac.
 *
 * \param a  Pointer to an almanac structure for the satellite of interest.
 * \param t    GPS time at which to calculate the az/el.
 * \param ref  ECEF coordinates of the reference point from which the azimuth
 *             and elevation is to be determined, passed as [X, Y, Z], all in
 *             meters.
 * \param az   Pointer to where to store the calculated azimuth output [rad].
 * \param el   Pointer to where to store the calculated elevation output [rad].
 * \return  0 on success,
 *         -1 if almanac is not valid or too old
 */
s8 calc_sat_az_el_almanac(const almanac_t* a, const gps_time_t *t,
                          const double ref[3], double* az, double* el)
{
  double sat_pos[3];
  double sat_vel[3];
  double clock_err, clock_rate_err;
  s8 ret = calc_sat_state_almanac(a, t, sat_pos, sat_vel, &clock_err, &clock_rate_err);
  if (ret != 0) {
    return ret;
  }
  wgsecef2azel(sat_pos, ref, az, el);

  return 0;
}

/** Calculate the Doppler shift of a satellite as observed at a reference
 * position given the satellite almanac.
 *
 * \param a  Pointer to an almanac structure for the satellite of interest.
 * \param t    GPS time at which to calculate the az/el.
 * \param ref  ECEF coordinates of the reference point from which the
 *             Doppler is to be determined, passed as [X, Y, Z], all in
 *             meters.
 * \param doppler The Doppler shift [Hz].
 * \return  0 on success,
 *         -1 if almanac is not valid or too old
 */
s8 calc_sat_doppler_almanac(const almanac_t* a, const gps_time_t *t,
                            const double ref[3], double *doppler)
{
  double sat_pos[3];
  double sat_vel[3];
  double clock_err, clock_rate_err;
  double vec_ref_sat[3];

  s8 ret = calc_sat_state_almanac(a, t, sat_pos, sat_vel, &clock_err, &clock_rate_err);
  if (ret != 0) {
    return ret;
  }

  /* Find the vector from the reference position to the satellite. */
  vector_subtract(3, sat_pos, ref, vec_ref_sat);

  /* Find the satellite velocity projected on the line of sight vector from the
   * reference position to the satellite. */
  double radial_velocity = vector_dot(3, vec_ref_sat, sat_vel) / \
                           vector_norm(3, vec_ref_sat);

  /* Return the Doppler shift. */
  *doppler = GPS_L1_HZ * radial_velocity / GPS_C;

  return 0;
}

/** Is this almanac usable?
 *
 * \todo This should actually be more than just the "valid" flag.
 *       When we write an is_usable() function, lets use that instead
 *       of just es[prn].valid.
 *
 * \param a Almanac struct
 * \param t Current time
 * \return 1 if the ephemeris is valid and not too old.
 *         0 otherwise.
 */
u8 almanac_good(const almanac_t *a, gps_time_t *t)
{
  /* Seconds from the time from almanac reference epoch (toe) */
  double dt = gpsdifftime(t, &a->toa);

  /* TODO: this doesn't exclude almanacs older than a week so could be made
   * better. */
  return (a->valid && fabs(dt) < ((u32)a->fit_interval)*60*60);
}

static bool almanac_xyz_equal(const almanac_xyz_t *a,
                              const almanac_xyz_t *b)
{
  return (memcmp(a->pos, b->pos, sizeof(a->pos)) == 0) &&
         (memcmp(a->vel, b->vel, sizeof(a->vel)) == 0) &&
         (memcmp(a->acc, b->acc, sizeof(a->acc)) == 0);
}

static bool almanac_kepler_equal(const almanac_kepler_t *a,
                                 const almanac_kepler_t *b)
{
  return (a->m0 == b->m0) &&
         (a->ecc == b->ecc) &&
         (a->sqrta == b->sqrta) &&
         (a->omega0 == b->omega0) &&
         (a->omegadot == b->omegadot) &&
         (a->w == b->w) &&
         (a->inc == b->inc) &&
         (a->af0 == b->af0) &&
         (a->af1 == b->af1);
}

/** Are the two almanacs the same?
 *
 * \param a First almanac
 * \param b Second almanac
 * \return true if they are equal
 */
bool almanac_equal(const almanac_t *a, const almanac_t *b)
{
  if (!sid_is_equal(a->sid, b->sid) ||
      (a->ura != b->ura) ||
      (a->fit_interval != b->fit_interval) ||
      (a->valid != b->valid) ||
      (a->healthy != b->healthy) ||
      (a->toa.wn != b->toa.wn) ||
      (a->toa.tow != b->toa.tow))
    return false;

  switch (sid_to_constellation(a->sid)) {
  case CONSTELLATION_GPS:
    return almanac_kepler_equal(&a->kepler, &b->kepler);
  case CONSTELLATION_SBAS:
    return almanac_xyz_equal(&a->xyz, &b->xyz);
  default:
    assert(!"Unsupported constellation");
    return false;
  }
}

/** \} */
