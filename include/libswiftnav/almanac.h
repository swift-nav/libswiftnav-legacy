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

#ifndef LIBSWIFTNAV_ALMANAC_H
#define LIBSWIFTNAV_ALMANAC_H

#include <libswiftnav/common.h>
#include <libswiftnav/signal.h>
#include <libswiftnav/time.h>

/** \addtogroup almanac
 * \{ */

/** Structure containing the SBAS almanac for one satellite. */
typedef struct {
  double m0;       /**< Mean anomaly at reference time [semi-circles] */
  double ecc;      /**< Eccentricity. */
  double sqrta;    /**< Square root of the semi-major axis [sqrt(m)] */
  double omega0;   /**< Longitude of ascending node
                        of orbit plane at weekly epoch [semi-circles] */
  double omegadot; /**< Rate of right ascension [semi-circles/s] */
  double w;        /**< Argument of perigee [semi-circles] */
  double inc;      /**< Inclindation angle at reference time [semi-circles]
                        This must include the 0.3 offset from NAV delta_i. */
  double af0;      /**< Time offset of the sat clock [s] **/
  double af1;      /**< Drift of the sat clock [s/s] **/
} almanac_kepler_t;

/** Structure containing the SBAS almanac for one satellite. */
typedef struct {
  double pos[3]; /**< Position of the GEO at time toe [m] */
  double vel[3]; /**< velocity of the GEO at time toe [m/s] */
  double acc[3]; /**< velocity of the GEO at time toe [m/s^2] */
} almanac_xyz_t;

/** Structure containing the almanac for one satellite. */
typedef struct {
  gnss_signal_t sid; /**< Signal ID. */
  gps_time_t toa;    /**< Reference time of almanac. */
  float ura;         /**< User range accuracy [m] */
  u8 fit_interval;   /**< Curve fit interval. */
  u8 valid;          /**< Almanac is valid. */
  u8 healthy;        /**< Satellite health status. */
  union {
    almanac_kepler_t kepler; /**< Parameters specific to GPS. */
    almanac_xyz_t xyz;       /**< Parameters specific to SBAS. */
  };
} almanac_t;

/** \} */

s8 calc_sat_state_almanac(const almanac_t* a, const gps_time_t *t,
                            double pos[3], double vel[3],
                            double *clock_err, double *clock_rate_err);
s8 calc_sat_az_el_almanac(const almanac_t* a, const gps_time_t *t,
                          const double ref[3], double* az, double* el);
s8 calc_sat_doppler_almanac(const almanac_t* a, const gps_time_t *t,
                            const double ref[3], double *doppler);
u8 almanac_good(const almanac_t *a, gps_time_t *t);
bool almanac_equal(const almanac_t *a, const almanac_t *b);

#endif /* LIBSWIFTNAV_ALMANAC_H */
