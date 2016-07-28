/*
 * Copyright (C) 2010, 2016 Swift Navigation Inc.
 * Contact: Henry Hallam <henry@swift-nav.com>
 *          Fergus Noble <fergus@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef LIBSWIFTNAV_EPHEMERIS_H
#define LIBSWIFTNAV_EPHEMERIS_H

#include <libswiftnav/signal.h>
#include <libswiftnav/time.h>
#include <libswiftnav/common.h>

#define INVALID_GPS_URA_INDEX -1
#define INVALID_GPS_URA_VALUE -1.0f
#define MAX_ALLOWED_GPS_URA_IDX 15

/** \addtogroup ephemeris
 * \{ */

/** Structure containing the GPS ephemeris for one satellite. */
typedef struct {
  double tgd;      /**< Group delay between L1 and L2 [s] */
  double crc;      /**< Amplitude of the cosine harmonic correction term
                        to the orbit radius [m] */
  double crs;      /**< Amplitude of the sine harmonic correction term
                        to the orbit radius [m] */
  double cuc;      /**< Amplitude of the cosine harmonic correction term
                        to the argument of latitude [rad] */
  double cus;      /**< Amplitude of the sine harmonic correction term
                        to the argument of latitude [rad] */
  double cic;      /**< Amplitude of the cosine harmonic correction term
                        to the angle of inclination [rad] */
  double cis;      /**< Amplitude of the sine harmonic correction term
                        to the angle of inclination [rad] */
  double dn;       /**< Mean motion difference from computed value
                        [semi-circles/s] */
  double m0;       /**< Mean anomaly at reference time [semi-circles] */
  double ecc;      /**< Eccentricity. */
  double sqrta;    /**< Square root of the semi-major axis [sqrt(m)] */
  double omega0;   /**< Longitude of ascending node
                        of orbit plane at weekly epoch [semi-circles] */
  double omegadot; /**< Rate of right ascension [semi-circles/s] */
  double w;        /**< Argument of perigee [semi-circles] */
  double inc;      /**< Inclindation angle at reference time [semi-circles] */
  double inc_dot;  /**< Rate of inclination angle [semi-circles/s] */
  double af0;      /**< Time offset of the sat clock [s] **/
  double af1;      /**< Drift of the sat clock [s/s] **/
  double af2;      /**< Acceleration of the sat clock [s/s^2] **/
  gps_time_t toc;  /**< Reference time of clock. */
  u16 iodc;        /**< Issue of data clock. */
  u8 iode;         /**< Issue of data ephemeris. */
} ephemeris_kepler_t;

/** Structure containing the SBAS ephemeris for one satellite. */
typedef struct {
  double pos[3]; /**< Position of the GEO at time toe [m] */
  double vel[3]; /**< velocity of the GEO at time toe [m/s] */
  double acc[3]; /**< velocity of the GEO at time toe [m/s^2] */
  double a_gf0;  /**< Time offset of the GEO clock w.r.t. SNT [s] */
  double a_gf1;  /**< Drift of the GEO clock w.r.t. SNT [s/s] */
} ephemeris_xyz_t;

/** Structure containing the GLONASS ephemeris for one satellite. */
typedef struct {
  double gamma;     /**< Relative deviation of predicted carrier frequency
                         from nominal value, dimensionless */
  double tau;       /**< Correction to the SV time, seconds*/
  double pos[3];    /**< Position of the SV at tb in PZ-90.02 coordinates
                         system, meters */
  double vel[3];    /**< Velocity vector of the SV at tb in PZ-90.02
                         coordinates system, m/s */
  double acc[3];    /**< Acceleration vector of the SV at tb in PZ-90.02
                         coordinates system, m/s^2 */
} ephemeris_glo_t;

/** Structure containing the ephemeris for one satellite. */
typedef struct {
  gnss_signal_t sid; /**< Signal ID. */
  gps_time_t toe;    /**< Reference time of ephemeris. */
  float ura;         /**< User range accuracy [m] */
  u32 fit_interval;  /**< Curve fit interval [s] */
  u8 valid;          /**< Ephemeris is valid. */
  u8 health_bits;    /**< Satellite health status. */
  union {
    ephemeris_kepler_t kepler; /**< Parameters specific to GPS. */
    ephemeris_xyz_t xyz;       /**< Parameters specific to SBAS. */
    ephemeris_glo_t glo;       /**< Parameters specific to GLONASS. */
  };
} ephemeris_t;

/** \} */

s8 calc_sat_state(const ephemeris_t *e, const gps_time_t *t,
                  double pos[3], double vel[3],
                  double *clock_err, double *clock_rate_err);
s8 calc_sat_az_el(const ephemeris_t *e, const gps_time_t *t,
                  const double ref[3], double *az, double *el);
s8 calc_sat_doppler(const ephemeris_t *e, const gps_time_t *t,
                    const double ref[3], double *doppler);

u8 ephemeris_valid(const ephemeris_t *e, const gps_time_t *t);
u8 ephemeris_params_valid(const u8 valid, const u32 fit_interval,
                      const gps_time_t* toe, const gps_time_t *t);
u8 signal_healthy(const ephemeris_t *e, code_t code);

void decode_ephemeris(u32 frame_words[4][8], ephemeris_t *e);
bool ephemeris_equal(const ephemeris_t *a, const ephemeris_t *b);

float decode_ura_index(const u8 index);

#endif /* LIBSWIFTNAV_EPHEMERIS_H */
