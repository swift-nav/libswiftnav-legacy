/*
 * Copyright (C) 2015 Swift Navigation Inc.
 * Contact: Leith Bade <leith@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#include <float.h>
#include <math.h>
#include <stdio.h>

#include <libswiftnav/constants.h>
#include <libswiftnav/ionosphere.h>

/** \defgroup ionosphere Ionospheric models
 * Implemenations of ionoshperic delay correction models.
 * \{ */

/** Calculate ionospheric delay using Klobuchar model.
 *
 * References:
 *   -# IS-GPS-200H, Section 20.3.3.5.2.5 and Figure 20-4
 *
 * \param t_gps GPS time at which to calculate the ionospheric delay
 * \param lat_u Latitude of the receiver [rad]
 * \param lon_u Longitude of the receiver [rad]
 * \param a Azimuth of the satellite, clockwise positive from North [rad]
 * \param e Elevation of the satellite [rad]
 * \param i Ionosphere parameters struct from GPS NAV data
 *
 * \return  Ionospheric delay distance for GPS L1 frequency [m]
 */
double calc_ionosphere(const gps_time_t *t_gps,
                       double lat_u, double lon_u,
                       double a, double e,
                       const ionosphere_t *i)
{
  /* Convert inputs from radians to semicircles */
  /* All calculations are in semicircles */
  lat_u = lat_u / M_PI;
  lon_u = lon_u / M_PI;
  /* a can remain in radians */
  e = e / M_PI;

  /* Calculate the earth-centered angle */
  double psi = 0.0137 / (e + 0.11) - 0.022;

  /* Compute the latitude of the Ionospheric Pierce Point */
  double lat_i = lat_u + psi * cos(a);
  if (lat_i > 0.416) {
    lat_i = 0.416;
  }
  if (lat_i < -0.416) {
    lat_i = -0.416;
  }

  /* Compute the longitude of the IPP */
  double lon_i = lon_u + (psi * sin(a)) / cos(lat_i * M_PI);

  /* Find the geomagnetic latitude of the IPP */
  double lat_m = lat_i + 0.064 * cos((lon_i - 1.617) * M_PI);

  /* Find the local time at the IPP */
  double t = 43200.0 * lon_i + t_gps->tow;
  t = fmod(t, DAY_SECS);
  if (t > DAY_SECS) {
    t -= DAY_SECS;
  }
  if (t < 0.0) {
    t += DAY_SECS;
  }

  /* Compute the amplitude of ionospheric delay */
  double amp = i->a0 + lat_m * (i->a1 + lat_m * (i->a2 + i->a3 * lat_m));
  if (amp < 0.0) {
    amp = 0.0;
  }

  /* Compute the period of ionospheric delay */
  double per = i->b0 + lat_m * (i->b1 + lat_m * (i->b2 + i->b3 * lat_m));
  if (per < 72000.0) {
    per = 72000.0;
  }

  /* Compute the phase of ionospheric delay */
  double x = 2.0 * M_PI * (t - 50400.0) / per;

  /* Compute the slant factor */
  double temp = 0.53 - e;
  double sf = 1.0 + 16.0 * temp * temp * temp;

  /* Compute the ionospheric time delay */
  double d_l1;
  if (x >= 1.57) {
    d_l1 = sf * 5e-9;
  } else {
    double x_2 = x * x;
    d_l1 = sf * (5e-9 + amp * (1.0 - x_2 / 2.0 + x_2 * x_2 / 24.0));
  }

  d_l1 *= GPS_C;

  return d_l1;
}

/** \} */
