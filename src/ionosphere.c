/*
 * Copyright (C) 2015, 2016 Swift Navigation Inc.
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
#include <string.h>

#include <libswiftnav/constants.h>
#include <libswiftnav/ionosphere.h>
#include <libswiftnav/logging.h>

/** \defgroup ionosphere Ionospheric models
 * Implemenations of ionospheric delay correction models.
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
 * \return Ionospheric delay distance for GPS L1 frequency [m]
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
  if (fabs(x) >= 1.57) {
    d_l1 = sf * 5e-9;
  } else {
    double x_2 = x * x;
    d_l1 = sf * (5e-9 + amp * (1.0 - x_2 / 2.0 + x_2 * x_2 / 24.0));
  }

  d_l1 *= GPS_C;

  return d_l1;
}

/** The function decodes ionospheric parameters
 * \param subframe4_words pointer to received frame word,
 *        Note: Ionospheric parmeters are passed in subframe 4,
 *        pass function parameters accordingly.
 * \param iono pointer to ionosphere_t where decoded data should be stored
 */
void decode_iono_parameters(const u32 *subframe4_words, ionosphere_t *iono)
{
  union {
    s8 s8;
    u8 u8;
  } onebyte;

  onebyte.u8 = subframe4_words[3-3] >> (30-16) & 0xff; /* alfa 0 */
  iono->a0 = onebyte.s8 * pow(2, -30);
  onebyte.u8 = subframe4_words[3-3] >> (30-24) & 0xff;   /* alfa 1 */
  iono->a1 = onebyte.s8 * pow(2, -27);
  onebyte.u8 = subframe4_words[4-3] >> (30-8) & 0xff;    /* alfa 2 */
  iono->a2 = onebyte.s8 * pow(2, -24);
  onebyte.u8 = subframe4_words[4-3] >> (30-16) & 0xff;    /* alfa 3 */
  iono->a3 = onebyte.s8 * pow(2, -24);
  onebyte.u8 = subframe4_words[4-3] >> (30-24) & 0xff;   /* beta 0 */
  iono->b0 = onebyte.s8 * pow(2, 11);
  onebyte.u8 = subframe4_words[5-3] >> (30-8) & 0xff;    /* beta 1 */
  iono->b1 = onebyte.s8 * pow(2, 14);
  onebyte.u8 = subframe4_words[5-3] >> (30-16) & 0xff;    /* beta 2 */
  iono->b2 = onebyte.s8 * pow(2, 16);
  onebyte.u8 = subframe4_words[5-3] >> (30-24) & 0xff;   /* beta 3 */
  iono->b3 = onebyte.s8 * pow(2, 16);
}

/** \} */
