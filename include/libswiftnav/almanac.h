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

#include "common.h"
#include "almanac.h"

/** \addtogroup almanac
 * \{ */

/** Structure containing the GPS almanac for one satellite. */
typedef struct {
  double ecc;   /**< Eccentricity in radians. */
  double toa;   /**< Time of Applicability in seconds since Sunday. */
  double inc;   /**< Inclination in radians. */
  double rora;  /**< Rate of Right Ascension in radians/sec. */
  double a;     /**< Semi-major axis in meters. */
  double raaw;  /**< Right Ascension at Week in radians. */
  double argp;  /**< Argument of Perigee in radians. */
  double ma;    /**< Mean Anomaly at Time of Applicability in radians. */
  double af0;   /**< 0-order clock correction in seconds. */
  double af1;   /**< 1-order clock correction in seconds/second. */
  u16 week;     /**< GPS week number, modulo 1024. */
  u8 prn;       /**< PRN number of the satellite. */
  u8 healthy;   /**< Satellite health status. */
  u8 valid;     /**< Almanac is valid. */
} almanac_t;

/** \} */

void calc_sat_state_almanac(almanac_t* alm, double t, s16 week,
                            double pos[3], double vel[3]);
void calc_sat_az_el_almanac(almanac_t* alm, double t, s16 week,
                            double ref[3], double* az, double* el);
double calc_sat_doppler_almanac(almanac_t* alm, double t, s16 week,
                                double ref[3]);

#endif /* LIBSWIFTNAV_ALMANAC_H */
