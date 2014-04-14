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

#ifndef LIBSWIFTNAV_CONSTANTS_H
#define LIBSWIFTNAV_CONSTANTS_H

#include <math.h>

/** \defgroup constants Constants
 * Useful constants.
 * \{ */

#define MAX_CHANNELS 11 /**< Maximum sats we can track */
#define MAX_SATS 32 /**< Maximum sats in the universe */

#define R2D (180.0 / M_PI) /**< Conversion factor from radians to degrees. */
#define D2R (M_PI / 180.0) /**< Conversion factor from degrees to radians. */

/** \defgroup gps_constants GPS
 * Constants related to the Global Positioning System.
 * See ICD-GPS-200C.
 * \{ */

/** The official GPS value of Pi.
 * This is the value used by the CS to curve fit ephemeris parameters and
 * should be used in all ephemeris calculations. */
#define GPS_PI 3.14159265358979323846

/** The GPS L1 center frequency in Hz. */
#define GPS_L1_HZ 1.57542e9

/** Earth's rotation rate as defined in the ICD in rad / s
 * \note This is actually not identical to the usual WGS84 definition. */
#define GPS_OMEGAE_DOT 7.2921151467e-5

/** Earth’s Gravitational Constant as defined in the ICD in m^3 / s^2
 * \note This is actually not identical to the usual WGS84 definition. */
#define GPS_GM 3.986005e14

/** The official GPS value of the speed of light in m / s. 
 * \note This is the exact value of the speed of light in vacuum (by the definition of meters). */
#define GPS_C 299792458.0

/** The speed of light in air at standard temperature and pressure.
 * \note This is GPS_C / mu where mu is 1.0002926 */
#define GPS_C_NO_VAC (GPS_C / 1.0002926)

/** The wavelength of L1 in a vacuum.
 * \note This is GPS_C / GPS_L1_HZ. */
#define GPS_L1_LAMBDA (GPS_C / GPS_L1_HZ)

/** The wavelength of L1 in air at standard temperature and pressure.
 * \note This is GPS_C_NO_VAC / GPS_L1_HZ. */
#define GPS_L1_LAMBDA_NO_VAC (GPS_C_NO_VAC / GPS_L1_HZ)

/** Approximate average distance to the GPS satellites in m. */
#define GPS_NOMINAL_RANGE 22.980e6

/** GPS C/A code chipping rate in Hz. */
#define GPS_CA_CHIPPING_RATE 1.023e6

/* \} */

/* \} */

#endif /* LIBSWIFTNAV_CONSTANTS_H */


