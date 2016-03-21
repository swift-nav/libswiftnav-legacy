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

#define R2D (180.0 / M_PI) /**< Conversion factor from radians to degrees. */
#define D2R (M_PI / 180.0) /**< Conversion factor from degrees to radians. */

/** \defgroup gps_constants GPS
 * Constants related to the Global Positioning System.
 * See ICD-GPS-200C.
 * \{ */

/** The official GPS value of Pi.
 * This is the value used by the CS to curve fit ephemeris parameters and
 * should be used in all ephemeris calculations. */
#define GPS_PI 3.1415926535898

/** The GPS L1 center frequency in Hz. */
#define GPS_L1_HZ 1.57542e9

/** The GPS L2 center frequency in Hz. */
#define GPS_L2_HZ 1.2276e9

/** Earth's rotation rate as defined in the ICD in rad / s
 * \note This is actually not identical to the usual WGS84 definition. */
#define GPS_OMEGAE_DOT 7.2921151467e-5

/** Earth’s Gravitational Constant as defined in the ICD in m^3 / s^2
 * \note This is actually not identical to the usual WGS84 definition. */
#define GPS_GM 3.986005e14

/** The official GPS value of the speed of light in m / s.
 * \note This is the exact value of the speed of light in vacuum (by the definition of meters). */
#define GPS_C 299792458.0

/** The official GPS value of the relativistic clock correction coefficient F. */
#define GPS_F -4.442807633e-10

/** The speed of light in air at standard temperature and pressure.
 * \note This is GPS_C / mu where mu is 1.0002926 */
#define GPS_C_NO_VAC (GPS_C / 1.0002926)

/** The wavelength of L1 in air at standard temperature and pressure.
 * \note This is GPS_C_NO_VAC / GPS_L1_HZ. */
#define GPS_L1_LAMBDA_NO_VAC (GPS_C_NO_VAC / GPS_L1_HZ)

/** The wavelength of L2 in air at standard temperature and pressure.
 * \note This is GPS_C_NO_VAC / GPS_L2_HZ. */
#define GPS_L2_LAMBDA_NO_VAC (GPS_C_NO_VAC / GPS_L2_HZ)

/** Approximate average distance to the GPS satellites in m. */
#define GPS_NOMINAL_RANGE 22.980e6

/** GPS C/A code chipping rate in Hz. */
#define GPS_CA_CHIPPING_RATE 1.023e6

/* \} */

/** \defgroup dgnss_constants DGNSS
 * Approximate variance values used by the KF and IAR hypothesis test.
 * \{ */

/** The default DD carrier phase variance to use in the hypothesis testing. */
#define DEFAULT_PHASE_VAR_TEST  (9e-4 * 16)
/** The default DD pseudorange variance to use in the hypothesis testing. */
#define DEFAULT_CODE_VAR_TEST   (100 * 400)
/** The default DD carrier phase variance to use in the Kalman filter. */
#define DEFAULT_PHASE_VAR_KF    (9e-4 * 16)
/** The default DD pseudorange variance to use in the Kalman filter. */
#define DEFAULT_CODE_VAR_KF     (100 * 400)
/** The default variance of the process noise Kalman filter. Its particular use
 * is different from that of a normal KF process noise. It's still a random
 * walk, but in a special space. Look at the code for its current usage.*/
#define DEFAULT_AMB_DRIFT_VAR   1e-8
/** The variance with which to initialize the Kalman Filter. */
#define DEFAULT_AMB_INIT_VAR    1e25
/** The variance with which to add new sats to the Kalman Filter.
 * TODO deprecate in lieu of amb_init_var once we do some tuning. */
#define DEFAULT_NEW_INT_VAR     1e25

/* \} */

/* \} */

#endif /* LIBSWIFTNAV_CONSTANTS_H */
