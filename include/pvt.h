/*
 * Copyright (C) 2010 Swift Navigation Inc.
 * Contact: Henry Hallam <henry@swift-nav.com>
 *          Matt Peddie <peddie@alum.mit.edu>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef LIBSWIFTNAV_PVT_H
#define LIBSWIFTNAV_PVT_H

#include "common.h"
#include "track.h"

/* GPS constants */
/** \todo These constants don't belong in this file. */
#define GPS_NUM_SATS  32
#define GPS_PI      3.14159265358979323846
#define GPS_L1_HZ     1.57542E9
#define NAV_OMEGAE_DOT 7.2921151467e-005
#define NAV_GM 3.986005e14
#define NAV_C 299792458.0
#define NAV_F -4.442807633e-010
#define NOMINAL_RANGE 22980E3 /* Average distance to satellite, meters. */
#define R2D (180.0/GPS_PI)    /* Radians to degrees. */
#define D2R (GPS_PI/180.0)    /* Degrees to radians. */

#define GNSS_MAX_OBS GPS_NUM_SATS
#define GNSS_MAX_RECEIVERS 1

#define PVT_MAX_ITERATIONS 20

typedef struct {
  double pdop;
  double gdop;
  double tdop;
  double hdop;
  double vdop;
} dops_t;

typedef struct __attribute__((packed)) {
  /*
   * Be careful of stuct packing to avoid (very mild) slowness,
   * try to keep all the types aligned i.e. put the 64bit
   * things together at the top, then the 32bit ones etc.
   */
  /** Receiver position latitude [deg], longitude [deg], altitude [m] */
  double pos_llh[3];
  /** Receiver position ECEF XYZ [m] */
  double pos_ecef[3];
  /** Receiver velocity in NED [m/s] */
  double vel_ned[3];
  /** Receiver velocity in ECEF XYZ [m/s] */
  double vel_ecef[3];

  /* This is the row-first upper diagonal matrix of error covariances
   * in x, y, z (all receiver clock covariance terms are ignored).  So
   * it goes like so:
   *
   *    0  1  2
   *    _  3  4
   *    _  _  5
   *
   *    Index 6 is the GDOP.
   */
  double err_cov[7];

  double clock_offset;
  double clock_bias;

  /* GPS time */
  gps_time_t time;

  /* 0 = invalid, 1 = code phase */
  u8 valid;
  /* Number of channels used in the soluton. */
  u8 n_used;
} gnss_solution;

u8 calc_PVT(const u8 n_used,
            const navigation_measurement_t const nav_meas[n_used],
            gnss_solution *soln,
            dops_t *dops);

#endif /* LIBSWIFTNAV_PVT_H */

