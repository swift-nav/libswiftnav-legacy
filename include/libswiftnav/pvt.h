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

#include <libswiftnav/common.h>
#include <libswiftnav/track.h>

#define PVT_MAX_ITERATIONS 10

typedef struct {
  double pdop;
  double gdop;
  double tdop;
  double hdop;
  double vdop;
} dops_t;

extern const char *pvt_err_msg[7];

#define PVT_CONVERGED_NO_RAIM        2
#define PVT_CONVERGED_RAIM_REPAIR    1
#define PVT_CONVERGED_RAIM_OK        0
#define PVT_PDOP_TOO_HIGH           -1
#define PVT_BAD_ALTITUDE            -2
#define PVT_VELOCITY_LOCKOUT        -3
#define PVT_RAIM_REPAIR_FAILED      -4
#define PVT_RAIM_REPAIR_IMPOSSIBLE  -5
#define PVT_UNCONVERGED             -6
#define PVT_INSUFFICENT_MEAS        -7

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

s8 calc_PVT(const u8 n_used,
            const navigation_measurement_t nav_meas[n_used],
            bool disable_raim,
            gnss_solution *soln,
            dops_t *dops);

#endif /* LIBSWIFTNAV_PVT_H */
