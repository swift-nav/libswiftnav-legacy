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

/*****************************************************************************
 * Automatically generated from sbp.yaml with generate.py, do not hand edit! *
 *****************************************************************************/

#ifndef LIBSWIFTNAV_SBP_H
#define LIBSWIFTNAV_SBP_H

#include "common.h"


/** GPS Time
 * GPS Time.
 */
#define SBP_GPS_TIME 0x0100

typedef struct __attribute__((packed)) {
  u16 wn; /**< GPS week number */
  u32 tow; /**< GPS Time of Week rounded to the nearest ms */
  s32 ns; /**< Nanosecond remainder of rounded tow */
  u8 flags; /**< Status flags (reserved) */
} sbp_gps_time_t;


/** Dilution of Precision
 * Dilution of Precision.
 */
#define SBP_DOPS 0x0206

typedef struct __attribute__((packed)) {
  u32 tow; /**< GPS Time of Week */
  u16 gdop; /**< Geometric Dilution of Precision */
  u16 pdop; /**< Position Dilution of Precision */
  u16 tdop; /**< Time Dilution of Precision */
  u16 hdop; /**< Horizontal Dilution of Precision */
  u16 vdop; /**< Vertical Dilution of Precision */
} sbp_dops_t;


/** Position in ECEF
 * Position solution in absolute Earth Centered Earth Fixed (ECEF) coordinates.
 */
#define SBP_POS_ECEF 0x0200

typedef struct __attribute__((packed)) {
  u32 tow; /**< GPS Time of Week */
  double x; /**< ECEF X coordinate */
  double y; /**< ECEF Y coordinate */
  double z; /**< ECEF Z coordinate */
  u16 accuracy; /**< Position accuracy estimate */
  u8 n_sats; /**< Number of satellites used in solution */
  u8 flags; /**< Status flags */
} sbp_pos_ecef_t;


/** Geodetic Position
 * Geodetic position solution.
 */
#define SBP_POS_LLH 0x0201

typedef struct __attribute__((packed)) {
  u32 tow; /**< GPS Time of Week */
  double lat; /**< Latitude */
  double lon; /**< Longitude */
  double height; /**< Height */
  u16 h_accuracy; /**< Horizontal position accuracy estimate */
  u16 v_accuracy; /**< Vertical position accuracy estimate */
  u8 n_sats; /**< Number of satellites used in solution */
  u8 flags; /**< Status flags */
} sbp_pos_llh_t;


/** Baseline in ECEF
 * Baseline in Earth Centered Earth Fixed (ECEF) coordinates.
 */
#define SBP_BASELINE_ECEF 0x0202

typedef struct __attribute__((packed)) {
  u32 tow; /**< GPS Time of Week */
  s32 x; /**< Baseline ECEF X coordinate */
  s32 y; /**< Baseline ECEF Y coordinate */
  s32 z; /**< Baseline ECEF Z coordinate */
  u16 accuracy; /**< Position accuracy estimate */
  u8 n_sats; /**< Number of satellites used in solution */
  u8 flags; /**< Status flags (reserved) */
} sbp_baseline_ecef_t;


/** Baseline in NED
 * Baseline in local North East Down (NED) coordinates.
 */
#define SBP_BASELINE_NED 0x0203

typedef struct __attribute__((packed)) {
  u32 tow; /**< GPS Time of Week */
  s32 n; /**< Baseline North coordinate */
  s32 e; /**< Baseline East coordinate */
  s32 d; /**< Baseline Down coordinate */
  u16 h_accuracy; /**< Horizontal position accuracy estimate */
  u16 v_accuracy; /**< Vertical position accuracy estimate */
  u8 n_sats; /**< Number of satellites used in solution */
  u8 flags; /**< Status flags (reserved) */
} sbp_baseline_ned_t;


/** Velocity in ECEF
 * Velocity in Earth Centered Earth Fixed (ECEF) coordinates.
 */
#define SBP_VEL_ECEF 0x0204

typedef struct __attribute__((packed)) {
  u32 tow; /**< GPS Time of Week */
  s32 x; /**< Velocity ECEF X coordinate */
  s32 y; /**< Velocity ECEF Y coordinate */
  s32 z; /**< Velocity ECEF Z coordinate */
  u16 accuracy; /**< Velocity accuracy estimate */
  u8 n_sats; /**< Number of satellites used in solution */
  u8 flags; /**< Status flags (reserved) */
} sbp_vel_ecef_t;


/** Velocity in NED
 * Velocity in local North East Down (NED) coordinates.
 */
#define SBP_VEL_NED 0x0205

typedef struct __attribute__((packed)) {
  u32 tow; /**< GPS Time of Week */
  s32 n; /**< Velocity North coordinate */
  s32 e; /**< Velocity East coordinate */
  s32 d; /**< Velocity Down coordinate */
  u16 h_accuracy; /**< Horizontal velocity accuracy estimate */
  u16 v_accuracy; /**< Vertical velocity accuracy estimate */
  u8 n_sats; /**< Number of satellites used in solution */
  u8 flags; /**< Status flags (reserved) */
} sbp_vel_ned_t;


#endif /* LIBSWIFTNAV_SBP_H */

