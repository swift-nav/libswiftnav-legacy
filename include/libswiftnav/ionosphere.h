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

#ifndef LIBSWIFTNAV_IONOSHPERE_H
#define LIBSWIFTNAV_IONOSHPERE_H

#include "common.h"
#include "gpstime.h"

typedef struct {
  double a0, a1, a2, a3;
  double b0, b1, b2, b3;
} ionosphere_t;

double calc_ionosphere(gps_time_t t_gps,
                       double lat_u, double lon_u,
                       double a, double e,
                       const ionosphere_t *i);

#endif /* LIBSWIFTNAV_IONOSHPERE_H */
