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

#ifndef LIBSWIFTNAV_RTCM3_H
#define LIBSWIFTNAV_RTCM3_H

#include "common.h"

#include "gpstime.h"
#include "track.h"

s16 rtcm3_check_frame(u8 *buff);
s8 rtcm3_write_frame(u16 len, u8 *buff);

void rtcm3_write_header(u8 *buff, u16 type, u16 id, gps_time_t t,
                        u8 sync, u8 n_sat, u8 div_free, u8 smooth);
void rtcm3_read_header(u8 *buff, u16 *type, u16 *id, double *tow,
                       u8 *sync, u8 *n_sat, u8 *div_free, u8 *smooth);

u16 rtcm3_encode_1002(u8 *buff, u16 id, gps_time_t t, u8 n_sat,
                      navigation_measurement_t *nm, u8 sync);
s8 rtcm3_decode_1002(u8 *buff, u16 *id, double *tow, u8 *n_sat,
                     navigation_measurement_t *nm, u8 *sync);

#endif /* LIBSWIFTNAV_RTCM3_H */
