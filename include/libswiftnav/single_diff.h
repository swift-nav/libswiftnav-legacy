/*
 * Copyright (C) 2014 Swift Navigation Inc.
 * Contact: Fergus Noble <fergus@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef LIBSWIFTNAV_SINGLE_DIFF_H
#define LIBSWIFTNAV_SINGLE_DIFF_H

#include "common.h"
#include "track.h"
#include "almanac.h"
#include "ephemeris.h"
#include "gpstime.h"

typedef struct {
  double pseudorange;
  double carrier_phase;
  double doppler;
  double sat_pos[3];
  double sat_vel[3];
  double snr;
  u8 prn;
} sdiff_t;

u8 propagate(u8 n, double ref_ecef[3],
             navigation_measurement_t *m_in_base, gps_time_t *t_base,
             navigation_measurement_t *m_in_rover, gps_time_t *t_rover,
             navigation_measurement_t *m_out_base);

u8 single_diff(u8 n_a, navigation_measurement_t *m_a,
               u8 n_b, navigation_measurement_t *m_b,
               sdiff_t *sds);
void double_diff(u8 n, sdiff_t *sds, sdiff_t *dds, u8 ref_idx);

int sdiff_search_prn(const void *a, const void *b);

u8 make_propagated_sdiffs(u8 n_local, navigation_measurement_t *m_local,
                          u8 n_remote, navigation_measurement_t *m_remote,
                          double *remote_dists, double remote_pos_ecef[3],
                          ephemeris_t *es, gps_time_t t,
                          sdiff_t *sds);

void almanacs_to_single_diffs(u8 n, almanac_t *alms, gps_time_t timestamp, sdiff_t *sdiffs);

s8 copy_sdiffs_put_ref_first(u8 ref_prn, u8 num_sdiffs, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first);

u8 filter_sdiffs(u8 num_sdiffs, sdiff_t *sdiffs, u8 num_sats_to_drop, u8 *sats_to_drop);

#endif /* LIBSWIFTNAV_SINGLE_DIFF_H */

