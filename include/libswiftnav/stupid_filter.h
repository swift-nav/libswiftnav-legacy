/*
 * Copyright (C) 2014 Swift Navigation Inc.
 * Contact: Ian Horn <ian@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef LIBSWIFTNAV_STUPID_FILTER_H
#define LIBSWIFTNAV_STUPID_FILTER_H

#include "common.h"
#include "constants.h"
#include "single_diff.h"

typedef struct {
  s32 N[MAX_CHANNELS];
} stupid_filter_state_t;

void amb_from_baseline(u8 num_sats, double *DE, double *dd_meas,
                       double b[3], s32 *N);
void init_stupid_filter(stupid_filter_state_t *s, u8 num_sats, sdiff_t *sdiffs, double *dd_measurements, double b[3], double ref_ecef[3]);
void rebase_stupid_filter(stupid_filter_state_t *s, u8 num_sats, u8 *old_prns, u8 *new_prns);
void update_sats_stupid_filter(stupid_filter_state_t *s, u8 num_old, u8 *old_prns, u8 num_new, sdiff_t *sdiffs, double *dd_measurements, double ref_ecef[3]);
void update_stupid_filter(stupid_filter_state_t *s, u8 num_sats, sdiff_t *sdiffs,
                        double *dd_measurements, double b[3], double ref_ecef[3]);
void lesq_solution(u8 num_dds, double *dd_meas, s32 *N, double *DE, double b[3], double *resid);

#endif /* LIBSWIFTNAV_STUPID_FILTER_H */

