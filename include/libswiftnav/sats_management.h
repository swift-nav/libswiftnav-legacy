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

#ifndef LIBSWIFTNAV_SATS_MANAGEMENT_H
#define LIBSWIFTNAV_SATS_MANAGEMENT_H

#include "constants.h"
#include "single_diff.h"

#define OLD_REF 0
#define NEW_REF 1
#define NEW_REF_START_OVER -1

#define INTERSECTION_SATS_THRESHOLD_SIZE 2

typedef struct {
  u8 num_sats;
  u8 prns[MAX_CHANNELS];
} sats_management_t;

void init_sats_management(sats_management_t *sats_management,
                          u8 num_sats, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first);
void print_sats_management(sats_management_t *sats_management);
s8 rebase_sats_management(sats_management_t *sats_management,
                          u8 num_sats, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first);
void update_sats_sats_management(sats_management_t *sats_management, u8 num_non_ref_sdiffs, sdiff_t *non_ref_sdiffs);

void set_reference_sat_of_prns(u8 ref_prn, u8 num_sats, u8 *prns);

#endif /* LIBSWIFTNAV_SATS_MANAGEMENT_H */
