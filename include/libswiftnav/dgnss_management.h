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


#include "float_kf.h"

#define PHASE_VAR 9e-4
#define CODE_VAR 100
#define POS_TRANS_VAR 1e-1
#define VEL_TRANS_VAR 1e-5
#define INT_TRANS_VAR 1e-8
#define POS_INIT_VAR 1e2
#define VEL_INIT_VAR 4e2
#define INT_INIT_VAR 1e4

void make_measurements(u8 num_diffs, sdiff_t *sdiffs, double *raw_measurements);
void dgnss_init(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3], double dt);
void dgnss_update(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3], double dt);
u8 update_filter(kf_t *kf, u8 num_sats, sdiff_t *sdiffs);
