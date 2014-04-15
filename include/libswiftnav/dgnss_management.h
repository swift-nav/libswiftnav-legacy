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
#include "sats_management.h"

#define PHASE_VAR 9e-4 * 9
#define CODE_VAR 100 * 16
#define POS_TRANS_VAR 1e-1
#define VEL_TRANS_VAR 1e-5
#define INT_TRANS_VAR 1e-8
#define POS_INIT_VAR 1e2
#define VEL_INIT_VAR 4e2
#define INT_INIT_VAR 1e4
#define NEW_INT_VAR 1e10

void make_measurements(u8 num_diffs, sdiff_t *sdiffs, double *raw_measurements);
void dgnss_init(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3], double dt);
void dgnss_update(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3], double dt);
void dgnss_rebase_ref(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3], double dt, u8 old_prns[MAX_CHANNELS], sdiff_t *corrected_sdiffs);
kf_t * get_dgnss_kf(void);
s32 * get_stupid_filter_ints(void);
sats_management_t * get_sats_management(void);

s8 dgnss_iar_resolved(void);
void dgnss_float_baseline(u8 *num_used, double b[3]);
void dgnss_fixed_baseline(u8 n, sdiff_t *sdiffs, double ref_ecef[3],
                          u8 *num_used, double b[3]);

