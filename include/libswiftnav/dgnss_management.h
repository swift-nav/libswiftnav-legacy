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
#include "amb_kf.h"
#include "sats_management.h"

#define DEFAULT_PHASE_VAR_TEST  (9e-4 * 16)
#define DEFAULT_CODE_VAR_TEST   (100 * 400)
#define DEFAULT_PHASE_VAR_KF    (9e-4 * 16)
#define DEFAULT_CODE_VAR_KF     (100 * 400)
#define DEFAULT_POS_TRANS_VAR   1e-1
#define DEFAULT_VEL_TRANS_VAR   1e-5
#define DEFAULT_INT_TRANS_VAR   1e-8
#define DEFAULT_POS_INIT_VAR    1e2
#define DEFAULT_VEL_INIT_VAR    4e2
#define DEFAULT_AMB_INIT_VAR    1e8
#define DEFAULT_NEW_INT_VAR     1e10

typedef struct {
  double phase_var_test;
  double code_var_test;
  double phase_var_kf;
  double code_var_kf;
  double pos_trans_var;
  double vel_trans_var;
  double int_trans_var;
  double pos_init_var;
  double vel_init_var;
  double amb_init_var;
  double new_int_var;
} dgnss_settings_t;

extern dgnss_settings_t dgnss_settings;

void make_measurements(u8 num_diffs, sdiff_t *sdiffs, double *raw_measurements);
void dgnss_init(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3], double dt);
void dgnss_update(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3], double dt);
void dgnss_rebase_ref(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3], u8 old_prns[MAX_CHANNELS], sdiff_t *corrected_sdiffs);
kf_t * get_dgnss_kf(void);
nkf_t * get_dgnss_nkf(void);
s32 * get_stupid_filter_ints(void);
sats_management_t * get_sats_management(void);

s8 dgnss_iar_resolved(void);
u32 dgnss_iar_num_hyps(void);
void dgnss_reset_iar(void);
void dgnss_init_known_baseline(u8 num_sats, sdiff_t *sdiffs, double receiver_ecef[3], double b[3]);
void dgnss_float_baseline(u8 *num_used, double b[3]);
void dgnss_new_float_baseline(u8 num_sats, sdiff_t *sdiffs, double ref_ecef[3], u8 *num_used, double b[3]);
void dgnss_fixed_baseline(u8 n, sdiff_t *sdiffs, double ref_ecef[3],
                          u8 *num_used, double b[3]);

