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

#ifndef LIBSWIFTNAV_DGNSS_MANAGEMENT_H
#define LIBSWIFTNAV_DGNSS_MANAGEMENT_H

#include "amb_kf.h"
#include "sats_management.h"
#include "ambiguity_test.h"
#include "baseline.h"
#include "constants.h"

typedef struct {
  double phase_var_test;
  double code_var_test;
  double phase_var_kf;
  double code_var_kf;
  double pos_trans_var;
  double vel_trans_var;
  double int_trans_var;
  double amb_drift_var;
  double pos_init_var;
  double vel_init_var;
  double amb_init_var;
  double new_int_var;
} dgnss_settings_t;

typedef struct {
  ambiguities_t fixed_ambs;
  ambiguities_t float_ambs;
} ambiguity_state_t;

typedef struct {
  dgnss_settings_t settings;
  nkf_t nkf;
  sats_management_t sats_management;  /* TODO: This really belongs inside of nkf */
  ambiguity_test_t ambiguity_test;  /* Note: ambiguity_test contains another, different
                                       sats_management (which should be a subset of the above) */
} dgnss_state_t;

void make_measurements(u8 num_diffs, const sdiff_t *sdiffs, double *raw_measurements);
void dgnss_init(dgnss_state_t *dgs, u8 num_sats, const sdiff_t *sdiffs, const double receiver_ecef[3]);
void dgnss_update(dgnss_state_t *dgs, u8 num_sats, const sdiff_t *sdiffs,
                  const double receiver_ecef[3],
                  bool disable_raim, double raim_threshold);
void dgnss_rebase_ref(dgnss_state_t *dgs, u8 num_sats, const sdiff_t *sdiffs, const double receiver_ecef[3],
                      u8 old_prns[MAX_CHANNELS], sdiff_t *corrected_sdiffs);

s8 dgnss_iar_resolved(const dgnss_state_t *dgs);
u32 dgnss_iar_num_hyps(const dgnss_state_t *dgs);
u32 dgnss_iar_num_sats(const dgnss_state_t *dgs);
s8 dgnss_iar_get_single_hyp(const dgnss_state_t *dgs, double *hyp);
void dgnss_reset_iar(dgnss_state_t *dgs); 
void dgnss_init_known_baseline(dgnss_state_t *dgs, u8 num_sats, sdiff_t *sdiffs, double receiver_ecef[3], double b[3]);
void dgnss_update_ambiguity_state(dgnss_state_t *dgs, ambiguity_state_t *s);
s8 dgnss_baseline(u8 num_sdiffs, const sdiff_t *sdiffs,
                  const double ref_ecef[3], const ambiguity_state_t *s,
                  u8 *num_used, double b[3],
                  bool disable_raim, double raim_threshold);
void measure_amb_kf_b(const dgnss_state_t *dgs, u8 num_sdiffs, const sdiff_t *sdiffs,
                      const double receiver_ecef[3], double *b);
void measure_b_with_external_ambs(const dgnss_state_t *dgs, u8 state_dim, const double *state_mean,
                                  u8 num_sdiffs, sdiff_t *sdiffs,
                                  const double receiver_ecef[3], double *b);
void measure_iar_b_with_external_ambs(const dgnss_state_t *dgs, double *state_mean,
                                      u8 num_sdiffs, const sdiff_t *sdiffs,
                                      const double receiver_ecef[3],
                                      double *b);
u8 get_amb_kf_de_and_phase(const dgnss_state_t *dgs, u8 num_sdiffs, sdiff_t *sdiffs,
                           double ref_ecef[3],
                           double *de, double *phase);
u8 get_iar_de_and_phase(const dgnss_state_t *dgs, u8 num_sdiffs, sdiff_t *sdiffs,
                        double ref_ecef[3],
                        double *de, double *phase);
u8 dgnss_iar_pool_contains(const dgnss_state_t *dgs, double *ambs);
double dgnss_iar_pool_ll(const dgnss_state_t *dgs, u8 num_ambs, double *ambs);
double dgnss_iar_pool_prob(const dgnss_state_t *dgs, u8 num_ambs, double *ambs);
u8 get_amb_kf_mean(const dgnss_state_t *dgs, double *ambs);
u8 get_amb_kf_cov(const dgnss_state_t *dgs, double *cov);
u8 get_amb_kf_prns(const dgnss_state_t *dgs, u8 *prns);
u8 get_amb_test_prns(const dgnss_state_t *dgs, u8 *prns);
u8 dgnss_iar_MLE_ambs(const dgnss_state_t *dgs, s32 *ambs);

#endif /* LIBSWIFTNAV_DGNSS_MANAGEMENT_H */
