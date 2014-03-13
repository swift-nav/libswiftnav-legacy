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

#ifndef LIBSWIFTNAV_AMBIGUITY_TEST_H
#define LIBSWIFTNAV_AMBIGUITY_TEST_H

#include "constants.h"
#include "common.h"
#include "memory_pool.h"
#include "sats_management.h"

#define MAX_HYPOTHESES 20000

typedef struct {
  s32 N[MAX_CHANNELS-1];
  float ll;
} hypothesis_t;

typedef struct {
  u32 res_dim;
  u8 null_space_dim;
  double null_projector[(MAX_CHANNELS-4) * (MAX_CHANNELS-1)];
  double half_res_cov_inv[(2*MAX_CHANNELS - 5) * (2*MAX_CHANNELS - 5)];
} residual_mtxs_t;

typedef struct {
  u8 num_dds;
  memory_pool_t *pool;
  residual_mtxs_t res_mtxs;
  sats_management_t sats;
} ambiguity_test_t;

void create_ambiguity_test(ambiguity_test_t *amb_test);
void destroy_ambiguity_test(ambiguity_test_t *amb_test);
void init_ambiguity_test(ambiguity_test_t *amb_test, u32 state_dim, u8 *prns, sdiff_t *sdiffs, 
                         double *float_mean, double *float_cov, double *DE_mtx, double *obs_cov);
void print_hyp(void *arg, element_t *elem);
s8 sats_match(ambiguity_test_t *amb_test, u8 num_sdiffs, sdiff_t *sdiffs);
void ambiguity_update_reference(ambiguity_test_t *amb_test, u8 num_sdiffs, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first);
void update_ambiguity_test(ambiguity_test_t *amb_test, u32 state_dim, sats_management_t *float_sats, sdiff_t *sdiffs,
                           double *float_mean, double *float_cov, double *DE_mtx, double *obs_cov);
void add_sats(ambiguity_test_t *amb_test, 
              u32 num_added_dds, u8 *added_dd_prns,
              double *float_mean, double *float_cov_diag,
              s32 *Z_inv, double *new_DE_mtx, double *new_obs_cov);
void init_residual_matrices(residual_mtxs_t *res_mtxs, u8 num_dds, double *DE_mtx, double *obs_cov);
void assign_phase_obs_null_basis(u8 num_dds, double *DE_mtx, double *q);
void assign_residual_covariance_inverse(u8 num_dds, double *obs_cov, double *q, double *r_cov_inv);
void assign_r_vec(residual_mtxs_t *res_mtxs, u8 num_dds, double *dd_measurements, double *r_vec);
void assign_r_mean(residual_mtxs_t *res_mtxs, u8 num_dds, double *hypothesis, double *r_mean);
double get_quadratic_term(residual_mtxs_t *res_mtxs, u8 num_dds, double *hypothesis, double *r_vec);

#endif /* LIBSWIFTNAV_AMBIGUITY_TEST_H */