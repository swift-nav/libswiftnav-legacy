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

#include "memory_pool.h"
#include "sats_management.h"

#define MAX_HYPOTHESES 1000

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
  u8 initialized;
  u8 num_matching_ndxs;
  u8 matching_ndxs[MAX_CHANNELS-1];
  s32 ambs[MAX_CHANNELS-1];
} unanimous_amb_check_t; //NOTE maybe do this in a semi-decorrelated space, where more should match sooner.

typedef struct {
  u8 num_dds;
  memory_pool_t *pool;
  residual_mtxs_t res_mtxs;
  sats_management_t sats;
  unanimous_amb_check_t amb_check;
} ambiguity_test_t;

typedef s64 z_t;

/* See doc string above inclusion_loop_body in ambiguity_test.c for info on
 * matrices and lower/upper bounds fields: */
typedef struct {
  u16 intersection_size;
  z_t *Z;
  z_t *Z1;
  z_t *Z1_inv;
  z_t *Z2;
  z_t *Z2_inv;
  z_t *counter;          /* Current position within itr_box bounds. */
  z_t *zimage;           /* Current point in V1 being tested. */
  u8 new_dim;            /* Dimension of new satellite space V2. */
  u8 old_dim;            /* Dimension of old satellite space. */
  z_t *itr_lower_bounds;
  z_t *itr_upper_bounds;
  z_t *box_lower_bounds;
  z_t *box_upper_bounds;
} intersection_count_t; 

typedef struct {
  intersection_count_t *x;
  u8 ndxs_of_old_in_new[MAX_CHANNELS-1];
  u8 ndxs_of_added_in_new[MAX_CHANNELS-1];
  z_t *Z_new_inv;
} generate_hypothesis_state_t2;

void print_s32_mtx_diff(u32 m, u32 n, s32 *Z_inv1, s32 *Z_inv2);
s8 get_single_hypothesis(ambiguity_test_t *amb_test, s32 *hyp_N);
void create_empty_ambiguity_test(ambiguity_test_t *amb_test);
void create_ambiguity_test(ambiguity_test_t *amb_test);
void reset_ambiguity_test(ambiguity_test_t *amb_test);
void destroy_ambiguity_test(ambiguity_test_t *amb_test);
s8 sats_match(const ambiguity_test_t *amb_test, const u8 num_sdiffs, const sdiff_t *sdiffs);
u8 ambiguity_update_reference(ambiguity_test_t *amb_test, const u8 num_sdiffs, const sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first);
void update_ambiguity_test(double ref_ecef[3], double phase_var, double code_var,
                           ambiguity_test_t *amb_test, u8 state_dim, sdiff_t *sdiffs,
                           u8 changed_sats);
void update_unanimous_ambiguities(ambiguity_test_t *amb_test);
u32 ambiguity_test_n_hypotheses(ambiguity_test_t *amb_test);
u8 ambiguity_test_pool_contains(ambiguity_test_t *amb_test, double *ambs);
double ambiguity_test_pool_ll(ambiguity_test_t *amb_test, u8 num_ambs, double *ambs);
double ambiguity_test_pool_prob(ambiguity_test_t *amb_test, u8 num_ambs, double *ambs);
void ambiguity_test_MLE_ambs(ambiguity_test_t *amb_test, s32 *ambs);
void test_ambiguities(ambiguity_test_t *amb_test, double *ambiguity_dd_measurements);
u8 ambiguity_update_sats(ambiguity_test_t *amb_test, const u8 num_sdiffs,
                         const sdiff_t *sdiffs, const sats_management_t *float_sats,
                         const double *float_mean, const double *float_cov_U,
                         const double *float_cov_D);
u8 find_indices_of_intersection_sats(const ambiguity_test_t *amb_test, const u8 num_sdiffs, const sdiff_t *sdiffs_with_ref_first, u8 *intersection_ndxs);
u8 ambiguity_iar_can_solve(ambiguity_test_t *ambiguity_test);
s8 make_dd_measurements_and_sdiffs(u8 ref_prn, u8 *non_ref_prns, u8 num_dds,
                                   u8 num_sdiffs, sdiff_t *sdiffs,
                                   double *ambiguity_dd_measurements, sdiff_t *amb_sdiffs);
s8 make_ambiguity_resolved_dd_measurements_and_sdiffs(ambiguity_test_t *amb_test,
            u8 num_sdiffs, sdiff_t *sdiffs,
            double *ambiguity_dd_measurements, sdiff_t *amb_sdiffs);
s8 make_ambiguity_dd_measurements_and_sdiffs(ambiguity_test_t *amb_test, u8 num_sdiffs, sdiff_t *sdiffs,
                                               double *ambiguity_dd_measurements, sdiff_t *amb_sdiffs);
s8 make_dd_measurements_and_sdiffs(u8 ref_prn, u8 *non_ref_prns, u8 num_dds,
                                   u8 num_sdiffs, sdiff_t *sdiffs,
                                   double *ambiguity_dd_measurements, sdiff_t *amb_sdiffs);
u8 ambiguity_sat_projection(ambiguity_test_t *amb_test, const u8 num_dds_in_intersection, const u8 *dd_intersection_ndxs);
// TODO(dsk) delete
u8 ambiguity_sat_inclusion_old(ambiguity_test_t *amb_test, u8 num_dds_in_intersection,
                               sats_management_t *float_sats, double *float_mean,
                               double *float_cov_U, double *float_cov_D);
u8 ambiguity_sat_inclusion(ambiguity_test_t *amb_test, const u8 num_dds_in_intersection,
                            const sats_management_t *float_sats, const double *float_mean,
                            const double *float_cov_U, const double *float_cov_D);
z_t float_to_decor(const double *addible_float_cov,
                   const double *addible_float_mean,
                   u8 num_addible_dds,
                   u8 num_dds_to_add,
                   z_t *lower_bounds, z_t *upper_bounds,
                   z_t *Z, z_t *Z_inv);
// TODO(dsk) delete
s8 determine_sats_addition(ambiguity_test_t *amb_test,
                           double *float_N_cov, u8 num_float_dds, double *float_N_mean,
                           z_t *lower_bounds, z_t *upper_bounds, u8 *num_dds_to_add,
                           z_t *Z_inv);
// TODO(dsk) delete
void add_sats_old(ambiguity_test_t *amb_test,
                  u8 ref_prn,
                  u32 num_added_dds, u8 *added_prns,
                  z_t *lower_bounds, z_t *upper_bounds,
                  z_t *Z_inv);
void init_residual_matrices(residual_mtxs_t *res_mtxs, u8 num_dds, double *DE_mtx, double *obs_cov);
void assign_residual_covariance_inverse(u8 num_dds, double *obs_cov, double *q, double *r_cov_inv);
void assign_r_vec(residual_mtxs_t *res_mtxs, u8 num_dds, double *dd_measurements, double *r_vec);
void assign_r_mean(residual_mtxs_t *res_mtxs, u8 num_dds, double *hypothesis, double *r_mean);
double get_quadratic_term(residual_mtxs_t *res_mtxs, u8 num_dds, double *hypothesis, double *r_vec);

#endif /* LIBSWIFTNAV_AMBIGUITY_TEST_H */
