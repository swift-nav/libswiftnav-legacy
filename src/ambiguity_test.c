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

#include <clapack.h>
#include <cblas.h>
#include <stdio.h>
#include <string.h>
#include "ambiguity_test.h"
#include "common.h"
#include "constants.h"
#include "linear_algebra.h"
#include "single_diff.h"
#include "float_kf.h"
#include "lambda.h"

#define RAW_PHASE_BIAS_VAR 0
#define DECORRELATED_PHASE_BIAS_VAR 0
#define NUM_SEARCH_STDS 8

void create_ambiguity_test(ambiguity_test_t *amb_test)
{
  amb_test->pool = memory_pool_new(MAX_HYPOTHESES, sizeof(hypothesis_t));
}

void destroy_ambiguity_test(ambiguity_test_t *amb_test)
{
  memory_pool_destroy(amb_test->pool);
}

void init_ambiguity_test(ambiguity_test_t *amb_test, u32 state_dim, u8 *float_prns, sdiff_t *sdiffs, double *float_mean,
                         double *float_cov, double *DE_mtx, double *obs_cov)
{
  u8 num_dds = state_dim-6;
  double decor_float_cov_diag[num_dds];
  double float_cov_N[num_dds * num_dds];
  double Z[num_dds * num_dds];
  double Z_inv_[num_dds * num_dds];

  /* Re-shape float_cov to contain just ambiguity states. */
  for (u8 i=0; i<num_dds; i++) {
    for (u8 j=0; j<num_dds; j++) {
      float_cov_N[i*num_dds + j] = float_cov[(i+6)*state_dim + (j+6)];
      #if RAW_PHASE_BIAS_VAR != 0
      if (i == j) {
        float_cov_N[i*num_dds + j] += RAW_PHASE_BIAS_VAR;
      }
      #endif
    }
  }

  lambda_reduction(num_dds, float_cov_N, Z);

  memset(decor_float_cov_diag, 0, num_dds * sizeof(double));
  for (u8 i=0; i < num_dds; i++) {
    for (u8 j=0; j < num_dds; j++) {
      for (u8 k=0; k < num_dds; k++) {
        decor_float_cov_diag[i] += Z[i*num_dds + j] * float_cov_N[j * num_dds + k] * Z[i*num_dds + k];
      }
    }
    #if DECORRELATED_PHASE_BIAS_VAR != 0
    decor_float_cov_diag[i] += DECORRELATED_PHASE_BIAS_VAR;
    #endif
  }

  matrix_inverse(num_dds, Z, Z_inv_);

  MAT_PRINTF(Z, num_dds, num_dds);
  MAT_PRINTF(Z_inv_, num_dds, num_dds);

  s32 Z_inv[num_dds * num_dds];
  for (u8 i=0; i < num_dds; i++) {
    for (u8 j=0; j < num_dds; j++) {
      Z_inv[i*num_dds + j] = (s32) lround(Z_inv_[i*num_dds + j]);
    }
  }

  double decor_float_mean[num_dds];
  memset(decor_float_mean, 0, num_dds * sizeof(double));
  for (u8 i=0; i < num_dds; i++) {
    for (u8 j=0; j < num_dds; j++) {
      decor_float_mean[i] += Z[i*num_dds + j] * float_mean[6+j];
    }
    printf("decor_float_mean[%u] = %f\n", i, decor_float_mean[i]);
  }



  /* Initialize pool with single element with num_dds = 0, i.e.
   * zero length N vector, i.e. no satellites. When we take the
   * product of this single element with the set of new satellites
   * we will just get a set of elements corresponding to the new sats. */
  hypothesis_t *empty_element = (hypothesis_t *)memory_pool_add(amb_test->pool);
  /* Start with ll = 0, just for the sake of argument. */
  empty_element->ll = 0;
  amb_test->sats.num_sats = 0;
  add_sats(amb_test, num_dds, &float_prns[1], decor_float_mean, decor_float_cov_diag, Z_inv);
  /* Update the rest of the amb_test state with the new sats. */
  init_residual_matrices(&amb_test->res_mtxs, num_dds + 1, DE_mtx, obs_cov);
  init_sats_management(&amb_test->sats, num_dds + 1, sdiffs, 0);
}

// void add_sats(ambiguity_test_t *amb_test, 
//               u32 num_added_dds, u8 *added_prns,
//               double *float_mean, double *float_cov_diag,
//               s32 *Z_inv, double *new_DE_mtx, double *new_obs_cov)


void update_ambiguity_test(ambiguity_test_t *amb_test, u32 state_dim, sats_management_t *float_sats, sdiff_t *sdiffs,
                           double *float_mean, double *float_cov)
{
  (void) amb_test;
  (void) state_dim;
  (void) sdiffs;
  (void) float_mean;
  (void) float_cov;
  (void) float_sats;
  u8 num_sdiffs = state_dim-5;
  u8 changed_sats = ambiguity_update_sats(amb_test, num_sdiffs, sdiffs,
                                          float_mean, float_cov);
  if (1 == 1 || changed_sats == 1) {
    double DE_mtx[(amb_test->sats.num_sats-1) * 3];
    double obs_cov[(amb_test->sats.num_sats-1) * (amb_test->sats.num_sats-1)];
    init_residual_matrices(&amb_test->res_mtxs, amb_test->sats.num_sats, DE_mtx, obs_cov);
  }
  
  double ambiguity_dd_measurements[amb_test->sats.num_sats-1];
  make_ambiguity_dd_measurements(amb_test, num_sdiffs, sdiffs, ambiguity_dd_measurements);

  
  // test_ambiguities(amb_test, ambiguity_dd_measurements);
}

// void test_ambiguities(ambiguity_test_t *amb_test, double *ambiguity_dd_measurements) {
//   double 
// }

void make_ambiguity_dd_measurements(ambiguity_test_t *amb_test, u8 num_sdiffs, sdiff_t *sdiffs, double *ambiguity_dd_measurements)
{
  u8 ref_prn = amb_test->sats.prns[0];
  double ref_phase;
  double ref_pseudorange;
  u8 i=0;
  u8 j=0;
  u8 *amb_test_non_dd_prns = &amb_test->sats.prns[1];
  u8 num_dds = amb_test->sats.num_sats-1;
  u8 found_ref = 0; //DEBUG
  while (i < num_dds) {
    if (amb_test_non_dd_prns[i] == sdiffs[i].prn) {
      ambiguity_dd_measurements[i] = sdiffs[i].carrier_phase;
      ambiguity_dd_measurements[i+num_dds] = sdiffs[i].pseudorange;
      i++;
      j++;
    } else if (ref_prn == sdiffs[i].prn) {
      ref_phase =  sdiffs[i].carrier_phase;
      ref_pseudorange = sdiffs[i].pseudorange;
      found_ref = 1; //DEBUG
    }
    // else {
    //   j++;
    // }
    else if (amb_test_non_dd_prns[i] > sdiffs[i].prn) {
      j++;
    } else { //DEBUG
      //then the ref_prn wasn't in the sdiffs and something has gone wrong in setting up/rebasing amb_test->sats
      printf("there is either disorder in amb_test->sats or it contains a sat not in sdiffs. amb_test->sats must be a subset of sdiffs by this point.\n");
      printf("amb_test->sats.prns = [");
      for (u8 j=0; j < amb_test->sats.num_sats; j++) {
        printf("%d, ", amb_test->sats.prns[i]);
      }
      printf("]\n");
      printf("sdiffs.prns = [");
      for (u8 j=0; j < num_sdiffs; j++) {
        printf("%d, ", sdiffs[i].prn);
      }
      printf("]\n");
    }
  }
  if (found_ref == 0) { //DEBUG
    //then the ref_prn wasn't in the sdiffs and something has gone wrong in setting up/rebasing amb_test->sats
    printf("amb_test->sats's reference wasn't foung in the sdiffs, but it should have already been rebased.\n");
    printf("amb_test->sats.prns = [");
    for (u8 j=0; j < amb_test->sats.num_sats; j++) {
      printf("%d, ", amb_test->sats.prns[i]);
    }
    printf("]\n");
    printf("sdiffs.prns = [");
    for (u8 j=0; j < num_sdiffs; j++) {
      printf("%d, ", sdiffs[i].prn);
    }
    printf("]\n");
  }
  for (u8 i=0; i < num_dds; i++) {
    ambiguity_dd_measurements[i] -= ref_phase;
    ambiguity_dd_measurements[i+num_dds] -= ref_pseudorange;
  }
}

s8 sats_match(ambiguity_test_t *amb_test, u8 num_sdiffs, sdiff_t *sdiffs)
{
  if (amb_test->sats.num_sats != num_sdiffs) {
    return 0;
  }
  u8 *prns = amb_test->sats.prns;
  u8 amb_ref = amb_test->sats.prns[0];
  u8 j=0;
  for (u8 i = 1; i<amb_test->sats.num_sats-1; i++) {
    if (prns[i] == sdiffs[j].prn) {
      j++;
    }
    else if (amb_ref == sdiffs[j].prn) {
      j++;
      i--;
    }
    else {
      return 0;
    } 
  }
  return 1;
}

typedef struct {
  u8 num_sats;
  u8 old_prns[MAX_CHANNELS];
  u8 new_prns[MAX_CHANNELS];
} rebase_prns_t;

void rebase_hypothesis(void *arg, element_t *elem) //TODO make it so it doesn't have to do all these lookups every time
{
  rebase_prns_t *prns = (rebase_prns_t *) arg;
  u8 num_sats = prns->num_sats;
  u8 *old_prns = prns->old_prns;
  u8 *new_prns = prns->new_prns;

  hypothesis_t *hypothesis = (hypothesis_t *)elem;

  u8 old_ref = old_prns[0];
  u8 new_ref = new_prns[0];

  s32 new_N[num_sats-1];
  s32 index_of_new_ref_in_old = find_index_of_element_in_u8s(num_sats, new_ref, &old_prns[1]);
  s32 val_for_new_ref_in_old_basis = hypothesis->N[index_of_new_ref_in_old];
  for (u8 i=0; i<num_sats-1; i++) {
    u8 new_prn = new_prns[1+i];
    if (new_prn == old_ref) {
      new_N[i] = - val_for_new_ref_in_old_basis;
    }
    else {
      s32 index_of_this_sat_in_old_basis = find_index_of_element_in_u8s(num_sats, new_prn, &old_prns[1]);
      new_N[i] = hypothesis->N[index_of_this_sat_in_old_basis] - val_for_new_ref_in_old_basis;
    }
  }
  memcpy(hypothesis->N, new_N, (num_sats-1) * sizeof(s32));
}

void print_sats_man(sats_management_t *sats_man) {
  for (u8 i=0; i<sats_man->num_sats; i++) {
    printf("%d,", sats_man->prns[i]);
  }
  printf("\n");
}

u8 ambiguity_update_reference(ambiguity_test_t *amb_test, u8 num_sdiffs, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  u8 changed_ref = 0;
  u8 old_prns[amb_test->sats.num_sats];
  memcpy(old_prns, amb_test->sats.prns, amb_test->sats.num_sats * sizeof(u8));

  print_sats_man(&amb_test->sats);
  s8 sats_management_code = rebase_sats_management(&amb_test->sats, num_sdiffs, sdiffs, sdiffs_with_ref_first);
  print_sats_man(&amb_test->sats);
  if (sats_management_code != OLD_REF) {
    changed_ref = 1;
    u8 new_prns[amb_test->sats.num_sats];
    memcpy(new_prns, amb_test->sats.prns, amb_test->sats.num_sats * sizeof(u8));

    rebase_prns_t prns = {.num_sats = amb_test->sats.num_sats};
    memcpy(prns.old_prns, old_prns, amb_test->sats.num_sats * sizeof(u8));
    memcpy(prns.new_prns, new_prns, amb_test->sats.num_sats * sizeof(u8));
    memory_pool_map(amb_test->pool, &prns, &rebase_hypothesis);
  }
  return changed_ref;
}

typedef struct {
  u8 num_ndxs;
  u8 intersection_ndxs[MAX_CHANNELS-1];
} intersection_ndxs_t;

s32 projection_comparator(void *arg, element_t *a, element_t *b)
{
  intersection_ndxs_t *intersection_struct = (intersection_ndxs_t *) arg;
  u8 num_ndxs = intersection_struct->num_ndxs;
  u8 *intersection_ndxs = intersection_struct->intersection_ndxs;

  hypothesis_t *hyp_a = (hypothesis_t *) a;
  hypothesis_t *hyp_b = (hypothesis_t *) b;

  for (u8 i=0; i<num_ndxs; i++) {
    u8 ndxi = intersection_ndxs[i];
    s32 ai = hyp_a->N[ndxi];
    s32 bi = hyp_b->N[ndxi];
    if (ai == bi) {
      continue;
    }
    if (ai < bi) {
      return -1;
    }
    if (ai > bi) {
      return 1;
    }
  }
  return 0;
}

void projection_aggregator(element_t *new_, void *x_, u32 n, element_t *elem_)
{
  intersection_ndxs_t *x = (intersection_ndxs_t *)x_;
  hypothesis_t *new = (hypothesis_t *)new_;
  hypothesis_t *elem = (hypothesis_t *)elem_;

  if (n == 0) {
    for (u8 i=0; i<x->num_ndxs; i++) {
      u8 ndxi = x->intersection_ndxs[i];
      new->N[i] = elem->N[ndxi];
    }
    new->ll = elem->ll;
  }
  else {
    new->ll += log(1 + exp(elem->ll - new->ll));
    // new->ll = MAX(new->ll, elem->ll);
  }
  
}

u8 ambiguity_sat_projection(ambiguity_test_t *amb_test, u8 num_dds_in_intersection, u8 *dd_intersection_ndxs)
{
  if (amb_test->sats.num_sats == num_dds_in_intersection) {
    return 0;
  }

  intersection_ndxs_t intersection = {.num_ndxs = num_dds_in_intersection};
  memcpy(intersection.intersection_ndxs, dd_intersection_ndxs, num_dds_in_intersection * sizeof(u8));
  memory_pool_group_by(amb_test->pool, 
                       &intersection, &projection_comparator, 
                       &intersection, sizeof(intersection), 
                       &projection_aggregator);
  u8 work_prns[MAX_CHANNELS];
  memcpy(work_prns, amb_test->sats.prns, amb_test->sats.num_sats * sizeof(u8));
  for (u8 i=0; i<num_dds_in_intersection; i++) {
    amb_test->sats.prns[i] = work_prns[dd_intersection_ndxs[i]];
  }
  amb_test->sats.num_sats = num_dds_in_intersection;
  return 1;
}


u8 ambiguity_sat_inclusion(ambiguity_test_t *amb_test, u8 num_dds_in_intersection, u8 *dd_intersection_ndxs,
                             sats_management_t *float_sats, double *float_mean, double *float_cov)
{
  (void) dd_intersection_ndxs;
  if (amb_test->sats.num_sats == num_dds_in_intersection) {
    return 0;
  }
  double N_cov[(float_sats->num_sats-1) * (float_sats->num_sats-1)];
  u32 state_dim = float_sats->num_sats + 5;
  for (u8 i = 0; i < float_sats->num_sats-1; i++) {
    memcpy(&N_cov[i*(float_sats->num_sats-1)],
           &float_cov[(6+i)*state_dim + 6],
           (float_sats->num_sats - 1) * sizeof(double));
  }
  u8 float_prns[float_sats->num_sats];
  memcpy(float_prns, float_sats->prns, float_sats->num_sats * sizeof(u8));
  double N_mean[float_sats->num_sats-1];
  if (amb_test->sats.prns[0] != float_sats->prns[0]) {
    u8 old_prns[float_sats->num_sats];
    memcpy(old_prns, float_sats->prns, float_sats->num_sats * sizeof(u8));
    memcpy(N_mean, &float_mean[6], (float_sats->num_sats-1) * sizeof(double));
    set_reference_sat_of_prns(amb_test->sats.prns[0], float_sats->num_sats, float_prns);
    rebase_mean(N_mean, float_sats->num_sats, old_prns, float_prns);

    double rebase_mtx[(float_sats->num_sats-1) * (float_sats->num_sats-1)];
    assign_state_rebase_mtx(float_sats->num_sats, old_prns, float_prns, rebase_mtx);
    double intermediate_cov[(float_sats->num_sats-1) * (float_sats->num_sats-1)];
    //TODO make more efficient via structure of rebase_mtx
    cblas_dsymm(CblasRowMajor, CblasRight, CblasUpper, //CBLAS_ORDER, CBLAS_SIDE, CBLAS_UPLO
              float_sats->num_sats-1, float_sats->num_sats-1, // int M, int N
              1, N_cov, float_sats->num_sats-1, // double alpha, double *A, int lda
              rebase_mtx, float_sats->num_sats-1, // double *B, int ldb
              0, intermediate_cov, float_sats->num_sats-1); // double beta, double *C, int ldc
    // MAT_PRINTF(intermediate_cov, state_dim, state_dim);

    //TODO make more efficient via the structure of rebase_mtx
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, // CBLAS_ORDER, CBLAS_TRANSPOSE transA, cBLAS_TRANSPOSE transB
              float_sats->num_sats-1, float_sats->num_sats-1, float_sats->num_sats-1, // int M, int N, int K
              1, intermediate_cov, float_sats->num_sats-1, // double alpha, double *A, int lda
              rebase_mtx, float_sats->num_sats-1, //double *B, int ldb
              0, N_cov, float_sats->num_sats-1); //beta, double *C, int ldc
  }
  //now float_prns has the correct reference, as does N_cov

  //next we add the new sats
  //loop through the new sat sets in decreasing number (for now, in prn order), adding when it suits us

  /* first set up the largest prn lists of added sats
   * then get those prns covariances
   * then loop through:
   *    test if we'll add these
   *    if so:
   *      add them
   *      break loop
   *    else:
   *      decrease size of cov set
   */

  //first get all the prns we might add, and their covariances
  u8 i = 1;
  u8 j = 1;
  u8 num_addible_dds = 0;
  u8 ndxs_of_new_dds_in_float[MAX_CHANNELS-1];
  u8 new_dd_prns[MAX_CHANNELS-1];
  while (j < float_sats->num_sats) {
    if (i < amb_test->sats.num_sats && amb_test->sats.prns[i] == float_prns[j]) {
      i++;
      j++;
    } else { //else float_sats[j] is a new one
      ndxs_of_new_dds_in_float[num_addible_dds] = j;
      new_dd_prns[num_addible_dds] = float_prns[j];
      num_addible_dds++;
    }
  }

  double addible_float_cov[num_addible_dds * num_addible_dds];
  for (i=0; i < num_addible_dds; i++) {
    for (j=0; j < num_addible_dds; j++) {
      addible_float_cov[i*num_addible_dds + j] = N_cov[ndxs_of_new_dds_in_float[i]*(float_sats->num_sats-1) + ndxs_of_new_dds_in_float[j]];
    }
  }

  //now loop through, adding once we can
  u8 num_dds_to_add = num_addible_dds;
  double added_float_cov[num_dds_to_add * num_dds_to_add];
  double decor_float_cov_diag[num_dds_to_add];
  double Z[num_dds_to_add * num_dds_to_add];
  double rough_multiplier;
  double decor_float_mean[num_dds_to_add];
  double Z_inv_[num_dds_to_add * num_dds_to_add];
  s32 Z_inv[num_dds_to_add * num_dds_to_add];
  while (num_dds_to_add > 0) {
    //TODO just allow stride in lambda_reduction, instead of creating intermediate matrices

    for (i=0; i<num_dds_to_add; i++) {
      for (j=0; j<num_dds_to_add; j++) {
        added_float_cov[i*num_dds_to_add + j] = addible_float_cov[i*num_addible_dds + j];
        #if RAW_PHASE_BIAS_VAR != 0
        if (i == j) {
          added_float_cov[i*num_dds_to_add + j] += RAW_PHASE_BIAS_VAR;
        }
        #endif
      }
    }

    lambda_reduction(num_dds_to_add, added_float_cov, Z);

    memset(decor_float_cov_diag, 0, num_dds_to_add * sizeof(double));
    rough_multiplier = 1;
    for (i=0; i < num_dds_to_add; i++) {
      for (j=0; j < num_dds_to_add; j++) {
        for (u8 k=0; k < num_dds_to_add; k++) {
          decor_float_cov_diag[i] += Z[i*num_dds_to_add + j] * added_float_cov[j * num_dds_to_add + k] * Z[i*num_dds_to_add + k];
        }
      }
      #if DECORRELATED_PHASE_BIAS_VAR != 0
      decor_float_cov_diag[i] += DECORRELATED_PHASE_BIAS_VAR;
      #endif
      rough_multiplier *= MIN(2, 2 * NUM_SEARCH_STDS * sqrt(decor_float_cov_diag[i]));
    }

    if (rough_multiplier < 999999999) { //TODO make this an actual proper test
      //add the sats
      matrix_inverse(num_dds_to_add, Z, Z_inv_);

      MAT_PRINTF(Z, num_dds_to_add, num_dds_to_add);
      MAT_PRINTF(Z_inv_, num_dds_to_add, num_dds_to_add);

      for (i=0; i < num_dds_to_add; i++) {
        for (j=0; j < num_dds_to_add; j++) {
          Z_inv[i*num_dds_to_add + j] = (s32) lround(Z_inv_[i*num_dds_to_add + j]);
        }
      }

      memset(decor_float_mean, 0, num_dds_to_add * sizeof(double));
      for (i=0; i < num_dds_to_add; i++) {
        for (j=0; j < num_dds_to_add; j++) {
          decor_float_mean[i] += Z[i*num_dds_to_add + j] * N_mean[j];
        }
        printf("decor_float_mean[%u] = %f\n", i, decor_float_mean[i]);
      }

      add_sats(amb_test, num_dds_to_add, new_dd_prns, decor_float_mean, decor_float_cov_diag, Z_inv);
      //TODO then still need to update amb_test->sats
      break;
    }
    else {
      //try to add fewer sats
      num_dds_to_add -= 1;
    }

    
  }
  return 1;





}

u8 ambiguity_update_sats(ambiguity_test_t *amb_test, u8 num_sdiffs, sdiff_t *sdiffs,
                           double *float_mean, double *float_cov)
{
  //if the sats are the same, we're good
  u8 changed_sats = 0;
  if (!sats_match(amb_test, num_sdiffs, sdiffs)) {
    sdiff_t sdiffs_with_ref_first[num_sdiffs];
    if (ambiguity_update_reference(amb_test, num_sdiffs, sdiffs, sdiffs_with_ref_first)) {
      changed_sats=1;
    }
    u8 intersection_ndxs[num_sdiffs];
    u8 num_dds_in_intersection = 0; //TODO do it proper
    // u8 num_dds_in_intersection = ambiguity_order_sdiffs_with_intersection(amb_test, sdiffs, float_cov, intersection_ndxs);
    if (ambiguity_sat_projection(amb_test, num_dds_in_intersection, intersection_ndxs)) {
      changed_sats = 1;
    }
    sats_management_t float_sats;
    if (ambiguity_sat_inclusion(amb_test, num_dds_in_intersection, intersection_ndxs,
                            &float_sats, float_mean, float_cov)) {
      changed_sats = 1;
    }
    //on exit, ambiguity_sdiffs should match amb_test->sats
  }
  return changed_sats;
  /* TODO: Should we order by 'goodness'? Perhaps by volume of hyps? */

  /* Find out which sats have been added and construct float_mean and float_cov
   * for just the added ones. Rebase float_mean to use our reference (KF may be different!)*/

}

typedef struct {
  s32 upper_bounds[MAX_CHANNELS-1];
  s32 lower_bounds[MAX_CHANNELS-1];
  s32 counter[MAX_CHANNELS-1];
  u8 ndxs_of_old_in_new[MAX_CHANNELS-1];
  u8 ndxs_of_added_in_new[MAX_CHANNELS-1];
  u8 num_added_dds;
  u8 num_old_dds;
  s32 Z_inv[(MAX_CHANNELS-1) * (MAX_CHANNELS-1)];
} generate_hypothesis_state_t;

s8 generate_next_hypothesis(void *x_, u32 n)
{
  (void) n;
  generate_hypothesis_state_t *x = (generate_hypothesis_state_t *)x_;

  printf("[");
  for (u8 i=0; i<x->num_added_dds; i++) {
    printf("%d, ", x->counter[i]);
  }
  printf("]\n");

  if (memcmp(x->upper_bounds, x->counter, x->num_added_dds * sizeof(s32)) == 0) {
    /* counter has reached upper_bound, terminate iteration. */
    return 0;
  }

  for (u8 i=0; i<x->num_added_dds; i++) {
    x->counter[i]++;
    if (x->counter[i] > x->upper_bounds[i]) {
      /* counter[i] has reached maximum, reset counter[i]
       * to lower[i] and 'carry' to next 'digit' */
      x->counter[i] = x->lower_bounds[i];
    } else {
      /* Incremented, so now we have the next counter value. */
      break;
    }
  }

  return 1;
}

void hypothesis_prod(element_t *new_, void *x_, u32 n, element_t *elem_)
{
  (void) elem_;
  (void) n;
  generate_hypothesis_state_t *x = (generate_hypothesis_state_t *)x_;
  hypothesis_t *new = (hypothesis_t *) new_;
  
  u8 *ndxs_of_old_in_new = x->ndxs_of_old_in_new;
  u8 *ndxs_of_added_in_new = x->ndxs_of_added_in_new;

  s32 old_N[MAX_CHANNELS-1];
  memcpy(old_N, new->N, x->num_old_dds * sizeof(s32));

  // printf("counter: [");
  // for (u8 i=0; i < x->num_added_dds; i++) {
  //   printf("%d, ", x->counter[i]);
  // }
  // printf("]\n");
  // printf("old_N:ndx [");
  // for (u8 i=0; i < x->num_old_dds; i++) {
  //   printf("%d, ", ndxs_of_old_in_new[i]);
  // }
  // printf("]\n");
  // printf("old_N:    [");
  // for (u8 i=0; i < x->num_old_dds; i++) {
  //   printf("%d, ", old_N[i]);
  // }
  // printf("]\n");

  // printf("added_N:ndx [");
  // for (u8 i=0; i < x->num_added_dds; i++) {
  //   printf("%d, ", ndxs_of_added_in_new[i]);
  // }
  // printf("\n");
  // printf("added_N:    [");
  // for (u8 i=0; i < x->num_added_dds; i++) {
  //   printf("%d, ", x->counter[i]);
  // }
  // printf("]\n");


  for (u8 i=0; i < x->num_old_dds; i++) {
    new->N[ndxs_of_old_in_new[i]] = old_N[i];
  }
  for (u8 i=0; i<x->num_added_dds; i++) {
    new->N[ndxs_of_added_in_new[i]] = 0;
    for (u8 j=0; j<x->num_added_dds; j++) {
      new->N[ndxs_of_added_in_new[i]] += x->Z_inv[i*x->num_added_dds + j] * x->counter[x->num_old_dds + j];
    }
  }
  // for (u8 i=0; i < x->num_added_dds; i++) {
  //   new->N[ndxs_of_added_in_new[i]] = x->counter[i];
  // }
  // printf("new_N:    [");
  // for (u8 i=0; i < x->num_added_dds + x->num_old_dds; i++) {
  //   printf("%d, ", new->N[i]);
  // }
  // printf("]\n\n");

  /* NOTE: new->ll remains the same as elem->ll as p := exp(ll) is invariant under a
   * constant multiplicative factor common to all hypotheses. TODO: reference^2 document (currently lives in page 3/5.6/2014 of ian's notebook) */
}

typedef struct {
  u8 num_added_dds;
  u8 num_old_dds;
  s32 Z_inv[(MAX_CHANNELS-1)*(MAX_CHANNELS-1)];
} recorrelation_params_t;

void recorrelate_added_sats(void *arg, element_t *elem_)
{
  hypothesis_t *elem = (hypothesis_t *) elem_;
  recorrelation_params_t *params = (recorrelation_params_t *)arg;

  /* elem->N <- [[I 0] [0 Z_inv]] . elem->N
   * where Z_inv is the inverse Lambda reduction matrix having
   * shape (num_added_dds, num_added_dds) and I has shape
   * (num_old_dds, num_old_dds) */

  s32 recorrelated_N[params->num_added_dds];
  memset(recorrelated_N, 0, params->num_added_dds * sizeof(s32));
  for (u8 i=0; i<params->num_added_dds; i++) {
    for (u8 j=0; j<params->num_added_dds; j++) {
      recorrelated_N[i] += params->Z_inv[i*params->num_added_dds + j] * elem->N[params->num_old_dds + j];
    }
  }
  memcpy(&elem->N[params->num_old_dds], recorrelated_N, params->num_added_dds * sizeof(s32));
}

void print_hyp(void *arg, element_t *elem)
{
  u8 num_dds = *( (u8 *) arg );

  hypothesis_t *hyp = (hypothesis_t *)elem;
  printf("[");
  for (u8 i=0; i< num_dds; i++) {
    printf("%d, ", hyp->N[i]);
  }
  printf("]: %f\n", hyp->ll);
}

double add_to_truncate_properly(double x) {
  if (x < 0) {
    return x-0.5;
  }
  else {
    return x + 0.5;
  }
}

void add_sats(ambiguity_test_t *amb_test, 
              u32 num_added_dds, u8 *added_prns,
              double *float_mean, double *float_cov_diag,
              s32 *Z_inv)
{
  /* Make a generator that iterates over the new hypotheses. */
  generate_hypothesis_state_t x0;
  for (u8 i=0; i<num_added_dds; i++) {
    double search_distance = NUM_SEARCH_STDS * sqrt(float_cov_diag[i]);
    // x0.upper_bounds[i] = MAX(floor(float_mean[i] + search_distance), ceil(float_mean[i]));
    // x0.lower_bounds[i] = MIN(ceil(float_mean[i] - search_distance), floor(float_mean[i]));
    x0.upper_bounds[i] = add_to_truncate_properly(ceil(float_mean[i] + search_distance));
    x0.lower_bounds[i] = add_to_truncate_properly(floor(float_mean[i] - search_distance));
    x0.counter[i] = x0.lower_bounds[i];
  }
  printf("float_mean = [");
  for (u8 i=0; i<num_added_dds; i++) {
    printf("%f, ", float_mean[i]);
  }
  printf("]\n");
  printf("float_dist = [");
  for (u8 i=0; i<num_added_dds; i++) {
    printf("%f, ", NUM_SEARCH_STDS * sqrt(float_cov_diag[i]));
  }
  printf("]\n");
  printf("upper = [");
  for (u8 i=0; i<num_added_dds; i++) {
    printf("%d, ", x0.upper_bounds[i]);
  }
  printf("]\n");
  printf("lower = [");
  for (u8 i=0; i<num_added_dds; i++) {
    printf("%d, ", x0.lower_bounds[i]);
  }
  printf("]\n");

  x0.num_added_dds = num_added_dds;
  x0.num_old_dds = MAX(0,amb_test->sats.num_sats-1);

  //then construct the mapping from the old prn indices into the new, and from the added prn indices into the new
  u8 i = 0;
  u8 j = 0;
  u8 k = 0;
  u8 old_prns[x0.num_old_dds];
  memcpy(old_prns, amb_test->sats.prns, x0.num_old_dds * sizeof(u8));
  while (k < amb_test->sats.num_sats + num_added_dds) {
    if (j == x0.num_added_dds || old_prns[i] < added_prns[j]) {
      x0.ndxs_of_old_in_new[i] = k;
      i++;
      k++;
    } else if (i == x0.num_old_dds || old_prns[i] > added_prns[j]) {
      x0.ndxs_of_added_in_new[j] = k;
      j++;
      k++;
    } else {
      printf("This method is being used improperly. This shouldn't happen.\n");
      printf("old_prns = [");
      for (u8 ii=0; ii < x0.num_old_dds; ii++) {
        printf("%d, ",old_prns[ii]);
      }
      printf("]\n");
      printf("added_prns = [");
      for (u8 jj=0; jj < x0.num_old_dds; jj++) {
        printf("%d, ", added_prns[jj]);
      }
      printf("]\n");
      break;
    }
  }

  memcpy(x0.Z_inv, Z_inv, num_added_dds * num_added_dds * sizeof(s32));
  /* Take the product of our current hypothesis state with the generator, recorrelating the new ones as we go. */
  memory_pool_product_generator(amb_test->pool, &x0, MAX_HYPOTHESES, sizeof(x0),
                                &generate_next_hypothesis, &hypothesis_prod);

  // /* Recorrelate the new ambiguities we just added. */
  // recorrelation_params_t params;
  // params.num_added_dds = num_added_dds;
  // params.num_old_dds = MAX(0,amb_test->sats.num_sats-1);
  // memcpy(params.Z_inv, Z_inv, num_added_dds * num_added_dds * sizeof(s32));
  // memory_pool_map(amb_test->pool, &params, &recorrelate_added_sats);
  
  u8 num_all_dds = x0.num_added_dds + x0.num_old_dds;
  memory_pool_map(amb_test->pool, &num_all_dds, &print_hyp);

}

void init_residual_matrices(residual_mtxs_t *res_mtxs, u8 num_dds, double *DE_mtx, double *obs_cov)
{
  res_mtxs->res_dim = 2 * num_dds - 3;
  res_mtxs->null_space_dim = num_dds - 3;
  assign_phase_obs_null_basis(num_dds, DE_mtx, res_mtxs->null_projector);
  assign_residual_covariance_inverse(num_dds, obs_cov, res_mtxs->null_projector, res_mtxs->half_res_cov_inv);
}


void QR_part1(integer m, integer n, double *A, double *tau)
{
  double w[1];
  integer lwork = -1;
  integer info;
  integer jpvt[3];
  memset(jpvt, 0, 3 * sizeof(integer));
  dgeqp3_(&m, &n,
          A, &m,
          jpvt,
          tau,
          w, &lwork, &info);
  lwork = round(w[0]);
  double work[lwork];
  dgeqp3_(&m, &n,
          A, &m,
          jpvt,
          tau,
          work, &lwork, &info); //set A = QR(A)
}

void QR_part2(integer m, integer n, double *A, double *tau)
{
  double w[1];
  integer lwork = -1;
  integer info;
  dorgqr_(&m, &m, &n,
          A, &m,
          tau,
          w, &lwork, &info);
  lwork = round(w[0]);
  double work[lwork];
  dorgqr_(&m, &m, &n,
          A, &m,
          tau,
          work, &lwork, &info);
}

void assign_phase_obs_null_basis(u8 num_dds, double *DE_mtx, double *q)
{
  // use either GEQRF or GEQP3. GEQP3 is the version with pivoting
  // int dgeqrf_(__CLPK_integer *m, __CLPK_integer *n, __CLPK_doublereal *a, __CLPK_integer *
  //       lda, __CLPK_doublereal *tau, __CLPK_doublereal *work, __CLPK_integer *lwork, __CLPK_integer *info)
  // int dgeqp3_(__CLPK_integer *m, __CLPK_integer *n, __CLPK_doublereal *a, __CLPK_integer *
  //       lda, __CLPK_integer *jpvt, __CLPK_doublereal *tau, __CLPK_doublereal *work, __CLPK_integer *lwork,
  //        __CLPK_integer *info)

  //DE is num_sats-1 by 3, need to transpose it to column major
  double A[num_dds * num_dds];
  for (u8 i=0; i<num_dds; i++) {
    for (u8 j=0; j<3; j++) {
      A[j*num_dds + i] = DE_mtx[i*3 + j]; //set A = Transpose(DE_mtx)
    }
  }
  integer m = num_dds;
  integer n = 3;
  double tau[3];
  QR_part1(m, n, A, tau);
  QR_part2(m, n, A, tau);
  memcpy(q, &A[3*num_dds], (num_dds-3) * num_dds * sizeof(double));
}

void assign_residual_covariance_inverse(u8 num_dds, double *obs_cov, double *q, double *r_cov_inv) //TODO make this more efficient (e.g. via page 3/6.2-3/2014 of ian's notebook)
{
  integer dd_dim = 2*num_dds;
  integer res_dim = dd_dim - 3;
  u32 nullspace_dim = num_dds-3;
  double q_tilde[res_dim * dd_dim];
  memset(q_tilde, 0, res_dim * dd_dim * sizeof(double));
  // MAT_PRINTF(obs_cov, dd_dim, dd_dim);
  
  // MAT_PRINTF(q, nullspace_dim, num_dds);
  for (u8 i=0; i<nullspace_dim; i++) {
    memcpy(&q_tilde[i*dd_dim], &q[i*num_dds], num_dds * sizeof(double));
  }
  // MAT_PRINTF(q_tilde, res_dim, dd_dim);
  for (u8 i=0; i<num_dds; i++) {
    q_tilde[(i+nullspace_dim)*dd_dim + i] = 1;
    q_tilde[(i+nullspace_dim)*dd_dim + i+num_dds] = -1 / GPS_L1_LAMBDA_NO_VAC;
  }
  // MAT_PRINTF(q_tilde, res_dim, dd_dim);
  
  //TODO make more efficient via the structure of q_tilde, and it's relation to the I + 1*1^T structure of the obs cov mtx
  double QC[res_dim * dd_dim];
  cblas_dsymm(CblasRowMajor, CblasRight, CblasUpper, //CBLAS_ORDER, CBLAS_SIDE, CBLAS_UPLO
              res_dim, dd_dim, // int M, int N
              1, obs_cov, dd_dim, // double alpha, double *A, int lda
              q_tilde, dd_dim, // double *B, int ldb
              0, QC, dd_dim); // double beta, double *C, int ldc
  // MAT_PRINTF((double *) FC, kf->state_dim, kf->state_dim);
  // MAT_PRINTF(QC, res_dim, dd_dim);

  //TODO make more efficient via the structure of q_tilde, and it's relation to the I + 1*1^T structure of the obs cov mtx
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, // CBLAS_ORDER, CBLAS_TRANSPOSE transA, cBLAS_TRANSPOSE transB
              res_dim, res_dim, dd_dim, // int M, int N, int K
              2, QC, dd_dim, // double alpha, double *A, int lda
              q_tilde, dd_dim, //double *B, int ldb
              0, r_cov_inv, res_dim); //beta, double *C, int ldc
  // MAT_PRINTF(r_cov_inv, res_dim, res_dim);

  // dpotrf_(char *uplo, __CLPK_integer *n, __CLPK_doublereal *a, __CLPK_integer *
  //       lda, __CLPK_integer *info)
  //dpotri_(char *uplo, __CLPK_integer *n, __CLPK_doublereal *a, __CLPK_integer *
  //        lda, __CLPK_integer *info)
  char uplo = 'L'; //actually this makes it upper. the lapack stuff is all column major, but we're good since we're operiating on symmetric matrices
  integer info;
  dpotrf_(&uplo, &res_dim, r_cov_inv, &res_dim, &info);
  // printf("info: %i\n", (int) info);
  // MAT_PRINTF(r_cov_inv, res_dim, res_dim);
  dpotri_(&uplo, &res_dim, r_cov_inv, &res_dim, &info);
  // printf("info: %i\n", (int) info);
  // MAT_PRINTF(r_cov_inv, res_dim, res_dim);
}

void assign_r_vec(residual_mtxs_t *res_mtxs, u8 num_dds, double *dd_measurements, double *r_vec)
{
  cblas_dgemv(CblasRowMajor, CblasNoTrans,
              res_mtxs->null_space_dim, num_dds,
              1, res_mtxs->null_projector, num_dds,
              dd_measurements, 1,
              0, r_vec, 1);
  for (u8 i=0; i< num_dds; i++) {
    r_vec[i + res_mtxs->null_space_dim] = dd_measurements[i] - dd_measurements[i+num_dds] / GPS_L1_LAMBDA_NO_VAC;
  }
}

void assign_r_mean(residual_mtxs_t *res_mtxs, u8 num_dds, double *hypothesis, double *r_mean)
{
  cblas_dgemv(CblasRowMajor, CblasNoTrans,
              res_mtxs->null_space_dim, num_dds,
              1, res_mtxs->null_projector, num_dds,
              hypothesis, 1,
              0, r_mean, 1);
  memcpy(&r_mean[res_mtxs->null_space_dim], hypothesis, num_dds * sizeof(double));
}

double get_quadratic_term(residual_mtxs_t *res_mtxs, u8 num_dds, double *hypothesis, double *r_vec)
{
  double r[res_mtxs->res_dim];
  assign_r_mean(res_mtxs, num_dds, hypothesis, r);
  for (u32 i=0; i<res_mtxs->res_dim; i++) {
    r[i] = r_vec[i] - r[i];
  }
  double half_sig_dot_r[res_mtxs->res_dim];
  cblas_dsymv(CblasRowMajor, CblasUpper,
                 res_mtxs->res_dim,
                 1, res_mtxs->half_res_cov_inv, res_mtxs->res_dim,
                 r, 1,
                 0, half_sig_dot_r, 1);
  double quad_term = 0;
  for (u32 i=0; i<res_mtxs->res_dim; i++) {
    quad_term -= half_sig_dot_r[i] * r[i];
  }
  return quad_term;
}


