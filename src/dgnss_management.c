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

#include <string.h>
#include <stdio.h>
#include "amb_kf.h"
#include "float_kf.h"
#include "stupid_filter.h"
#include "single_diff.h"
#include "dgnss_management.h"
#include "linear_algebra.h"
#include "ambiguity_test.h"

nkf_t nkf;
kf_t kf;
stupid_filter_state_t stupid_state;
sats_management_t sats_management;
ambiguity_test_t ambiguity_test;

dgnss_settings_t dgnss_settings = {
  .phase_var_test = DEFAULT_PHASE_VAR_TEST,
  .code_var_test = DEFAULT_CODE_VAR_TEST,
  .phase_var_kf = DEFAULT_PHASE_VAR_KF,
  .code_var_kf = DEFAULT_CODE_VAR_KF,
  .pos_trans_var = DEFAULT_POS_TRANS_VAR,
  .vel_trans_var = DEFAULT_VEL_TRANS_VAR,
  .int_trans_var = DEFAULT_INT_TRANS_VAR,
  .pos_init_var = DEFAULT_POS_INIT_VAR,
  .vel_init_var = DEFAULT_VEL_INIT_VAR,
  .amb_init_var = DEFAULT_AMB_INIT_VAR,
  .new_int_var = DEFAULT_NEW_INT_VAR,
};

void make_measurements(u8 num_double_diffs, sdiff_t *sdiffs, double *raw_measurements)
{
  double phase0 = sdiffs[0].carrier_phase;
  double code0 = sdiffs[0].pseudorange;
  for (u8 i=0; i<num_double_diffs; i++) {
    raw_measurements[i] = sdiffs[i+1].carrier_phase - phase0;
    raw_measurements[i+num_double_diffs] = sdiffs[i+1].pseudorange - code0;
  }
}

bool prns_match(u8 *old_non_ref_prns, u8 num_non_ref_sdiffs, sdiff_t *non_ref_sdiffs)
{
  if (sats_management.num_sats-1 != num_non_ref_sdiffs) {
    // printf("len_doesn't match\n");
    return false;
  }
  for (u8 i=0; i<num_non_ref_sdiffs; i++) {
    //iterate through the non-reference_sats
    /*printf("old[%u]=%u, new[%u]=%u\n", i+1, old_non_ref_prns[i], i+1, non_ref_sdiffs[i].prn);*/
    if (old_non_ref_prns[i] != non_ref_sdiffs[i].prn) {
      return false;
    }
  }
  /*printf("prns_match\n");*/
  return true;
}

/** Finds the prns of the intersection between old prns and new measurements.
 * It returns the length of the intersection
 */
u8 dgnss_intersect_sats(u8 num_old_prns, u8 *old_prns,
                  u8 num_sdiffs, sdiff_t *sdiffs,
                  u8 *ndx_of_intersection_in_old,
                  u8 *ndx_of_intersection_in_new)
{
  u8 i, j, n = 0;
  /* Loop over old_prns and sdiffs and check if a PRN is present in both. */
  for (i=0, j=0; i<num_old_prns && j<num_sdiffs; i++, j++) {
    if (old_prns[i] < sdiffs[j].prn)
      j--;
    else if (old_prns[i] > sdiffs[j].prn)
      i--;
    else {
      ndx_of_intersection_in_old[n] = i;
      ndx_of_intersection_in_new[n] = j;
      n++;
    }
  }
  return n;
}

void dgnss_init(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3], double dt)
{
  sdiff_t corrected_sdiffs[num_sats];
  init_sats_management(&sats_management, num_sats, sdiffs, corrected_sdiffs);

  double dd_measurements[2*(num_sats-1)];
  make_measurements(num_sats-1, corrected_sdiffs, dd_measurements);

  get_kf(
    &kf,
    dgnss_settings.phase_var_kf, dgnss_settings.code_var_kf,
    dgnss_settings.pos_trans_var, dgnss_settings.vel_trans_var,
    dgnss_settings.int_trans_var,
    dgnss_settings.pos_init_var, dgnss_settings.vel_init_var,
    dgnss_settings.amb_init_var,
    num_sats, corrected_sdiffs, dd_measurements, reciever_ecef,
    dt
  );

  set_nkf(
    &nkf,
    dgnss_settings.phase_var_kf, dgnss_settings.code_var_kf,
    dgnss_settings.amb_init_var,
    num_sats, corrected_sdiffs, dd_measurements, reciever_ecef
  );

  create_ambiguity_test(&ambiguity_test);
}

void dgnss_rebase_ref(u8 num_sdiffs, sdiff_t *sdiffs, double reciever_ecef[3], u8 old_prns[MAX_CHANNELS], sdiff_t *corrected_sdiffs)
{
  (void)reciever_ecef;
  //all the ref sat stuff
  s8 sats_management_code = rebase_sats_management(&sats_management, num_sdiffs, sdiffs, corrected_sdiffs);
  if (sats_management_code == NEW_REF_START_OVER) {
    printf("====== START OVER =======\n"); // TODO WRITE WHAT GOTTA BE DONE, YO
    /*dgnss_init(num_sdiffs, sdiffs, reciever_ecef, dt); //TODO use current baseline state*/
    return;
  }
  else if (sats_management_code == NEW_REF) {
    // do everything related to changing the reference sat here
    rebase_kf(&kf, sats_management.num_sats, &old_prns[0], &sats_management.prns[0]);
    rebase_nkf(&nkf, sats_management.num_sats, &old_prns[0], &sats_management.prns[0]);
  }
}

void sdiffs_to_prns(u8 n, sdiff_t *sdiffs, u8 *prns)
{
  for (u8 i=0; i<n; i++) {
    prns[i] = sdiffs[i].prn;
  }
}

void dgnss_update_sats(u8 num_sdiffs, double reciever_ecef[3], sdiff_t *corrected_sdiffs,
                       double *dd_measurements, double dt)
{
  (void)dd_measurements;
  u8 new_prns[num_sdiffs];
  sdiffs_to_prns(num_sdiffs, corrected_sdiffs, new_prns);

  u8 old_prns[MAX_CHANNELS];
  memcpy(old_prns, sats_management.prns, sats_management.num_sats * sizeof(u8));

  if (!prns_match(&old_prns[1], num_sdiffs-1, &corrected_sdiffs[1])) {
    u8 ndx_of_intersection_in_old[sats_management.num_sats];
    u8 ndx_of_intersection_in_new[sats_management.num_sats];
    ndx_of_intersection_in_old[0] = 0;
    ndx_of_intersection_in_new[0] = 0;
    u8 num_intersection_sats = dgnss_intersect_sats(
        sats_management.num_sats-1, &old_prns[1],
        num_sdiffs-1, &corrected_sdiffs[1],
        &ndx_of_intersection_in_old[1],
        &ndx_of_intersection_in_new[1]) + 1;

    reset_kf_except_state(
      &kf,
      dgnss_settings.phase_var_kf, dgnss_settings.code_var_kf,
      dgnss_settings.pos_trans_var, dgnss_settings.vel_trans_var,
      dgnss_settings.int_trans_var,
      num_sdiffs, corrected_sdiffs, reciever_ecef, dt
    );

    set_nkf_matrices(
      &nkf,
      dgnss_settings.phase_var_kf, dgnss_settings.code_var_kf,
      num_sdiffs, corrected_sdiffs, reciever_ecef
    );

    if (num_intersection_sats < sats_management.num_sats) { //lost sats
      kalman_filter_state_projection(&kf,
                                     sats_management.num_sats-1,
                                     num_intersection_sats-1,
                                     &ndx_of_intersection_in_old[1]);
      nkf_state_projection(&nkf,
                           sats_management.num_sats-1,
                           num_intersection_sats-1,
                           &ndx_of_intersection_in_old[1]);
    }
    if (num_intersection_sats < num_sdiffs) { //gained sats
      kalman_filter_state_inclusion(&kf,
                                    num_intersection_sats-1,
                                    num_sdiffs-1,
                                    &ndx_of_intersection_in_new[1],
                                    dgnss_settings.new_int_var);
      nkf_state_inclusion(&nkf,
                          num_intersection_sats-1,
                          num_sdiffs-1,
                          &ndx_of_intersection_in_new[1],
                          dgnss_settings.new_int_var);
    }

    /*print_sats_management(&sats_management);*/
    update_sats_sats_management(&sats_management, num_sdiffs-1, &corrected_sdiffs[1]);
    /*print_sats_management(&sats_management);*/

  }
  else {
    reset_kf_except_state(
      &kf,
      dgnss_settings.phase_var_kf, dgnss_settings.code_var_kf,
      dgnss_settings.pos_trans_var, dgnss_settings.vel_trans_var,
      dgnss_settings.int_trans_var,
      num_sdiffs, corrected_sdiffs, reciever_ecef, dt
    );
    set_nkf_matrices(
      &nkf,
      dgnss_settings.phase_var_kf, dgnss_settings.code_var_kf,
      num_sdiffs, corrected_sdiffs, reciever_ecef
    );
  }

}

void dgnss_incorporate_observation(sdiff_t *sdiffs, double * dd_measurements,
                                   double *reciever_ecef, double dt)
{
  (void) dt;

  double b2[3];
  least_squares_solve_b(&nkf, sdiffs, dd_measurements, reciever_ecef, b2);

  double ref_ecef[3];

  ref_ecef[0] = reciever_ecef[0] + 0.5 * b2[0];
  ref_ecef[1] = reciever_ecef[1] + 0.5 * b2[0];
  ref_ecef[2] = reciever_ecef[2] + 0.5 * b2[0];

  /* TODO: make a common DE and use it instead. */
  assign_decor_obs_mtx(sats_management.num_sats, sdiffs, ref_ecef,
                       kf.decor_mtx, kf.decor_obs_mtx);

  set_nkf_matrices(&nkf,
                   dgnss_settings.phase_var_kf, dgnss_settings.code_var_kf,
                   sats_management.num_sats, sdiffs, ref_ecef);

  kalman_filter_update(&kf, dd_measurements);
  nkf_update(&nkf, dd_measurements);
}

void dvec_printf(double *v, u32 n)
{
  for (u32 i = 0; i < n; i++)
    printf(", %f", v[i]);
}

void dgnss_update(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3], double dt)
{
  sdiff_t corrected_sdiffs[num_sats];

  u8 old_prns[MAX_CHANNELS];
  memcpy(old_prns, sats_management.prns, sats_management.num_sats * sizeof(u8));
  //rebase globals to a new reference sat (permutes corrected_sdiffs accordingly)
  dgnss_rebase_ref(num_sats, sdiffs, reciever_ecef, old_prns, corrected_sdiffs);

  double dd_measurements[2*(num_sats-1)];
  make_measurements(num_sats-1, corrected_sdiffs, dd_measurements);

  //all the added/dropped sat stuff
  dgnss_update_sats(num_sats, reciever_ecef, corrected_sdiffs, dd_measurements, dt);
  /*printf("done updating sats\n");*/
  /*MAT_PRINTF(kf.decor_obs_mtx, kf.obs_dim, kf.state_dim);*/

  if (num_sats >= 5) {
    // update for observation
    dgnss_incorporate_observation(corrected_sdiffs, dd_measurements, reciever_ecef, dt);
  }

  double b2[3];
  least_squares_solve_b(&nkf, corrected_sdiffs, dd_measurements, reciever_ecef, b2);

  double ref_ecef[3];
  ref_ecef[0] = reciever_ecef[0] + 0.5 * b2[0];
  ref_ecef[1] = reciever_ecef[1] + 0.5 * b2[1];
  ref_ecef[2] = reciever_ecef[2] + 0.5 * b2[2];

  /*
  update_ambiguity_test(ref_ecef,
                        dgnss_settings.phase_var_test,
                        dgnss_settings.code_var_test,
                        &ambiguity_test,
                        kf.state_dim, &sats_management, sdiffs,
                        kf.state_mean, kf.state_cov_U, kf.state_cov_D);
  */
  update_ambiguity_test(ref_ecef,
                        dgnss_settings.phase_var_test,
                        dgnss_settings.code_var_test,
                        &ambiguity_test, nkf.state_dim, &sats_management,
                        sdiffs, nkf.state_mean, nkf.state_cov_U,
                        nkf.state_cov_D);
}

u32 dgnss_iar_num_hyps(void)
{
  return ambiguity_test_n_hypotheses(&ambiguity_test);
}


void dgnss_float_baseline(u8 *num_used, double b[3])
{
  memcpy(b, kf.state_mean, 3 * sizeof(double));
  *num_used = sats_management.num_sats;
}

void dgnss_new_float_baseline(u8 num_sats, sdiff_t *sdiffs, double receiver_ecef[3], u8 *num_used, double b[3])
{
  sdiff_t corrected_sdiffs[num_sats];

  u8 old_prns[MAX_CHANNELS];
  memcpy(old_prns, sats_management.prns, sats_management.num_sats * sizeof(u8));
  //rebase globals to a new reference sat (permutes corrected_sdiffs accordingly)
  dgnss_rebase_ref(num_sats, sdiffs, receiver_ecef, old_prns, corrected_sdiffs);

  double dd_measurements[2*(num_sats-1)];
  make_measurements(num_sats-1, corrected_sdiffs, dd_measurements);

  least_squares_solve_b(&nkf, corrected_sdiffs, dd_measurements, receiver_ecef, b);
  *num_used = sats_management.num_sats;
}

void dgnss_fixed_baseline(u8 n, sdiff_t *sdiffs, double ref_ecef[3],
                          u8 *num_used, double b[3])
{
  if (dgnss_iar_resolved()) {
    sdiff_t ambiguity_sdiffs[ambiguity_test.sats.num_sats];
    double dd_meas[2*(ambiguity_test.sats.num_sats-1)];
    make_ambiguity_dd_measurements_and_sdiffs(&ambiguity_test, n, sdiffs,
        dd_meas, ambiguity_sdiffs);
    double DE[(ambiguity_test.sats.num_sats-1) * 3];
    assign_de_mtx(ambiguity_test.sats.num_sats, ambiguity_sdiffs, ref_ecef, DE);
    hypothesis_t *hyp = (hypothesis_t*)ambiguity_test.pool->allocated_nodes_head->elem;
    *num_used = ambiguity_test.sats.num_sats;
    lesq_solution(ambiguity_test.sats.num_sats-1, dd_meas, hyp->N, DE, b, 0);
  } else {
    dgnss_new_float_baseline(n, sdiffs, ref_ecef, num_used, b);
  }
}

void dgnss_reset_iar()
{
  create_ambiguity_test(&ambiguity_test);
}

void dgnss_init_known_baseline(u8 num_sats, sdiff_t *sdiffs, double receiver_ecef[3], double b[3])
{
  double ref_ecef[3];
  ref_ecef[0] = receiver_ecef[0] + 0.5 * b[0];
  ref_ecef[1] = receiver_ecef[1] + 0.5 * b[1];
  ref_ecef[2] = receiver_ecef[2] + 0.5 * b[2];

  sdiff_t corrected_sdiffs[num_sats];

  u8 old_prns[MAX_CHANNELS];
  memcpy(old_prns, sats_management.prns, sats_management.num_sats * sizeof(u8));
  //rebase globals to a new reference sat (permutes corrected_sdiffs accordingly)
  dgnss_rebase_ref(num_sats, sdiffs, ref_ecef, old_prns, corrected_sdiffs);

  double dds[2*(num_sats-1)];
  make_measurements(num_sats-1, corrected_sdiffs, dds);

  double DE[(num_sats-1)*3];
  assign_de_mtx(num_sats, corrected_sdiffs, ref_ecef, DE);

  dgnss_reset_iar();

  memcpy(&ambiguity_test.sats, &sats_management, sizeof(sats_management));
  hypothesis_t *hyp = (hypothesis_t *)memory_pool_add(ambiguity_test.pool);
  hyp->ll = 0;
  amb_from_baseline(num_sats, DE, dds, b, hyp->N);

  double obs_cov[(num_sats-1) * (num_sats-1) * 4];
  memset(obs_cov, 0, (num_sats-1) * (num_sats-1) * 4 * sizeof(double));
  u8 num_dds = num_sats-1;
  for (u8 i=0; i<num_dds; i++) {
    for (u8 j=0; j<num_dds; j++) {
      u8 i_ = i+num_dds;
      u8 j_ = j+num_dds;
      if (i==j) {
        obs_cov[i*2*num_dds + j] = dgnss_settings.phase_var_test * 2;
        obs_cov[i_*2*num_dds + j_] = dgnss_settings.code_var_test * 2;
      }
      else {
        obs_cov[i*2*num_dds + j] = dgnss_settings.phase_var_test;
        obs_cov[i_*2*num_dds + j_] = dgnss_settings.code_var_test;
      }
    }
  }

  init_residual_matrices(&ambiguity_test.res_mtxs, num_sats-1, DE, obs_cov);

  /*printf("Known Base: [");*/
  /*for (u8 i=0; i<num_sats-1; i++)*/
    /*printf("%d, ", hyp->N[i]);*/
  /*printf("]\n");*/
}

void dgnss_init_known_baseline2(u8 num_sats, sdiff_t *sdiffs, double receiver_ecef[3], double b[3])
{
  double ref_ecef[3];
  ref_ecef[0] = receiver_ecef[0] + 0.5 * b[0];
  ref_ecef[1] = receiver_ecef[1] + 0.5 * b[1];
  ref_ecef[2] = receiver_ecef[2] + 0.5 * b[2];

  sdiff_t corrected_sdiffs[num_sats];

  u8 old_prns[MAX_CHANNELS];
  memcpy(old_prns, sats_management.prns, sats_management.num_sats * sizeof(u8));
  //rebase globals to a new reference sat (permutes corrected_sdiffs accordingly)
  dgnss_rebase_ref(num_sats, sdiffs, ref_ecef, old_prns, corrected_sdiffs);

  double dds[2*(num_sats-1)];
  make_measurements(num_sats-1, corrected_sdiffs, dds);

  double DE[(num_sats-1)*3];
  s32 N[num_sats-1];
  assign_de_mtx(num_sats, corrected_sdiffs, ref_ecef, DE);
  amb_from_baseline(num_sats, DE, dds, b, N);

  printf("Known Base: [");
  for (u8 i=0; i<num_sats-1; i++)
    printf("%d, ", N[i]);
  printf("]\n");

  /* Construct fake state means. */
  double state_mean[num_sats-1+6];
  memcpy(&state_mean[0], b, 3 * sizeof(double));
  memset(&state_mean[3], 0, 3 * sizeof(double));
  for (u8 i=0; i<num_sats-1; i++) {
    state_mean[i+6] = N[i];
  }

  /* Construct fake covariance U factor (just identity). */
  double state_cov_U[(num_sats-1+6)*(num_sats-1+6)];
  eye(num_sats-1+6, state_cov_U);

  double state_cov_D[num_sats-1+6];
  memset(state_cov_D, 0, 6 * sizeof(double));
  for (u32 i=0; i<num_sats-1; i++) {
    state_cov_D[i+6] = 1.0 / 64.0;
  }

  dgnss_reset_iar();

  update_ambiguity_test(ref_ecef,
                        dgnss_settings.phase_var_test,
                        dgnss_settings.code_var_test,
                        &ambiguity_test,
                        num_sats-1+6, &sats_management, corrected_sdiffs,
                        state_mean, state_cov_U, state_cov_D);
}

s8 dgnss_iar_resolved()
{
  return ambiguity_test_n_hypotheses(&ambiguity_test) == 1;
}

kf_t * get_dgnss_kf()
{
  return &kf;
}

nkf_t * get_dgnss_nkf()
{
  return &nkf;
}

s32 * get_stupid_filter_ints()
{
  return stupid_state.N;
}

sats_management_t * get_sats_management()
{
  return &sats_management;
}


