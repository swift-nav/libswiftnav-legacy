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

#include <math.h>
#include <string.h>
#include <cblas.h>
#include <stdio.h>

#include "constants.h"
#include "stupid_filter.h"
#include "float_kf.h"
#include "linear_algebra.h"

#include <lapacke.h>

/** Estimate the integer ambiguity vector from a double difference measurement
 * and a given baseline.
 *
 * Given the double difference carrier phase measurement equation:
 * \f[
 *    \Delta \nabla \phi_i = N_i + \frac{1}{\lambda} (\mathbf{e}_i - \mathbf{e}_r) \cdot \mathbf{b} + \epsilon
 * \f]
 * where \f$ \Delta \nabla \phi_i \f$ is the double differenced carrier phase
 * between satellite \f$i\f$ and reference satellite \f$r\f$, \f$N_i \in
 * \mathbb{R}\f$ is the corresponding integer ambiguity, \f$\mathbf{e}_i\f$ is the
 * unit vector to the \f$i\f$th satellite and \f$\mathbf{b}\f$ is the baseline
 * vector between the reover and reference station.
 *
 * We can estimate \f$N_i\f$ given the baseline \f$\mathbf{b}\f$ as follows:
 * \f[
 *    \tilde{N_i} = \mathrm{round}\left(\Delta \nabla \phi_i - \frac{1}{\lambda} [\mathbf{DE} \cdot \mathbf{b}]_i + \epsilon\right)
 * \f]
 * where the \f$\mathbf{DE}\f$ matrix is defined as:
 * \f[
 *    \mathbf{DE}_i = \mathbf{e}_i - \mathbf{e}_r
 * \f]
 *
 * \param num_sats Number of satellites used
 * \param DE Double differenced matrix of unit vectors to the satellites
 * \param dd_meas Double differenced carrier phase measurements in cycles,
 *                length `num_sats - 1`
 * \param b Baseline vector in meters
 * \param N Vector where integer ambiguity estimate will be stored
 */
void amb_from_baseline(u8 num_sats, double *DE, double *dd_meas,
                       double b[3], s32 *N)
{
  double N_float[num_sats-1];
  /* Solve for ambiguity vector using the observation equation, i.e.
   *   N_float = dd_meas - DE . b / lambda
   * where N_float is a real valued vector */

  /* N_float <= dd_meas
   * alpha <= - 1.0 / GPS_L1_LAMBDA_NO_VAC
   * beta <= 1.0
   * N_float <= beta * N_float + alpha * (DE . b)
   */
  memcpy(N_float, dd_meas, (num_sats-1) * sizeof(double));
  cblas_dgemv(
    CblasRowMajor, CblasNoTrans, num_sats-1, 3,
    -1.0 / GPS_L1_LAMBDA_NO_VAC, DE, 3, b, 1,
    1.0, N_float, 1
  );

  /* Round the values of N_float to estimate the integer valued ambiguities. */
  for (u8 i=0; i<num_sats-1; i++) {
    N[i] = (s32)lround(N_float[i]);
  }
}

void init_stupid_filter(stupid_filter_state_t *s, u8 num_sats, sdiff_t *sdiffs,
                        double *dd_meas, double b[3], double ref_ecef[3])
{
  double DE[(num_sats-1)*3];

  /* Calculate DE matrix */
  assign_de_mtx(num_sats, sdiffs, ref_ecef, DE);

  amb_from_baseline(num_sats, DE, dd_meas, b, s->N);
}


static s32 find_index_of_element_in_u8s(u32 num_elements, u8 x, u8 *list)
{
  for (u32 i=0; i<num_elements; i++) {
    if (x == list[i]) {
      return i;
    }
  }
  return -1;
}


static void rebase_N(s32 *N, u8 num_sats, u8 *old_prns, u8 *new_prns)
{
  u8 old_ref = old_prns[0];
  u8 new_ref = new_prns[0];

  s32 new_N[num_sats-1];
  s32 index_of_new_ref_in_old = find_index_of_element_in_u8s(num_sats, new_ref, &old_prns[1]);
  s32 val_for_new_ref_in_old_basis = N[index_of_new_ref_in_old];
  for (u8 i=0; i<num_sats-1; i++) {
    u8 new_prn = new_prns[1+i];
    if (new_prn == old_ref) {
      new_N[i] = - val_for_new_ref_in_old_basis;
    }
    else {
      s32 index_of_this_sat_in_old_basis = find_index_of_element_in_u8s(num_sats, new_prn, &old_prns[1]);
      new_N[i] = N[index_of_this_sat_in_old_basis] - val_for_new_ref_in_old_basis;
    }
  }
  memcpy(N, new_N, (num_sats-1) * sizeof(s32));
}

void rebase_stupid_filter(stupid_filter_state_t *s, u8 num_sats, u8 *old_prns, u8 *new_prns)
{
  rebase_N(&(s->N[0]), num_sats, old_prns, new_prns);
}

//assumes both sets are ordered
u8 intersect_o_tron(u8 num_sats1, u8 num_sats2, u8 *sats1, sdiff_t *sdiffs, double *dd_measurements,
                  sdiff_t *intersection_sats, double *intersection_dd_measurements,
                  s32* N, s32 *intersection_N)
{
  u8 i, j, n = 0;

  /* Loop over sats1 and sdffs and check if a PRN is present in both. */
  for (i=0, j=0; i<num_sats1 && j<num_sats2; i++, j++) {
    if (sats1[i] < sdiffs[j].prn)
      j--;
    else if (sats1[i] > sdiffs[j].prn)
      i--;
    else {
      memcpy(&intersection_sats[n], &sdiffs[j], sizeof(sdiff_t));
      intersection_dd_measurements[n] = dd_measurements[j];
      intersection_N[n] = N[i];
      n++;
    }
  }
  return n;
}

void update_sats_stupid_filter(stupid_filter_state_t *s, u8 num_old, u8 *old_prns, u8 num_new,
                               sdiff_t *sdiffs, double *dd_measurements, double ref_ecef[3])
{
  // find intersection of old and new

  double intersection_dd_measurements[num_new];
  sdiff_t intersection_sats[num_new];
  memcpy(&intersection_sats[0], &sdiffs[0], sizeof(sdiff_t));
  s32 intersection_N[MAX_CHANNELS];
  u8 n_intersection = intersect_o_tron(num_old-1, num_new-1, &old_prns[1], &sdiffs[1], dd_measurements, &intersection_sats[1], intersection_dd_measurements, s->N, intersection_N);

  if ((num_old-1) == (num_new-1) && (num_old-1) == n_intersection) {
    u8 flag = 0;
    for (u8 i=0; i<n_intersection+1; i++) {
      if (intersection_sats[i].prn != old_prns[i])
        flag = 1;
    }
    if (flag == 0)
      // No changes!
      return;
  }

  printf("====== UPDATE =======\n");

  // ok to overwite s->N because we are going to call init in a second anyway
  memcpy(&(s->N), intersection_N, n_intersection*sizeof(s32));
  // calc least sq. b from intersection sat using new data
  double b[3];
  update_stupid_filter(s, n_intersection, intersection_sats, intersection_dd_measurements, b, ref_ecef);

  // call init with the new sat set and b as initial baseline
  init_stupid_filter(s, num_new, sdiffs, dd_measurements, b, ref_ecef);
}

void update_stupid_filter(stupid_filter_state_t *s, u8 num_sats, sdiff_t *sdiffs,
                        double *dd_measurements, double b[3], double ref_ecef[3])
{
  double DE[(num_sats-1)*3];

  /* Calculate DE matrix */
  assign_de_mtx(num_sats, sdiffs, ref_ecef, DE);

  /* Solve for b via least squares, i.e.
   * dd_meas = DE . b + N 
   *  =>  DE . b = (dd_meas - N) * lambda */

  /* min | A.x - b | wrt x
   * A <= DE
   * x <= b
   * b <= (dd_meas - N) * lambda
   */
   double rhs[MAX(num_sats-1, 3)];
   for (u8 i=0; i<num_sats-1; i++) {
     rhs[i] = (dd_measurements[i] - s->N[i]) * GPS_L1_LAMBDA_NO_VAC;
   }
  VEC_PRINTF(rhs, (u32) num_sats-1);
   int jpvt[3] = {0, 0, 0};
   int rank;
   LAPACKE_dgelsy(LAPACK_ROW_MAJOR, num_sats-1, 3,
                  1, DE, 3,
                  rhs, 1, jpvt,
                  -1, &rank);
   memcpy(b, rhs, 3*sizeof(double));

  /* Calculate DE matrix */
  assign_de_mtx(num_sats, sdiffs, ref_ecef, DE);

  /* Solve for ambiguity vector, i.e.
   * N = dd_meas - DE . b / lambda */
  double b_dot_DE[num_sats-1];
  cblas_dgemv(CblasRowMajor, CblasNoTrans, // CBLAS_ORDER, CBLAS_TRANSPOSE
            num_sats-1, 3, // int M, int N,
            1, DE, 3, // double alpha, double *A, int lda
            rhs, 1, // double *X, int incX
            0, b_dot_DE, 1); // double beta, double *Y, int incY
  VEC_PRINTF(b_dot_DE, (u32) num_sats-1);

  for (u8 i=0; i<num_sats-1; i++) {
    b_dot_DE[i] -= (dd_measurements[i] - s->N[i]) * GPS_L1_LAMBDA_NO_VAC;
    b_dot_DE[i] = 100 * b_dot_DE[i];
  }
  VEC_PRINTF(b_dot_DE, (u32) num_sats-1);
}











