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
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <cblas.h>
#include <clapack.h>

#include "logging.h"
#include "constants.h"
#include "baseline.h"
#include "amb_kf.h"
#include "linear_algebra.h"

/** \defgroup baseline Baseline calculations
 * Functions for relating the baseline vector with carrier phase observations
 * and ambiguities.
 * \{ */

/** Predict carrier phase double difference observations from baseline and
 * ambiguity vectors.
 *
 * The carrier phase double difference observations can be modelled as:
 * \f[
 *    \nabla \Delta \phi_i = \frac{1}{\lambda} \mathbf{DE}_i \cdot \mathbf{b} +
 *                           N_i
 * \f]
 *
 * where \f$ \nabla \Delta \phi_i \f$ is the double differenced carrier phase
 * between satellite \f$i\f$ and reference satellite \f$r\f$, \f$N_i \in
 * \mathbb{R}\f$ is the corresponding carrier phase ambiguity.
 *
 * The \f$\mathbf{DE}\f$ matrix is defined as:
 * \f[
 *    \mathbf{DE}_i = \mathbf{e}_i - \mathbf{e}_r
 * \f]
 *
 * where \f$\mathbf{e}_i\f$ is the unit vector to the \f$i\f$th satellite and
 * \f$\mathbf{b}\f$ is the baseline vector between the rover and reference
 * station.
 *
 * \param num_dds Number of double difference observations
 * \param DE Double differenced matrix of unit vectors to the satellites,
 *           length `3 * num_dds`
 * \param N Carrier phase ambiguity vector, length `num_dds`
 * \param dd_obs Double differenced carrier phase observations in cycles,
 *                length `num_dds`
 * \param b Baseline vector in meters
 */
void predict_carrier_obs(u8 num_dds, const double *N, const double *DE,
                         const double b[3], double *dd_obs)
{
  assert(N != NULL);
  assert(DE != NULL);
  assert(b != NULL);
  assert(dd_obs != NULL);

  for (u8 i=0; i<num_dds; i++) {
    dd_obs[i] = vector_dot(3, &DE[3*i], b) / GPS_L1_LAMBDA_NO_VAC + N[i];
  }
}

/** Estimate the integer ambiguity vector from a double difference carrier
 * phase measurement and a given baseline.
 *
 * Given the double difference carrier phase measurement equation:
 *
 * \f[
 *    \nabla \Delta \phi_i = N_i +
 *      \frac{1}{\lambda} (\mathbf{e}_i - \mathbf{e}_r) \cdot \mathbf{b} +
 *      \epsilon
 * \f]
 *
 * where \f$ \nabla \Delta \phi_i \f$ is the double differenced carrier phase
 * between satellite \f$i\f$ and reference satellite \f$r\f$, \f$N_i \in
 * \mathbb{R}\f$ is the corresponding integer ambiguity, \f$\mathbf{e}_i\f$ is
 * the unit vector to the \f$i\f$th satellite and \f$\mathbf{b}\f$ is the
 * baseline vector between the rover and reference station.
 *
 * We can estimate \f$N_i\f$ given the baseline \f$\mathbf{b}\f$ as follows:
 *
 * \f[
 *    \tilde{N_i} = \mathrm{round}\left(
 *      \nabla \Delta \phi_i -
 *      \frac{1}{\lambda} [\mathbf{DE} \cdot \mathbf{b}]_i
 *    \right)
 * \f]
 *
 * where the \f$\mathbf{DE}\f$ matrix is defined as:
 *
 * \f[
 *    \mathbf{DE}_i = \mathbf{e}_i - \mathbf{e}_r
 * \f]
 *
 * \param num_dds Number of double difference observations
 * \param DE Double differenced matrix of unit vectors to the satellites
 * \param dd_obs Double differenced carrier phase observations in cycles,
 *               length `num_sats - 1`
 * \param b Baseline vector in meters
 * \param N Vector where integer ambiguity estimate will be stored
 */
void amb_from_baseline(u8 num_dds, const double *DE, const double *dd_obs,
                       const double b[3], s32 *N)
{
  assert(N != NULL);
  assert(DE != NULL);
  assert(b != NULL);
  assert(dd_obs != NULL);

  /* Solve for ambiguity vector using the observation equation, i.e.
   *   N_float = dd_obs - DE . b / lambda
   * where N_float is a real valued vector */
  for (u8 i=0; i<num_dds; i++) {
    double N_float = dd_obs[i] -
                     vector_dot(3, &DE[3*i], b) / GPS_L1_LAMBDA_NO_VAC;
    /* Round the value of N_float to estimate the integer valued ambiguity. */
    N[i] = (s32)lround(N_float);
  }
}

/** Calculate least squares baseline solution from a set of double difference
 * carrier phase observations and carrier phase ambiguities.
 *
 * For more details, see lesq_solution_float().
 *
 * \note This function takes integer valued carrier phase ambiguities. For real
 *       valued ambiguities, see lesq_solution_float().
 *
 * \param num_dds Number of double difference observations
 * \param dd_obs  Double differenced carrier phase observations in cycles,
 *                length `num_dds`
 * \param N       Carrier phase ambiguity vector, length `num_dds`
 * \param DE      Double differenced matrix of unit vectors to the satellites,
 *                length `3 * num_dds`
 * \param b       The output baseline in meters.
 * \param resid   The output least squares residuals in cycles.
 * \return         0 on success,
 *                -1 if there were insufficient observations to calculate the
 *                   baseline (the solution was under-constrained),
 *                -2 if an error occurred
 */
s8 lesq_solution_int(u8 num_dds, const double *dd_obs, const s32 *N,
                     const double *DE, double b[3], double *resid)
{
  assert(N != NULL);
  double N_float[num_dds];
  for (u8 i=0; i<num_dds; i++) {
    N_float[i] = N[i];
  }
  return lesq_solution_float(num_dds, dd_obs, N_float, DE, b, resid);
}

/* TODO use the state covariance matrix for a better estimate:
 *    That is, decorrelate and scale the LHS of y = A * x before solving for x.
 *
 *    That is, we assume that y ~ N ( A * x, S), and want to make a maximum
 *    likelihood estimate (asymptotically optimal by the Cramer-Rao bound).
 *    This is the minimizer of (A * x - y)' * S^-1 * (A * x - y).
 *        = min_x  (x' * A' - y') * S^-1 * (A * x - y)
 *        = min_x  x' * A' * S^-1* A * x - x' * A' * S^-1 * y
 *                                       - y' * S^-1 * A * x  +  y' * S^-1 * y
 *    Differentiating by x, this becomes
 *    0   =   2 * A' * S^-1 * A * x - 2 * A' * S^-1 * y  ====>
 *    A' * S^-1 * y   =   A' * S^-1 * A * x              ====>
 *    (A'* S^-1 * A)^-1 * A' * S^-1 * y = x
 *
 *    Decomposing S^(-1) = U'^-1 * D^-1/2 * D^-1/2 * U^-1 and defining
 *    C = D^-1/2 * U^-1 * A, we get
 *    (C' * C)^-1 * C' * D^(-1/2) * U^(-1) * y = x
 *    We know that (C' * C)^-1 * C' is the pseudoinverse of C for all C, so
 *    we're really just solving
 *    D^(-1/2) * U^(-1) * y = C * x = D^(-1/2) * U^(-1) * A * x.
 *
 *
 *    Alternatively, we can show this as
 *    If y ~ N( A * x, S = U * D * U'),
 *        then D^(-1/2) * U^(-1) * y ~  N( D^(-1/2) * U^-1 * A * x, I).
 *
 *    This form brings more light on the fact that:
 *      z = U^(-1) * y ~  N(U^-1 * A * x, D).
 *    This form helps us decide how to solve when d_i ~~ 0.
 *    The z_I for I = {i | d_i ~~ 0} want to be solved exactly, because we are
 *    saying that there is no noise in the measurement. This can reduce the
 *    dimension of the problem we are solving, if we really want to.
 *
 *    Therefore we want to:
 *      Triangularly solve y = U * z, for z = U^(-1) * y
 *                         A = U * B, for C = U^(-1) * A
 *        (perhaps simultaneously)
 *      Scale w = D^(-1/2) * z, being mindful that D might not be full rank.
 *      Scale C = D^(-1/2) * B, being mindful that D might not be full rank.
 *      Perform an ordinary least squares solution for w = C * x.
 *          -If we have zeros (d_i below some threshold), we could just do an
 *              exact solution for those then least squares solve the rest.
 *          -Or we could set those d_i to be at threshold and then invert.
 *
 *    We then have that the covariance matrix for x is (C' * C)^-1.
 *
 *    We also need to determine if this is even worth the extra effort.
 */

/** Calculate least squares baseline solution from a set of double difference
 * carrier phase observations and carrier phase ambiguities.
 *
 * Given the double difference carrier phase measurement equation:
 *
 * \f[
 *    \nabla \Delta \phi_i = N_i +
 *      \frac{1}{\lambda} (\mathbf{e}_i - \mathbf{e}_r) \cdot \mathbf{b} +
 *      \epsilon
 * \f]
 *
 * where \f$ \nabla \Delta \phi_i \f$ is the double differenced carrier phase
 * between satellite \f$i\f$ and reference satellite \f$r\f$, \f$N_i \in
 * \mathbb{R}\f$ is the corresponding integer ambiguity, \f$\mathbf{e}_i\f$ is
 * the unit vector to the \f$i\f$th satellite and \f$\mathbf{b}\f$ is the
 * baseline vector between the rover and reference station.
 *
 * If there are 3 or more double difference carrier phase observations then the
 * baseline can be estimated using a least squares solution:
 *
 * \f[
 *    \tilde{\mathbf{b}} = \underset{\mathbf{b}}{\mathrm{argmin}}
 *      \left\|
 *        \frac{1}{\lambda} \mathbf{DE}_i \cdot \mathbf{b} -
 *        \nabla \Delta \phi_i + N_i
 *      \right\|_{\mathcal{l}_2}
 * \f]
 *
 * where:
 *
 * \f[
 *    \mathbf{DE}_i = \mathbf{e}_i - \mathbf{e}_r
 * \f]
 *
 * \note This function takes real valued carrier phase ambiguities. For integer
 *       valued ambiguities, see lesq_solution_int().
 *
 * \param num_dds Number of double difference observations
 * \param dd_obs  Double differenced carrier phase observations in cycles,
 *                length `num_dds`
 * \param N       Carrier phase ambiguity vector, length `num_dds`
 * \param DE      Double differenced matrix of unit vectors to the satellites,
 *                length `3 * num_dds`
 * \param b       The output baseline in meters.
 * \param resid   The output least squares residuals in cycles.
 * \return         0 on success,
 *                -1 if there were insufficient observations to calculate the
 *                   baseline (the solution was under-constrained),
 *                -2 if an error occurred
 */
s8 lesq_solution_float(u8 num_dds_u8, const double *dd_obs, const double *N,
                       const double *DE, double b[3], double *resid)
{
  assert(dd_obs != NULL);
  assert(N != NULL);
  assert(DE != NULL);
  assert(b != NULL);

  if (num_dds_u8 < 3) {
    return -1;
  }

  integer num_dds = num_dds_u8;
  double DET[num_dds * 3];
  matrix_transpose(num_dds, 3, DE, DET);

  double phase_ranges[MAX(num_dds,3)];
  for (u8 i=0; i< num_dds; i++) {
    phase_ranges[i] = dd_obs[i] - N[i];
  }

  s32 ldb = (s32) MAX(num_dds,3);
  integer jpvt[3] = {0, 0, 0};
  double rcond = 1e-12;
  s32 rank;
  s32 info;
  s32 three = 3;
  s32 one = 1;

  /* From LAPACK DGELSY documentation:
   * The unblocked strategy requires that:
   *   LWORK >= MAX( MN+3*N+1, 2*MN+NRHS )
   *   where MN = MIN( M, N )
   *
   * Therefore:
   *   M >= 3, N = 3, NRHS = 1
   *   MN = 3
   *   LWORK >= 13
   */
  s32 lwork = 13;
  double work[lwork];

  /* DGELSY solves:
   *   argmin || A.x - B ||
   * under the l2 norm, where
   *   A <- DE
   *   B <- phase_ranges = dd_obs - N
   *   M <- num_dds
   *   N <- 3
   *   NRHS <- 1
   *
   * the baseline result x is returned in the first 3 elements of phase_ranges.
   */
  dgelsy_(&num_dds, &three, &one, /* M, N, NRHS. */
          DET, &num_dds,          /* A, LDA. */
          phase_ranges, &ldb,     /* B, LDB. */
          jpvt, &rcond,           /* JPVT, RCOND. */
          &rank,                  /* RANK. */
          work, &lwork,           /* WORK, LWORK. */
          &info);                 /* INFO. */

  if (info != 0) {
    log_error("dgelsy returned error %"PRId32"\n", info);
    return -2;
  }

  b[0] = phase_ranges[0] * GPS_L1_LAMBDA_NO_VAC;
  b[1] = phase_ranges[1] * GPS_L1_LAMBDA_NO_VAC;
  b[2] = phase_ranges[2] * GPS_L1_LAMBDA_NO_VAC;

  if (resid) {
    /* Calculate Least Squares Residuals */

    /* resid <= dd_obs - N
     * alpha <= - 1.0 / GPS_L1_LAMBDA_NO_VAC
     * beta <= 1.0
     * resid <= beta * resid + alpha * (DE . b)
     */
    for (u8 i=0; i<num_dds; i++) {
      resid[i] = dd_obs[i] - N[i];
    }
    cblas_dgemv(
      CblasRowMajor, CblasNoTrans, num_dds, 3,
      -1.0 / GPS_L1_LAMBDA_NO_VAC, DE, 3, b, 1,
      1.0, resid, 1
    );
  }

  return 0;
}

/** A least squares solution for baseline from phases using the KF state.
 * This uses the current state of the KF and a set of phase observations to
 * solve for the current baseline.
 *
 * \param kf                    The Kalman filter struct.
 * \param sdiffs_with_ref_first A list of sdiffs. The first in the list must be
 *                              the reference sat of the KF, and the rest must
 *                              correspond to the KF's DD amb estimates' sats,
 *                              sorted in ascending PRN order.
 * \param dd_measurements       A vector of carrier phases. They must be double
 *                              differenced and ordered according to the sdiffs
 *                              and KF's sats (which must match each other).
 * \param ref_ecef              The reference position in ECEF frame, for
 *                              computing the sat direction vectors.
 * \param b                     The output baseline in meters.
 */
void least_squares_solve_b_external_ambs(u8 num_dds_u8, const double *state_mean,
         const sdiff_t *sdiffs_with_ref_first, const double *dd_measurements,
         const double ref_ecef[3], double b[3])
{
  DEBUG_ENTRY();

  integer num_dds = num_dds_u8;
  double DE[num_dds * 3];
  assign_de_mtx(num_dds+1, sdiffs_with_ref_first, ref_ecef, DE);

  s8 ret = lesq_solution_float(num_dds_u8, dd_measurements, state_mean,
                               DE, b, 0);
  (void)ret;
  DEBUG_EXIT();
}

/** \} */

