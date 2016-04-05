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

#include <libswiftnav/logging.h>
#include <libswiftnav/constants.h>
#include <libswiftnav/baseline.h>
#include <libswiftnav/amb_kf.h>
#include <libswiftnav/linear_algebra.h>
#include <libswiftnav/filter_utils.h>
#include <libswiftnav/set.h>
#include <libswiftnav/sats_management.h> /* choose_reference_sat */

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
 * \param num_dds_u8 Number of double difference observations
 * \param dd_obs     Double differenced carrier phase observations in cycles,
 *                   length `num_dds`
 * \param N          Carrier phase ambiguity vector, length `num_dds`
 * \param DE         Double differenced matrix of unit vectors to the satellites,
 *                   length `3 * num_dds`
 * \param b          The output baseline in meters.
 * \param resid      The output least squares residuals in cycles.
 * \return           BASELINE_SUCCESS on success,
 *                   BASELINE_NOT_ENOUGH_SATS_FLOAT if there were insufficient
                        observations to calculate the
 *                      baseline (the solution was under-constrained),
 *                   BASELINE_DGELSY_FAIL if dgelsy_() failed
 */
s8 lesq_solution_float(u8 num_dds_u8, const double *dd_obs,
                                   const double *N, const double *DE,
                                   double b[3], double *resid)
{
  assert(dd_obs != NULL);
  assert(N != NULL);
  assert(DE != NULL);
  assert(b != NULL);

  if (num_dds_u8 < 3) {
    return BASELINE_NOT_ENOUGH_SATS_FLOAT;
  }

  assert(num_dds_u8 < MAX_CHANNELS);

  for(int i = 0; i < num_dds_u8 * 3; i++) {
    assert(isfinite(DE[i]));
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
    log_error("dgelsy returned error %"PRId32"", info);
    return BASELINE_DGELSY_FAIL;
  }

  assert(num_dds == num_dds_u8);

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

  return BASELINE_SUCCESS;
}

static void drop_i(u32 index, u32 len, u32 size, const double *from, double *to)
{
  memcpy(to, from, index*size*sizeof(double));
  memcpy(to + index*size,
         from + (index + 1)*size,
         (len - index - 1)*size*sizeof(double));
}

static s8 lesq_without_i(u8 dropped_dd, u8 num_dds, const double *dd_obs,
                         const double *N, const double *DE, double b[3],
                         double *resid)
{
  assert(num_dds > 3);
  assert(num_dds < MAX_CHANNELS);

  u8 new_dds = num_dds - 1;
  double new_obs[new_dds];
  double new_N[new_dds];
  double new_DE[new_dds * 3];

  drop_i(dropped_dd, num_dds, 1, dd_obs, new_obs);
  drop_i(dropped_dd, num_dds, 1, N, new_N);
  drop_i(dropped_dd, num_dds, 3, DE, new_DE);

  return lesq_solution_float(new_dds, new_obs, new_N, new_DE, b, resid);
}

/** Approximate chi square test
 * Scales least squares residuals by estimated variance,
 *   compares against threshold.
 *
 * \param threshold The test threshold.
 *                  Value 5.5 recommended for float/fixed baseline residuals.
 * \param num_dds Length of residuals.
 * \param residuals Residual vector calculated by lesq_solution_float.
 * \param residual If not null, used to output double value of scaled residual.
 *
 * \return true if norm is below threshold.
 */
static bool chi_test(double threshold, u8 num_dds,
                     const double *residuals, double *residual)
{
  assert(num_dds < MAX_CHANNELS);

  double sigma = DEFAULT_PHASE_VAR_KF;
  double norm = vector_norm(num_dds, residuals) / sqrt(sigma);
  if (residual) {
    *residual = norm;
  }
  return norm < threshold;
}

/** Calculate lesq baseline with raim check/repair
 *
 * \param num_dds_u8 Number of double difference observations
 * \param dd_obs     Double differenced carrier phase observations in cycles,
 *                   length `num_dds_u8`
 * \param N          Carrier phase ambiguity vector, length `num_dds_u8`
 * \param DE         Double differenced matrix of unit vectors to the satellites,
 *                   length `3 * num_dds`
 * \param b          The output baseline in meters.
 * \param disable_raim if true raim check/repair will not be performed
 * \param raim_threshold test threshold for check
 * \param n_used if not null, outputs num obs used in solution
 * \param ret_residuals if not null, returns residual vector
 * \param removed_obs if not null and repair performed,
 *                    returns removed obs index
 * \return
 *    BASELINE_SUCCESS_NO_RAIM: solution ok, but raim check was not available
 *                              (exactly 3 dds, or explicitly disabled)
 *
 *    BASELINE_SUCCESS_RAIM_REPAIR: repaired solution,
 *                                  using one fewer observation
 *                                  returns index of removed observation if
 *                                  removed_obs ptr is passed
 *
 *    BASELINE_SUCCESS: solution with all dd's ok
 *
 *    BASELINE_RAIM_REPAIR_FAIL: raim check failed, repair failed
 *
 *    BASELINE_RAIM_FAIL_NOT_ENOUGH_SATS: raim check failed,
 *                                        not enough sats for repair
 *
 *    BASELINE_RAIM_REPAIR_MULTI_SOLNS: raim check failed,
 *                                      more then one solution
 *
 *    or one of the returns from lesq_solution_float()
 */
/* TODO(dsk) update all call sites to use n_used as calculated here.
 * TODO(dsk) add warn/info logging to call sites when repair occurs.
 * TODO(dsk) make this static
 */
s8 lesq_solve_raim(u8 num_dds_u8, const double *dd_obs,
                   const double *N, const double *DE, double b[3],
                   bool disable_raim, double raim_threshold,
                   u8 *n_used, double *ret_residuals, u8 *removed_obs)
{
  integer num_dds = num_dds_u8;
  double residuals[num_dds];
  double residual;

  assert(num_dds < MAX_CHANNELS);
  assert(num_dds_u8 < MAX_CHANNELS);

  s8 res = lesq_solution_float(num_dds_u8, dd_obs, N, DE, b, residuals);

  if (res != BASELINE_SUCCESS) {
    /* Not enough sats or other error returned by initial lesq solution. */
    return res;
  }

  if (disable_raim || chi_test(raim_threshold, num_dds, residuals, &residual)) {
    /* Solution using all sats ok. */
    if (ret_residuals) {
      memcpy(ret_residuals, residuals, num_dds * sizeof(double));
    }
    if (n_used) {
      *n_used = num_dds;
    }
    if (disable_raim || num_dds == 3) {
      return BASELINE_SUCCESS_NO_RAIM;
    }
    return BASELINE_SUCCESS;
  }

  if (num_dds < 5) {
    /* We have just enough sats for a solution; can't search for solution
     * after dropping one.
     * 5 are needed because a 3 dimensional system is exactly constrained,
     * so the bad measurement can't be detected.
     */
    if (n_used) {
      *n_used = 0;
    }
    return BASELINE_NOT_ENOUGH_SATS_RAIM;
  }

  u8 num_passing = 0;
  u8 bad_sat = -1;
  u8 new_dds = num_dds - 1;

  for (u8 i = 0; i < num_dds; i++) {
    if (0 == lesq_without_i(i, num_dds, dd_obs, N, DE, b, residuals)) {
      if (chi_test(raim_threshold, new_dds, residuals, &residual)) {
        num_passing++;
        bad_sat = i;
      }
    }
  }

  if (num_passing == 1) {
    /* bad_sat holds index of bad dd
     * Return solution without bad_sat. */
    /* Recalculate this solution. */
    (void)lesq_without_i(bad_sat, num_dds, dd_obs, N, DE, b, residuals);
    if (removed_obs) {
      *removed_obs = bad_sat;
    }
    if (ret_residuals) {
      memcpy(ret_residuals, residuals, (num_dds - 1) * sizeof(double));
    }
    if (n_used) {
      *n_used = num_dds - 1;
    }
    return BASELINE_SUCCESS_RAIM_REPAIR;
  } else if (num_passing == 0) {
    /* Ref sat is bad? */
    if (n_used) {
      *n_used = 0;
    }
    return BASELINE_RAIM_REPAIR_FAIL;
  } else {
    /* Had more than one acceptable solution.
     * TODO(dsk) should we return the best one? */
    if (n_used) {
      *n_used = 0;
    }
    return BASELINE_RAIM_REPAIR_MULTI_SOLNS;
  }
}

/** A least squares solution for baseline from phases using the KF state.
 * This uses the current state of the KF and a set of phase observations to
 * solve for the current baseline.
 *
 * \param num_dds_u8            State dimension
 * \param state_mean            KF estimated state mean
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
 * \param disable_raim          True disables raim check/repair
 * \param raim_threshold        Threshold for raim checks.
 * \return                      See lesq_solve_raim()
 */
s8 least_squares_solve_b_external_ambs(u8 num_dds_u8, const double *state_mean,
         const sdiff_t *sdiffs_with_ref_first, const double *dd_measurements,
         const double ref_ecef[3], double b[3],
         bool disable_raim, double raim_threshold)
{
  DEBUG_ENTRY();

  integer num_dds = num_dds_u8;
  double DE[num_dds * 3];
  assign_de_mtx(num_dds+1, sdiffs_with_ref_first, ref_ecef, DE);

  s8 code = lesq_solve_raim(num_dds_u8, dd_measurements, state_mean, DE, b,
                            disable_raim, raim_threshold, 0, 0, 0);
  DEBUG_EXIT();
  return code;
}

/** Comparison function for `ambiguity_t` by PRN.
 * See `cmp_fn`. */
int cmp_amb(const void *a_, const void *b_)
{
  const ambiguity_t *a = (const ambiguity_t *)a_;
  const ambiguity_t *b = (const ambiguity_t *)b_;
  return sid_compare(a->sid, b->sid);
}

int cmp_amb_sdiff(const void *a_, const void *b_)
{
  const ambiguity_t *a = (const ambiguity_t *)a_;
  const sdiff_t *b = (const sdiff_t *)b_;
  return sid_compare(a->sid, b->sid);
}

int cmp_amb_sid(const void *a_, const void *b_)
{
  const ambiguity_t *a = (const ambiguity_t *)a_;
  const gnss_signal_t *b = (const gnss_signal_t *)b_;
  return sid_compare(a->sid, *b);
}

/** Calculate least squares baseline solution from a set of single difference
 * observations and carrier phase ambiguities.
 *
 * \param num_sdiffs Number of single difference observations
 * \param sdiffs      Set of single differenced observations, length
 *                    `num_sdiffs`, sorted by PRN
 * \param ref_ecef    The reference position in ECEF frame, for computing the
 *                    sat direction vectors
 * \param num_ambs    Number of carrier phase ambiguities
 * \param single_ambs Array of (carrier phase ambiguity, prn) pairs.
 *                    length `num_ambs`
 * \param num_used    Pointer to where to store number of satellites used in the
 *                    baseline solution
 * \param b           The output baseline in meters
 * \param disable_raim   True disables raim check/repair
 * \param raim_threshold Threshold for raim checks.
 * \return            See baseline()
 */
s8 baseline_(u8 num_sdiffs, const sdiff_t *sdiffs, const double ref_ecef[3],
             u8 num_ambs, const ambiguity_t *single_ambs,
             u8 *num_used, double b[3],
             bool disable_raim, double raim_threshold)
{
  if (num_sdiffs < 4 || num_ambs < 4) {
    /* For a position solution, we need at least 4 sats. */
    return BASELINE_NOT_ENOUGH_SATS_ROVER;
  }

  assert(sdiffs != NULL);
  assert(ref_ecef != NULL);
  assert(single_ambs != NULL);
  assert(num_used != NULL);
  assert(b != NULL);

  assert(is_set(num_ambs, sizeof(ambiguity_t), single_ambs, cmp_amb));
  assert(is_set(num_sdiffs, sizeof(sdiff_t), sdiffs, cmp_sdiff));

  assert(num_sdiffs <= MAX_CHANNELS);

  /* Could use min(num_ambs, num_sdiffs) */
  ambiguity_t intersection_ambs[num_ambs];
  sdiff_t intersection_sdiffs[num_ambs];

  s32 intersection_size = intersection(
      num_ambs,   sizeof(ambiguity_t), single_ambs, intersection_ambs,
      num_sdiffs, sizeof(sdiff_t),     sdiffs,      intersection_sdiffs,
      cmp_amb_sdiff);

  if (intersection_size < 4) {
    /* For a position solution, we need at least 4 sats. */
    return BASELINE_NOT_ENOUGH_SATS_COMMON;
  }
  u8 num_dds = intersection_size - 1;

  /* Choose ref sat based on SNR. */
  gnss_signal_t ref_sid = choose_reference_sat(intersection_size, intersection_sdiffs);

  /* Calculate double differenced measurements. */
  sdiff_t sdiff_ref_first[intersection_size];
  u32 sdiff_ref_index = remove_element(intersection_size, sizeof(sdiff_t),
                                       intersection_sdiffs,
                                       &(sdiff_ref_first[1]),  /* New set */
                                       &ref_sid, cmp_sdiff_sid);
  memcpy(sdiff_ref_first, &intersection_sdiffs[sdiff_ref_index],
         sizeof(sdiff_t));

  double dd_meas[2 * num_dds];

  for (u32 i = 0; i < num_dds; i++) {
    dd_meas[i] =
      sdiff_ref_first[i+1].carrier_phase - sdiff_ref_first[0].carrier_phase;
    dd_meas[i + num_dds] =
      sdiff_ref_first[i+1].pseudorange - sdiff_ref_first[0].pseudorange;
  }

  double DE[num_dds * 3];
  assign_de_mtx(intersection_size, sdiff_ref_first, ref_ecef, DE);

  /* Calculate double differenced ambiguities. */
  double dd_ambs[num_dds];
  diff_ambs(ref_sid, intersection_size, intersection_ambs, dd_ambs);

  /* Compute least squares solution. */
  *num_used = intersection_size;
  return lesq_solve_raim(num_dds, dd_meas, dd_ambs, DE, b,
                         disable_raim, raim_threshold, 0, 0, 0);
}

void diff_ambs(gnss_signal_t ref_sid, u8 num_ambs, const ambiguity_t *amb_set,
               double *dd_ambs)
{
  u8 num_dds = num_ambs - 1;
  ambiguity_t amb_no_ref[num_dds];

  u32 amb_ref_index = remove_element(num_ambs, sizeof(ambiguity_t),
                                     amb_set,
                                     amb_no_ref,  /* New set */
                                     &ref_sid, cmp_amb_sid);
  for (u32 i = 0; i < num_dds; i++) {
    dd_ambs[i] = amb_no_ref[i].amb - amb_set[amb_ref_index].amb;
  }
}

/** Calculate least squares baseline solution from a set of single difference
 * observations and carrier phase ambiguities.
 *
 * \param num_sdiffs Number of single difference observations
 * \param sdiffs     Set of single differenced observations, length
 *                   `num_sdiffs`, sorted by PRN
 * \param ref_ecef   The reference position in ECEF frame, for computing the
 *                   sat direction vectors
 * \param ambs       Set of ambiguities as an `ambiguitites_t` structure
 * \param num_used   Pointer to where to store number of satellites used in the
 *                   baseline solution
 * \param b          The output baseline in meters
 * \param disable_raim True disables raim check/repair
 * \param raim_threshold Threshold for raim checks.
 *                       Value 5.5 has been extensively tested.
 * \return           BASELINE_SUCCESS on success,
 *                   BASELINE_NOT_ENOUGH_SATS_ROVER or
 *                   BASELINE_NOT_ENOUGH_SATS_COMMON
 *                     if there were insufficient observations to calculate the
 *                     baseline (the solution was under-constrained),
 *                   or one of the returns from lesq_solve_raim().
 */
s8 baseline(u8 num_sdiffs, const sdiff_t *sdiffs, const double ref_ecef[3],
            const ambiguities_t *ambs, u8 *num_used, double b[3],
            bool disable_raim, double raim_threshold)
{
  u8 num = ambs->n + 1;
  ambiguity_t ambts[ambs->n];
  ambiguity_t single_ambs[num];
  gnss_signal_t ref_sid = ambs->sids[0];
  ambiguity_t ref_amb = {.sid = ref_sid, .amb = 0};

  /* TODO(dsk) convert ambiguities_t to have a single-differenced ambiguity_t
   * array? */
  for (s32 i = 0; i < ambs->n; i ++) {
    /* ambs contains ambiguities for all non-ref prns */
    ambts[i].amb = ambs->ambs[i];
    /* prns contains ref prn followed by the rest */
    ambts[i].sid = ambs->sids[i+1];
  }

  insert_element(ambs->n, sizeof(ambiguity_t), ambts, single_ambs,
                 &ref_amb, cmp_amb);

  /* single_ambs is now an ambiguity_t set */
  return baseline_(num_sdiffs, sdiffs, ref_ecef,
                   num, single_ambs,
                   num_used, b,
                   disable_raim, raim_threshold);
}

/* Initialize a set of ambiguities.
 * Initializes to an empty set.
 *
 * \param ambs Pointer to set of ambiguities
 */
void ambiguities_init(ambiguities_t *ambs)
{
  ambs->n = 0;
}

/** \} */
