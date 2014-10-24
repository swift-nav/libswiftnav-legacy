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

/*
 * This is a Bierman-Thornton kalman filter implementation, as described in:
 *  [1] Gibbs, Bruce P. "Advanced Kalman Filtering, Least-Squares, and Modeling."
 *      John C. Wiley & Sons, Inc., 2011.
 */

#include <string.h>
#include <stdio.h>
#include <cblas.h>
#include <clapack.h>
// #include <lapacke.h>
#include <math.h>
#include <linear_algebra.h>
#include "constants.h"
#include "track.h"
#include "almanac.h"
#include "gpstime.h"
#include "amb_kf.h"

#define DEBUG_AMB_KF 0

/** \defgroup amb_kf Float Ambiguity Resolution
 * Preliminary integer ambiguity estimation with a Kalman Filter.
 * \{ */

/** In place updating of the state cov and k vec using a scalar observation
 * This is from section 10.2.1 of Gibbs [1], with some extra logic for handling
 *    singular matrices, dictating that zeros from cov_D dominate.
 *
 */
void incorporate_scalar_measurement(u32 state_dim, double *h, double R,
                               double *U, double *D, double *k)
{
  if (DEBUG_AMB_KF) {
    printf("<INCORPORATE_SCALAR_MEASUREMENT>\n");
    VEC_PRINTF(h, state_dim);
    printf("R = %.16f", R);
    if (abs(R) == 0) {
      printf(" \t (R == 0 exactly)\n");
    }
    else {
      printf("\n");
    }
  }

  double f[state_dim]; // f = U^T * h
  memcpy(f, h, state_dim * sizeof(double));
  cblas_dtrmv(CblasRowMajor, CblasUpper, CblasTrans, CblasUnit, //CBLAS_ORDER, CBLAS_UPLO, CBLAS_TRANSPOSE transA, CBLAS_DIAG
              state_dim, U, //int N, double *A
              state_dim, f, 1); // int lda, double *X, int incX


  double g[state_dim]; // g = diag(D) * f
  double alpha = R; // alpha = f * g + R = f^T * diag(D) * f + R
  for (u32 i=0; i<state_dim; i++) {
    g[i] = D[i] * f[i];
    alpha += f[i] * g[i];
  }
  if (DEBUG_AMB_KF) {
    VEC_PRINTF(f, state_dim);
    VEC_PRINTF(g, state_dim);
    printf("alpha = %.16f", alpha);
    if (abs(alpha) == 0) {
      printf(" \t (alpha == 0 exactly)\n");
    }
    else {
      printf("\n");
    }
  }

  double gamma[state_dim];
  double U_bar[state_dim * state_dim];
  double D_bar[state_dim];

  memset(gamma, 0,             state_dim * sizeof(double));
  memset(U_bar, 0, state_dim * state_dim * sizeof(double));
  memset(D_bar, 0,             state_dim * sizeof(double));
  memset(k,     0,             state_dim * sizeof(double));

  gamma[0] = R + g[0] * f[0];
  if (D[0] == 0 || R == 0) {
    D_bar[0] = 0; // this is just an expansion of the other branch with the proper 0 `div` 0 definitions
  }
  else {
    D_bar[0] = D[0] * R / gamma[0];
  }
  k[0] = g[0];
  U_bar[0] = 1;
  if (DEBUG_AMB_KF) {
    printf("gamma[0] = %f\n", gamma[0]);
    printf("D_bar[0] = %f\n", D_bar[0]);
    VEC_PRINTF(k, state_dim);
    printf("U_bar[:,0] = {");
    for (u32 i=0; i < state_dim; i++) {
      printf("%f, ", U_bar[i*state_dim]);
    }
    printf("}\n");
  }
  for (u32 j=1; j<state_dim; j++) {
    gamma[j] = gamma[j-1] + g[j] * f[j];
    if (D[j] == 0 || gamma[j-1] == 0) {
      D_bar[j] = 0; // this is just an expansion of the other branch with the proper 0 `div` 0 definitions
    }
    else {
      D_bar[j] = D[j] * gamma[j-1] / gamma[j];
    }
    double f_over_gamma = f[j] / gamma[j-1];
    for (u32 i=0; i<=j; i++) {
      if (k[i] == 0) {
        U_bar[i*state_dim + j] = U[i*state_dim + j]; // same as other branch with correct 0 `div` 0 definitions
      }
      else {
        U_bar[i*state_dim + j] = U[i*state_dim + j] - f_over_gamma * k[i]; // U_bar[:,j] = U[:,j] - f[j]/gamma[j-1] * k
      }
      k[i] += g[j] * U[i*state_dim + j]; // k = k + g[j] * U[:,j]
    }
    if (DEBUG_AMB_KF) {
      printf("gamma[%"PRIu32"] = %f\n", j, gamma[j]);
      printf("D_bar[%"PRIu32"] = %f\n", j, D_bar[j]);
      VEC_PRINTF(k, state_dim);
      printf("U_bar[:,%"PRIu32"] = {", j);
      for (u32 i=0; i < state_dim; i++) {
        printf("%f, ", U_bar[i*state_dim + j]);
      }
      printf("}\n");
    }
  }
  for (u32 i=0; i<state_dim; i++) {
    k[i] /= alpha;
  }
  memcpy(U, U_bar, state_dim * state_dim * sizeof(double));
  memcpy(D, D_bar,             state_dim * sizeof(double));
  if (DEBUG_AMB_KF) {
    MAT_PRINTF(U, state_dim, state_dim);
    VEC_PRINTF(D, state_dim);
    printf("</INCORPORATE_SCALAR_MEASUREMENT>\n");
  }
}

/** In place updating of the state mean and covariances to use the (decorrelated) observations
 * This is directly from section 10.2.1 of Gibbs [1]
 */
void incorporate_obs(nkf_t *kf, double *decor_obs)
{
  if (DEBUG_AMB_KF) {
    printf("<INCORPORATE_OBS>\n");
  }
  for (u32 i=0; i<kf->obs_dim; i++) {
    double *h = &kf->decor_obs_mtx[kf->state_dim * i]; //vector of length kf->state_dim
    double R = kf->decor_obs_cov[i]; //scalar
    double k[kf->state_dim]; // vector of length kf->state_dim
    // printf("i=%i\n", i);
    // VEC_PRINTF(h, kf->state_dim);

    incorporate_scalar_measurement(kf->state_dim, h, R, kf->state_cov_U, kf->state_cov_D, &k[0]); //updates cov and sets k
    // VEC_PRINTF(k, kf->state_dim);

    double predicted_obs = 0;
    for (u32 j=0; j<kf->state_dim; j++) {//TODO take advantage of sparsity of h
      predicted_obs += h[j] * kf->state_mean[j];
    }
    double obs_minus_predicted_obs = decor_obs[i] - predicted_obs;
    // printf("decor_obs = %f\n", decor_obs[i]);
    // printf("predi_obs = %f\n", predicted_obs);
    // printf(" diff_obs = %f\n", obs_minus_predicted_obs);

    for (u32 j=0; j<kf->state_dim; j++) {
      kf->state_mean[j] += k[j] * obs_minus_predicted_obs; // uses k to update mean
    }
    // VEC_PRINTF(intermediate_mean, kf->state_dim);
  }
  if (DEBUG_AMB_KF) {
    printf("</INCORPORATE_OBS>\n");
  }
}

// turns (phi, rho) into Q_tilde * (phi, rho)
void make_residual_measurements(nkf_t *kf, double *measurements, double *resid_measurements)
{
  u8 constraint_dim = MAX(kf->state_dim, 3) - 3;
  cblas_dgemv (CblasRowMajor, CblasNoTrans, //Order, TransA
               constraint_dim, kf->state_dim, // M, N
               1, kf->null_basis_Q, kf->state_dim, // alpha A, lda
               measurements, 1, // X, incX
               0, resid_measurements, 1); // beta, Y, incY
  for (u8 i=0; i< kf->state_dim; i++) {
    resid_measurements[i+constraint_dim] = measurements[i] - measurements[i+kf->state_dim] / GPS_L1_LAMBDA_NO_VAC;
  }
}

void diffuse_state(nkf_t *kf)
{
  for (u8 i=0; i< kf->state_dim; i++) {
    kf->state_cov_D[i] += kf->amb_drift_var; //TODO make this a tunable parameter defined at the right time
  }
}


/** In place updating of the state mean and covariance. Modifies measurements.
 */
void nkf_update(nkf_t *kf, double *measurements)
{
  if (DEBUG_AMB_KF) {
    printf("<NKF_UPDATE>\n");
  }
  double resid_measurements[kf->obs_dim];
  make_residual_measurements(kf, measurements, resid_measurements);
  // VEC_PRINTF(measurements, kf->obs_dim);
  // MAT_PRINTF(kf->decor_obs_mtx, kf->obs_dim, kf->state_dim);
  // MAT_PRINTF(kf->obs_cov_root_inv, kf->obs_dim, kf->obs_dim);

  // replaces residual measurements by their decorrelated version
  cblas_dtrmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasUnit, // Order, Uplo, TransA, Diag
              kf->obs_dim, kf->decor_mtx, kf->obs_dim, // N, A, lda
              resid_measurements, 1); // X, incX

  // predict_forward(kf);
  diffuse_state(kf);
  incorporate_obs(kf, resid_measurements);
  if (DEBUG_AMB_KF) {
    MAT_PRINTF(kf->state_cov_U, kf->state_dim, kf->state_dim);
    VEC_PRINTF(kf->state_cov_D, kf->state_dim);
    VEC_PRINTF(kf->state_mean, kf->state_dim);
    printf("</NKF_UPDATE>\n");
  }
}

// presumes that the first alm entry is the reference sat
void assign_de_mtx(u8 num_sats, sdiff_t *sats_with_ref_first, double ref_ecef[3], double *DE)
{
  if (DEBUG_AMB_KF) {
    printf("<ASSIGN_DE_MTX>\nnum_sats = %u\nsdiff prns&positions = {\n", num_sats);
    for (u8 i=0; i < num_sats; i++) {
      printf("i = %u, prn = %u, \tpos = [%f, \t%f, \t%f]\n",
             i,
             sats_with_ref_first[i].prn,
             sats_with_ref_first[i].sat_pos[0],
             sats_with_ref_first[i].sat_pos[1],
             sats_with_ref_first[i].sat_pos[2]);
    }
    printf("}\nref_ecef = {%f, \t%f, \t%f}\n", ref_ecef[0], ref_ecef[1], ref_ecef[2]);
  }
  if (num_sats <= 1) {
    if (DEBUG_AMB_KF) {
      printf("not enough sats\n</ASSIGN_DE_MTX>\n");
    }
    return;
  }

  memset(DE, 0, (num_sats - 1) * 3 * sizeof(double));
  double e0[3];
  double x0 = sats_with_ref_first[0].sat_pos[0] - ref_ecef[0];
  double y0 = sats_with_ref_first[0].sat_pos[1] - ref_ecef[1];
  double z0 = sats_with_ref_first[0].sat_pos[2] - ref_ecef[2];
  double norm0 = sqrt(x0*x0 + y0*y0 + z0*z0);
  e0[0] = x0 / norm0;
  e0[1] = y0 / norm0;
  e0[2] = z0 / norm0;
  // VEC_PRINTF(ref_ecef, 3);
  for (u8 i=1; i<num_sats; i++) {
    double x = sats_with_ref_first[i].sat_pos[0] - ref_ecef[0];
    double y = sats_with_ref_first[i].sat_pos[1] - ref_ecef[1];
    double z = sats_with_ref_first[i].sat_pos[2] - ref_ecef[2];
    double norm = sqrt(x*x + y*y + z*z);
    DE[3*(i-1)] = x / norm - e0[0];
    DE[3*(i-1) + 1] = y / norm - e0[1];
    DE[3*(i-1) + 2] = z / norm - e0[2];
  }
  if (DEBUG_AMB_KF) {
    MAT_PRINTF(DE, (num_sats-1), 3);
    printf("</ASSIGN_DE_MTX>\n");
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
/** A least squares solution for baseline from phases using the KF state.
 * This uses the current state of the KF and a set of phase observations to
 * solve for the current baseline.
 *
 * \param kf                    The Kalman filter struct.
 * \param sdiffs_with_ref_first A list of sdiffs. The first in the list must be
 *                              the reference sat of the KF, and the rest must
 *                              correspond to the KF's DD amb estimates' sats.
 * \param dd_measurements       A vector of carrier phases. They must be double
 *                              differenced and ordered according to the sdiffs
 *                              and KF's sats (which must match each other).
 * \param ref_ecef              The reference position in ECEF frame, for
 *                              computing the sat direction vectors.
 * \param b                     The output baseline in meters.
 */
void least_squares_solve_b(nkf_t *kf, sdiff_t *sdiffs_with_ref_first,
                      double *dd_measurements, double ref_ecef[3], double b[3])
{
  if (DEBUG_AMB_KF) {
    printf("<LEAST_SQUARES_SOLVE_B>\n");
  }
  integer num_dds = kf->state_dim;
  double DE[num_dds * 3];
  assign_de_mtx(num_dds+1, sdiffs_with_ref_first, ref_ecef, DE);
  double DET[num_dds * 3];
  for (u8 i=0; i<num_dds; i++) { //TODO this transposition is stupid and unnecessary
    DET[              i] = DE[i*3 + 0];
    DET[    num_dds + i] = DE[i*3 + 1];
    DET[2 * num_dds + i] = DE[i*3 + 2];
  }

  double phase_ranges[MAX(num_dds,3)];
  for (u8 i=0; i< num_dds; i++) {
    phase_ranges[i] = dd_measurements[i] - kf->state_mean[i];
  }

  if (DEBUG_AMB_KF) {
    printf("\tdd_measurements, \tkf->state_mean, \tdifferenced phase_ranges = {\n");
    for (u8 i=0; i< num_dds; i++) {
      printf("\t%f, \t%f, \t%f,\n", dd_measurements[i], kf->state_mean[i], phase_ranges[i]);
    }
    printf("\t}\n");
  }

  /* TODO could use plain old DGELS here */

  s32 ldb = (s32) MAX(num_dds,3);
  double s[3];
  double rcond = 1e-12;
  s32 rank;
  double w[1]; //try 25 + 10*num_sats
  s32 lwork = -1;
  s32 info;
  s32 three = 3;
  s32 one = 1;
  dgelss_(&num_dds, &three, &one, //M, N, NRHS
          DET, &num_dds, //A, LDA
          phase_ranges, &ldb, //B, LDB
          s, &rcond, // S, RCOND
          &rank, //RANK
          w, &lwork, // WORK, LWORK
          &info); //INFO
  lwork = round(w[0]);

  double work[lwork];
  dgelss_(&num_dds, &three, &one, //M, N, NRHS
          DET, &num_dds, //A, LDA
          phase_ranges, &ldb, //B, LDB
          s, &rcond, // S, RCOND
          &rank, //RANK
          work, &lwork, // WORK, LWORK
          &info); //INFO
  b[0] = phase_ranges[0] * GPS_L1_LAMBDA_NO_VAC;
  b[1] = phase_ranges[1] * GPS_L1_LAMBDA_NO_VAC;
  b[2] = phase_ranges[2] * GPS_L1_LAMBDA_NO_VAC;
  if (DEBUG_AMB_KF) {
    printf("b = {%f, %f, %f}\n", b[0]*100, b[1]*100, b[2]*100); // units --> cm
    printf("</LEAST_SQUARES_SOLVE_B>\n");
  }
}


/* presumes that the first alm entry is the reference sat.
 *TODO use the state covariance matrix for a better estimate:
 *  That is, decorrelate and scale the LHS of y = A * x before solving for x
 *TODO consolidate this with the one using the float solution
 */
void least_squares_solve_b_external_ambs(u8 num_dds_u8, double *ambs, sdiff_t *sdiffs_with_ref_first, double *dd_measurements, double ref_ecef[3], double b[3])
{
  integer num_dds = num_dds_u8;
  double DE[num_dds * 3];
  assign_de_mtx(num_dds+1, sdiffs_with_ref_first, ref_ecef, DE);
  double DET[num_dds * 3];
  for (u8 i=0; i<num_dds; i++) { //TODO this transposition is stupid and unnecessary
    DET[              i] = DE[i*3 + 0];
    DET[    num_dds + i] = DE[i*3 + 1];
    DET[2 * num_dds + i] = DE[i*3 + 2];
  }

  // double fake_ints[num_dds];
  // for (u8 i=0; i< num_dds; i++) {
  //   fake_ints[i] = (i+1)*10; //TODO REMOVE THIS version
  // }
  // double resid[num_dds];
  // lesq_solution_b(kf->state_dim, dd_measurements, kf->state_mean, DE, b, resid);
  // lesq_solution_b(kf->state_dim, dd_measurements, fake_ints, DE, b, resid);

  double phase_ranges[MAX(num_dds,3)];
  for (u8 i=0; i< num_dds; i++) {
    phase_ranges[i] = dd_measurements[i] - ambs[i];
    // phase_ranges[i] = dd_measurements[i] - (i+1)*10;
  }

  //TODO could use plain old DGELS here

  s32 ldb = (s32) MAX(num_dds,3);
  double s[3];
  double rcond = 1e-12;
  s32 rank;
  double w[1]; //try 25 + 10*num_sats
  s32 lwork = -1;
  s32 info;
  s32 three = 3;
  s32 one = 1;
  dgelss_(&num_dds, &three, &one, //M, N, NRHS
          DET, &num_dds, //A, LDA
          phase_ranges, &ldb, //B, LDB
          s, &rcond, // S, RCOND
          &rank, //RANK
          w, &lwork, // WORK, LWORK
          &info); //INFO
  lwork = round(w[0]);

  double work[lwork];
  dgelss_(&num_dds, &three, &one, //M, N, NRHS
          DET, &num_dds, //A, LDA
          phase_ranges, &ldb, //B, LDB
          s, &rcond, // S, RCOND
          &rank, //RANK
          work, &lwork, // WORK, LWORK
          &info); //INFO

  memcpy(b, phase_ranges, 3 * sizeof(double));
  b[0] = phase_ranges[0] * GPS_L1_LAMBDA_NO_VAC;
  b[1] = phase_ranges[1] * GPS_L1_LAMBDA_NO_VAC;
  b[2] = phase_ranges[2] * GPS_L1_LAMBDA_NO_VAC;
}


// initializes the ambiguity means and variances.
// Note that the covariance is  in UDU form, and U starts as identity.
void initialize_state(nkf_t *kf, double *dd_measurements, double init_var)
{
  u8 num_dds = kf->state_dim;
  for (u32 i=0; i<num_dds; i++) {
    // N = Expectation of [phi - rho / lambda]
    kf->state_mean[i] = dd_measurements[i] - dd_measurements[i + num_dds] / GPS_L1_LAMBDA_NO_VAC;
    // Sigma begins as a diagonal
    kf->state_cov_D[i] = init_var;
  }
  matrix_eye(num_dds, kf->state_cov_U);
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
  for (u8 i=0; i < num_dds; i++) {
    for (u8 j=0; j<3; j++) {
      A[j*num_dds + i] = DE_mtx[i*3 + j]; //set A = Transpose(DE_mtx)
    }
  }
  integer m = num_dds;
  integer n = 3;
  double tau[3];
  QR_part1(m, n, A, tau);
  QR_part2(m, n, A, tau);
  memcpy(q, &A[3*num_dds], (MAX(3, num_dds) - 3) * num_dds * sizeof(double));
}

void assign_dd_obs_cov(u8 num_dds, double phase_var, double code_var, double *dd_obs_cov) //TODO this could be made more efficient, if it matters
{
  for (u8 i = 0; i < num_dds; i++) {
    for (u8 j = 0; j < num_dds; j++) {
      if (i == j) {
        dd_obs_cov[i*2*num_dds + j] = 2 * phase_var;
        dd_obs_cov[(i+num_dds)*2*num_dds + num_dds + j] = 2 * code_var;
      }
      else {
        dd_obs_cov[i*2*num_dds + j] = phase_var;
        dd_obs_cov[(i+num_dds)*2*num_dds + num_dds + j] = code_var;
      }
    }
    memset(&dd_obs_cov[i*2*num_dds + num_dds], 0, num_dds * sizeof(double));
    memset(&dd_obs_cov[(i+num_dds)*2*num_dds], 0, num_dds * sizeof(double));
  }
}

void assign_residual_obs_cov(u8 num_dds, double phase_var, double code_var, double *q, double *r_cov) //TODO make this more efficient (e.g. via pages 3/6.2-3/2014 of ian's notebook)
{
  double dd_obs_cov[4 * num_dds * num_dds];
  assign_dd_obs_cov(num_dds, phase_var, code_var, dd_obs_cov);
  integer nullspace_dim = MAX(num_dds, 3) - 3;;
  integer dd_dim = 2*num_dds;
  integer res_dim = num_dds + nullspace_dim;
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
              1, dd_obs_cov, dd_dim, // double alpha, double *A, int lda
              q_tilde, dd_dim, // double *B, int ldb
              0, QC, dd_dim); // double beta, double *C, int ldc
  // MAT_PRINTF(QC, res_dim, dd_dim);

  //TODO make more efficient via the structure of q_tilde, and it's relation to the I + 1*1^T structure of the obs cov mtx
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, // CBLAS_ORDER, CBLAS_TRANSPOSE transA, cBLAS_TRANSPOSE transB
              res_dim, res_dim, dd_dim, // int M, int N, int K
              1, QC, dd_dim, // double alpha, double *A, int lda
              q_tilde, dd_dim, //double *B, int ldb
              0, r_cov, res_dim); //beta, double *C, int ldc
  // MAT_PRINTF(r_cov_inv, res_dim, res_dim);
}

void invert_U(u8 res_dim, double *U) // in place inversion of U
{
  char uplo = 'U'; //upper triangular
  char diag = 'U'; //unit triangular
  integer dim = res_dim;
  integer lda = res_dim;
  integer info;
  dtrtri_(&uplo, &diag,
          &dim, U, &lda, &info);
}

void assign_simple_sig(u8 num_dds, double var, double *simple_cov)
{
  for (u8 i = 0; i < num_dds; i++) {
    for (u8 j = 0; j < num_dds; j++) {
      if (i == j) {
        simple_cov[i*num_dds + j] = 2 * var;
      }
      else {
        simple_cov[i*num_dds + j] = var;
      }
    }
  }
}

// from get_kf_matrices:
// U^-1 * y == y'
//           = U^-1 * H * x
//          == H' * x
//
// H = ( Q )
//     ( I )
// where Q's rows form a basis for the left null space for DE
void assign_H_prime(u8 res_dim, u8 constraint_dim, u8 num_dds,
                    double *Q, double *U_inv, double *H_prime)
{
  // set the H_prime variable to equal H
  memcpy(H_prime, Q, constraint_dim * num_dds * sizeof(double));
  matrix_eye(num_dds, &H_prime[constraint_dim * num_dds]);

  // multiply H_prime by U_inv to make it the actual H_prime
  cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasUnit, // CBLAS_ORDER, CBLAS_SIDE, CBLAS_UPLO, CBLAS_TRANSPOSE, CBLAS_DIAG
              res_dim, num_dds, //M, N
              1, U_inv, res_dim, //alpha, A, lda
              H_prime, num_dds); //B, ldb
}

// y = H * x
// Var[y] = Sig = U * D * U^T
// ==> Var[U^-1 * y] = D
// U^-1 * y == y'
//           = U^-1 * H * x
//          == H' * x
//
// H = ( Q )
//     ( I )
// where Q's rows form a basis for the left null space for DE
//
// Sig = Q~ * Sig_v * Q~^T
//   where    Q~ = ( Q   0        )
//                 ( I  -I/lambda )
//     and Sig_v = ( D*D^T * var_phi   0               )
//                 ( 0                 D*D^T * var_rho )
//     and D*D^T = 1*1^T + I
//
// This function constructs D, U^-1, and H'
void get_kf_matrices(u8 num_sdiffs, sdiff_t *sdiffs_with_ref_first,
                     double ref_ecef[3],
                     double phase_var, double code_var,
                     double *null_basis_Q,
                     double *U_inv, double *D,
                     double *H_prime)
{
  u8 num_dds = MAX(1, num_sdiffs) - 1;
  u8 constraint_dim = MAX(num_dds, 3) - 3;;
  u8 res_dim = num_dds + constraint_dim;

  double Sig[res_dim * res_dim];

  //assign Sig and H
  if (constraint_dim > 0) {
    double DE[num_dds * 3];
    assign_de_mtx(num_sdiffs, sdiffs_with_ref_first, ref_ecef, DE);
    assign_phase_obs_null_basis(num_dds, DE, null_basis_Q);
    assign_residual_obs_cov(num_dds, phase_var, code_var, null_basis_Q, Sig);
    //TODO U seems to have that fancy blockwise structure we love so much. Use it
    matrix_udu(res_dim, Sig, U_inv, D); //U_inv holds U after this
    invert_U(res_dim, U_inv);
    //TODO this also has fancy structure
    assign_H_prime(res_dim, constraint_dim, num_dds, null_basis_Q, U_inv, H_prime);
  }
  else {
    assign_simple_sig(num_dds,
                      phase_var + code_var / (GPS_L1_LAMBDA_NO_VAC * GPS_L1_LAMBDA_NO_VAC),
                      Sig);
    matrix_udu(res_dim, Sig, U_inv, D); //U_inv holds U after this
    invert_U(res_dim, U_inv);

    // H = I in this case, so H' = U^-1 * H = U^-1
    memcpy(H_prime, U_inv, num_dds * num_dds * sizeof(double));
  }

}


void set_nkf(nkf_t *kf, double amb_drift_var, double phase_var, double code_var, double amb_init_var,
            u8 num_sdiffs, sdiff_t *sdiffs_with_ref_first, double *dd_measurements, double ref_ecef[3])
{
  if (DEBUG_AMB_KF) {
      printf("<SET_NKF>\n");
  }
  u32 state_dim = num_sdiffs - 1;
  u32 num_diffs = num_sdiffs - 1;
  kf->state_dim = state_dim;
  u32 constraint_dim = MAX(3, num_sdiffs) - 3;
  kf->obs_dim = num_diffs + constraint_dim;
  kf->amb_drift_var = amb_drift_var;

  get_kf_matrices(num_sdiffs, sdiffs_with_ref_first,
                  ref_ecef,
                  phase_var, code_var,
                  kf->null_basis_Q,
                  kf->decor_mtx, kf->decor_obs_cov,
                  kf->decor_obs_mtx);

  // given plain old measurements, initialize the state
  initialize_state(kf, dd_measurements, amb_init_var);
  if (DEBUG_AMB_KF) {
      printf("</SET_NKF>\n");
  }
}

void set_nkf_matrices(nkf_t *kf, double phase_var, double code_var,
                     u8 num_sdiffs, sdiff_t *sdiffs_with_ref_first, double ref_ecef[3])
{
  u32 state_dim = MAX(1, num_sdiffs) - 1;
  u32 num_diffs = MAX(1, num_sdiffs) - 1;
  kf->state_dim = state_dim;
  u32 constraint_dim = MAX(3, num_diffs) - 3;
  kf->obs_dim = num_diffs + constraint_dim;

  get_kf_matrices(num_sdiffs, sdiffs_with_ref_first,
                  ref_ecef,
                  phase_var, code_var,
                  kf->null_basis_Q,
                  kf->decor_mtx, kf->decor_obs_cov,
                  kf->decor_obs_mtx);
}

s32 find_index_of_element_in_u8s(u32 num_elements, u8 x, u8 *list)
{
  for (u32 i=0; i<num_elements; i++) {
    if (x == list[i]) {
      return i;
    }
  }
  return -1;
}

void rebase_mean_N(double *mean, u8 num_sats, u8 *old_prns, u8 *new_prns)
{
  u8 old_ref = old_prns[0];
  u8 new_ref = new_prns[0];

  double new_mean[num_sats-1];
  s32 index_of_new_ref_in_old = find_index_of_element_in_u8s(num_sats, new_ref, &old_prns[1]);
  double val_for_new_ref_in_old_basis = mean[index_of_new_ref_in_old];
  for (u8 i=0; i<num_sats-1; i++) {
    u8 new_prn = new_prns[1+i];
    if (new_prn == old_ref) {
      new_mean[i] = - val_for_new_ref_in_old_basis;
    }
    else {
      s32 index_of_this_sat_in_old_basis = find_index_of_element_in_u8s(num_sats, new_prn, &old_prns[1]);
      new_mean[i] = mean[index_of_this_sat_in_old_basis] - val_for_new_ref_in_old_basis;
    }
  }
  memcpy(mean, new_mean, (num_sats-1) * sizeof(double));
}

void assign_state_rebase_mtx(u8 num_sats, u8 *old_prns, u8 *new_prns, double *rebase_mtx)
{
  u8 state_dim = num_sats - 1;
  memset(rebase_mtx, 0, state_dim * state_dim * sizeof(double));
  u8 old_ref = old_prns[0];
  u8 new_ref = new_prns[0];

  s32 index_of_new_ref_in_old = find_index_of_element_in_u8s(num_sats-1, new_ref, &old_prns[1]);
  s32 index_of_old_ref_in_new = find_index_of_element_in_u8s(num_sats-1, old_ref, &new_prns[1]);
  for (u8 i=0; i<state_dim; i++) {
    rebase_mtx[i*state_dim + index_of_new_ref_in_old] = -1;
    if (i != (u8) index_of_old_ref_in_new) {
      s32 index_of_this_sat_in_old_basis = find_index_of_element_in_u8s(num_sats-1, new_prns[i+1], &old_prns[1]);
      rebase_mtx[i*state_dim + index_of_this_sat_in_old_basis] = 1;
    }
  }
  // MAT_PRINTF(rebase_mtx, state_dim, state_dim);
}


void rebase_covariance_sigma(double *state_cov, u8 num_sats, u8 *old_prns, u8 *new_prns)
{
  u8 state_dim = num_sats - 1;
  double rebase_mtx[state_dim * state_dim];
  assign_state_rebase_mtx(num_sats, old_prns, new_prns, rebase_mtx);

  double intermediate_cov[state_dim * state_dim];
  //TODO make more efficient via structure of rebase_mtx
  cblas_dsymm(CblasRowMajor, CblasRight, CblasUpper, //CBLAS_ORDER, CBLAS_SIDE, CBLAS_UPLO
              state_dim, state_dim, // int M, int N
              1, state_cov, state_dim, // double alpha, double *A, int lda
              rebase_mtx, state_dim, // double *B, int ldb
              0, intermediate_cov, state_dim); // double beta, double *C, int ldc
  // MAT_PRINTF(intermediate_cov, state_dim, state_dim);

  //TODO make more efficient via the structure of rebase_mtx
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, // CBLAS_ORDER, CBLAS_TRANSPOSE transA, cBLAS_TRANSPOSE transB
              state_dim, state_dim, state_dim, // int M, int N, int K
              1, intermediate_cov, state_dim, // double alpha, double *A, int lda
              rebase_mtx, state_dim, //double *B, int ldb
              0, state_cov, state_dim); //beta, double *C, int ldc
  // MAT_PRINTF(state_cov, state_dim, state_dim);
}

void rebase_covariance_udu(double *state_cov_U, double *state_cov_D, u8 num_sats, u8 *old_prns, u8 *new_prns)
{
  u8 state_dim = num_sats - 1;
  double state_cov[state_dim * state_dim];
  matrix_reconstruct_udu(state_dim, state_cov_U, state_cov_D, state_cov);
  rebase_covariance_sigma(state_cov, num_sats, old_prns, new_prns);
  matrix_udu(state_dim, state_cov, state_cov_U, state_cov_D);
}


void rebase_nkf(nkf_t *kf, u8 num_sats, u8 *old_prns, u8 *new_prns)
{
  rebase_mean_N(kf->state_mean, num_sats, old_prns, new_prns);
  rebase_covariance_udu(kf->state_cov_U, kf->state_cov_D, num_sats, old_prns, new_prns);
}

void nkf_state_projection(nkf_t *kf,
                                    u8 num_old_non_ref_sats,
                                    u8 num_new_non_ref_sats,
                                    u8 *ndx_of_new_sat_in_old)
{
  u8 old_state_dim = num_old_non_ref_sats;
  double old_cov[old_state_dim * old_state_dim];
  matrix_reconstruct_udu(old_state_dim, kf->state_cov_U, kf->state_cov_D, old_cov);
  /*MAT_PRINTF(old_cov, old_state_dim, old_state_dim);*/

  u8 new_state_dim = num_new_non_ref_sats;
  double new_cov[new_state_dim * new_state_dim];
  double new_mean[new_state_dim];

  for (u8 i=0; i<num_new_non_ref_sats; i++) {
    u8 ndxi = ndx_of_new_sat_in_old[i];
    new_mean[i] = kf->state_mean[ndxi];
    for (u8 j=0; j<num_new_non_ref_sats; j++) {
      u8 ndxj = ndx_of_new_sat_in_old[j];
      new_cov[i*new_state_dim + j] = old_cov[ndxi*old_state_dim + ndxj];
    }
  }
  /*MAT_PRINTF(new_cov, new_state_dim, new_state_dim);*/

  //put it all back into the kf
  memcpy(kf->state_mean, new_mean, new_state_dim * sizeof(double));
  matrix_udu(new_state_dim, new_cov, kf->state_cov_U, kf->state_cov_D);
  //NOTE: IT DOESN'T UPDATE THE OBSERVATION OR TRANSITION MATRICES, JUST THE STATE
}

void nkf_state_inclusion(nkf_t *kf,
                                   u8 num_old_non_ref_sats,
                                   u8 num_new_non_ref_sats,
                                   u8 *ndx_of_old_sat_in_new,
                                   double int_init_var)
{
  u8 old_state_dim = num_old_non_ref_sats;
  double old_cov[old_state_dim * old_state_dim];
  matrix_reconstruct_udu(old_state_dim, kf->state_cov_U, kf->state_cov_D, old_cov);

  u8 new_state_dim = num_new_non_ref_sats;
  double new_cov[new_state_dim * new_state_dim];
  memset(new_cov, 0, new_state_dim * new_state_dim * sizeof(double));
  double new_mean[new_state_dim];
  // initialize_state(kf, dd_measurements, amb_init_var); //TODO do we really want to initialize new states this way?
  memset(new_mean, 0, new_state_dim * sizeof(double));
  for (u8 i=0; i<num_new_non_ref_sats; i++) {
    new_cov[i*new_state_dim + i] = int_init_var;
  }

  for (u8 i=0; i<num_old_non_ref_sats; i++) {
    u8 ndxi = ndx_of_old_sat_in_new[i];
    new_mean[ndxi] = kf->state_mean[i];
    for (u8 j=0;j<num_old_non_ref_sats; j++) {
      u8 ndxj = ndx_of_old_sat_in_new[j];
      new_cov[ndxi*new_state_dim + ndxj] = old_cov[i*old_state_dim + j];
    }
  }
  matrix_udu(new_state_dim, new_cov, kf->state_cov_U, kf->state_cov_D);
  memcpy(kf->state_mean, new_mean, new_state_dim * sizeof(double));
}

/** \} */
