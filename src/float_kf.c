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
#include <cblas.h>
#include <clapack.h>
#include <math.h>
#include <linear_algebra.h>
#include "constants.h"
#include "track.h"
#include "almanac.h"
#include "gpstime.h"
#include "amb_kf.h"
#include "float_kf.h"

s8 udu(u32 n, double *M, double *U, double *D) //todo: replace with DSYTRF
{
  double alpha, beta;
  triu(n, M);
  eye(n, U);
  memset(D, 0, n * sizeof(double));

  for (u32 j=n; j>=2; j--) {
    D[j - 1] = M[(j-1)*n + j-1];
    if (D[j-1] > 0) {
      alpha = 1.0 / D[j-1];
    } else {
      alpha = 0.0;
    }
    for (u32 k=1; k<j; k++) {
      beta = M[(k-1)*n + j-1];
      U[(k-1)*n + j-1] = alpha * beta;
      for (u32 kk = 0; kk < k; kk++) {
        M[kk*n + k-1] = M[kk*n + k-1] - beta * U[kk*n + j-1];
      }
    }

  }
  D[0] = M[0];
  return 0;
}

void eye(u32 n, double *M)
{
  memset(M, 0, n * n * sizeof(double));
  for (u32 i=0; i<n; i++) {
    M[i*n + i] = 1;
  }
}

void triu(u32 n, double *M)
{
  for (u32 i=1; i<n; i++) {
    for (u32 j=0; j<i; j++) {
      M[i*n + j] = 0;
    }
  }
}

/** In place prediction of the next state mean and covariances
 */
void predict_forward(kf_t *kf)
{
  //TODO take advantage of sparsity in this function

  double x[kf->state_dim];
  memcpy(x, kf->state_mean, kf->state_dim * sizeof(double));

  //TODO make more efficient via the structure of the transition matrix
  cblas_dgemv(CblasRowMajor, CblasNoTrans, // CBLAS_ORDER, CBLAS_TRANSPOSE
              kf->state_dim, kf->state_dim, // int M, int N,
              1, (double *) kf->transition_mtx, kf->state_dim, // double 1, double *A, int lda
              x, 1, // double *X, int incX
              0, kf->state_mean, 1); // double beta, double *Y, int incY
  // VEC_PRINTF((double *) state_mean, kf->state_dim);

  double state_cov[kf->state_dim * kf->state_dim];
  reconstruct_udu(kf->state_dim, kf->state_cov_U, kf->state_cov_D, state_cov);
  // MAT_PRINTF((double *) state_cov, kf->state_dim, kf->state_dim);

  //TODO make more efficient via the structure of the transition matrix
  double FC[kf->state_dim * kf->state_dim];
  cblas_dsymm(CblasRowMajor, CblasRight, CblasUpper, //CBLAS_ORDER, CBLAS_SIDE, CBLAS_UPLO
              kf->state_dim, kf->state_dim, // int M, int N
              1, state_cov, kf->state_dim, // double alpha, double *A, int lda
              kf->transition_mtx, kf->state_dim, // double *B, int ldb
              0, FC, kf->state_dim); // double beta, double *C, int ldc
  // MAT_PRINTF((double *) FC, kf->state_dim, kf->state_dim);

  //TODO make more efficient via the structure of the transition matrix
  double FCF[kf->state_dim * kf->state_dim];
  memcpy(FCF, kf->transition_cov, kf->state_dim * kf->state_dim * sizeof(double));
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, // CBLAS_ORDER, CBLAS_TRANSPOSE transA, cBLAS_TRANSPOSE transB
              kf->state_dim, kf->state_dim, kf->state_dim, // int M, int N, int K
              1, FC, kf->state_dim, // double alpha, double *A, int lda
              kf->transition_mtx, kf->state_dim, //double *B, int ldb
              1, FCF, kf->state_dim); //beta, double *C, int ldc
  // MAT_PRINTF((double *) FCF, kf->state_dim, kf->state_dim);

  udu(kf->state_dim, FCF, kf->state_cov_U, kf->state_cov_D);
  // MAT_PRINTF((double *) state_cov_U, kf->state_dim, kf->state_dim);
  // VEC_PRINTF((double *) state_cov_D, kf->state_dim);
}

/** In place updating of the state mean and covariances to use the (decorrelated) observations
 */
static void incorporate_obs(kf_t *kf, double *decor_obs)
{

  for (u32 i=0; i<kf->obs_dim; i++) {
    double *h = &kf->decor_obs_mtx[kf->state_dim * i]; //vector of length kf->state_dim
    double R = kf->decor_obs_cov[i]; //scalar
    double k[kf->state_dim]; // vector of length kf->state_dim
    // printf("i=%i\n", i);
    // VEC_PRINTF(h, kf->state_dim);

    update_scalar_measurement(kf->state_dim, h, R, kf->state_cov_U, kf->state_cov_D, &k[0]); //updates cov and sets k
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
}

/** In place updating of the state covariances and k vec to use a single decorrelated observation
 */
void update_scalar_measurement(u32 state_dim, double *h, double R,
                               double *U, double *D, double *k)
{
  // VEC_PRINTF(h, state_dim);
  double f[state_dim]; // f = U^T * h
  memcpy(f, h, state_dim * sizeof(double));
  cblas_dtrmv(CblasRowMajor, CblasUpper, CblasTrans, CblasUnit, //CBLAS_ORDER, CBLAS_UPLO, CBLAS_TRANSPOSE transA, CBLAS_DIAG
              state_dim, U, //int N, double *A
              state_dim, f, 1); // int lda, double *X, int incX
  // VEC_PRINTF(f, state_dim);

  double g[state_dim]; // g = diag(D) * f
  double alpha = R; // alpha = f * g + R = f^T * diag(D) * f + R
  for (u32 i=0; i<state_dim; i++) {
    g[i] = D[i] * f[i];
    alpha += f[i] * g[i];
  }
  // VEC_PRINTF(g, state_dim);
  // printf("%f\n", alpha);

  double gamma[state_dim];
  double U_bar[state_dim * state_dim];
  double D_bar[state_dim];

  memset(gamma, 0,             state_dim * sizeof(double));
  memset(U_bar, 0, state_dim * state_dim * sizeof(double));
  memset(D_bar, 0,             state_dim * sizeof(double));
  memset(k,     0,             state_dim * sizeof(double));

  gamma[0] = R + g[0] * f[0];
  D_bar[0] = D[0] * R / gamma[0];
  k[0] = g[0];
  U_bar[0] = 1;
  for (u32 j=1; j<state_dim; j++) {
    gamma[j] = gamma[j-1] + g[j] * f[j];
    D_bar[j] = D[j] * gamma[j-1] / gamma[j];
    double f_over_gamma = f[j] / gamma[j-1];
    for (u32 i=0; i<=j; i++) {
      U_bar[i*state_dim + j] = U[i*state_dim + j] - f_over_gamma * k[i]; // U_bar[:,j] = U[:,j] - f[j]/gamma[j-1] * k
      k[i] += g[j] * U[i*state_dim + j]; // k = k + g[j] * U[:,j]
    }
  }
  for (u32 i=0; i<state_dim; i++) {
    k[i] /= alpha;
  }
  memcpy(U, U_bar, state_dim * state_dim * sizeof(double));
  memcpy(D, D_bar,             state_dim * sizeof(double));

}

static void decorrelate(kf_t *kf, double *measurements, double *decor_measurements)
{
  u8 n_diffs = kf->obs_dim/2;
  memcpy(decor_measurements, measurements, kf->obs_dim * sizeof(double));
  cblas_dtrmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasUnit,
              n_diffs, kf->decor_mtx,
              n_diffs, decor_measurements, 1); // replaces raw phase measurements by their decorrelated version
  cblas_dtrmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasUnit,
              n_diffs, kf->decor_mtx,
              n_diffs, &decor_measurements[n_diffs], 1); // replaces raw measurements by their decorrelated version
}

/** In place updating of the state mean and covariance. Modifies measurements.
 */
void kalman_filter_update(kf_t *kf, double *measurements)
{
  double measurements_copy[kf->obs_dim];
  memcpy(measurements_copy, measurements, kf->obs_dim * sizeof(double));
  // VEC_PRINTF(measurements, kf->obs_dim);
  // MAT_PRINTF(kf->decor_obs_mtx, kf->obs_dim, kf->state_dim);
  // MAT_PRINTF(kf->obs_cov_root_inv, kf->obs_dim, kf->obs_dim);
  u8 n_diffs = kf->obs_dim/2;

  cblas_dtrmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasUnit,
              n_diffs, kf->decor_mtx,
              n_diffs, measurements_copy, 1); // replaces raw phase measurements by their decorrelated version
  cblas_dtrmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasUnit,
              n_diffs, kf->decor_mtx,
              n_diffs, &measurements_copy[n_diffs], 1); // replaces raw measurements by their decorrelated version

  predict_forward(kf);
  incorporate_obs(kf, measurements_copy);
}

void assign_transition_mtx(u32 state_dim, double dt, double *transition_mtx)
{
  eye(state_dim, transition_mtx);
  transition_mtx[3] = dt;
  transition_mtx[state_dim + 4] = dt;
  transition_mtx[2 * state_dim + 5] = dt;
}

void assign_d_mtx(u8 num_sats, double *D)
{
  memset(D, 0, (num_sats - 1) * num_sats * sizeof(double));
  //shape = num_sats-1  x  num_sats
  for(u8 i=0; i<num_sats-1; i++) {
    D[i * num_sats] = -1;
    D[i * num_sats + i + 1] = 1;
  }
}

void assign_e_mtx(u8 num_sats, sdiff_t *sats_with_ref_first, double ref_ecef[3], double *E)
{
  memset(E, 0, num_sats * 3 * sizeof(double));
  for (u8 i=0; i<num_sats; i++) {
    double *e = sats_with_ref_first[i].sat_pos;

    double x = e[0] - ref_ecef[0];
    double y = e[1] - ref_ecef[1];
    double z = e[2] - ref_ecef[2];
    double norm = sqrt(x*x + y*y + z*z);
    E[3*i] = x / norm;
    E[3*i + 1] = y / norm;
    E[3*i + 2] = z / norm;
  }
}

void assign_e_mtx_from_alms(u8 num_sats, almanac_t *alms, gps_time_t timestamp, double ref_ecef[3], double *E)
{
  memset(E, 0, num_sats * 3 * sizeof(double));
  // VEC_PRINTF(ref_ecef, 3);
  for (u8 i=0; i<num_sats; i++) {
    double e[3];
    double v[3];
    // calc_sat_state_almanac(&alms[i], timestamp.tow, -1, e, v); //TODO change -1 back to timestamp.wn
    calc_sat_state_almanac(&alms[i], timestamp.tow, timestamp.wn, e, v);
    // printf("\nprn=%u\n", alms[i].prn);
    // VEC_PRINTF(e, 3);
    double x = e[0] - ref_ecef[0];
    double y = e[1] - ref_ecef[1];
    double z = e[2] - ref_ecef[2];
    double norm = sqrt(x*x + y*y + z*z);
    E[3*i] = x / norm;
    E[3*i + 1] = y / norm;
    E[3*i + 2] = z / norm;
  }
}

// presumes that the first alm entry is the reference sat
void assign_de_mtx_from_alms(u8 num_sats, almanac_t *alms, gps_time_t timestamp, double ref_ecef[3], double *DE)
{
  memset(DE, 0, (num_sats - 1) * 3 * sizeof(double));
  double e0[3];
  double v0[3];
  calc_sat_state_almanac(&alms[0], timestamp.tow, timestamp.wn, e0, v0);
  // printf("\nprn=%u\n", alms[i].prn);
  // VEC_PRINTF(e, 3);
  double x0 = e0[0] - ref_ecef[0];
  double y0 = e0[1] - ref_ecef[1];
  double z0 = e0[2] - ref_ecef[2];
  double norm0 = sqrt(x0*x0 + y0*y0 + z0*z0);
  e0[0] = x0 / norm0;
  e0[1] = y0 / norm0;
  e0[2] = z0 / norm0;
  // VEC_PRINTF(ref_ecef, 3);
  for (u8 i=1; i<num_sats; i++) {
    double e[3];
    double v[3];
    // calc_sat_state_almanac(&alms[i], timestamp.tow, -1, e, v); //TODO change -1 back to timestamp.wn
    calc_sat_state_almanac(&alms[i], timestamp.tow, timestamp.wn, e, v);
    // printf("\nprn=%u\n", alms[i].prn);
    // VEC_PRINTF(e, 3);
    double x = e[0] - ref_ecef[0];
    double y = e[1] - ref_ecef[1];
    double z = e[2] - ref_ecef[2];
    double norm = sqrt(x*x + y*y + z*z);
    DE[3*(i-1)] = x / norm - e0[0];
    DE[3*(i-1) + 1] = y / norm - e0[1];
    DE[3*(i-1) + 2] = z / norm - e0[2];
  }
}

void assign_obs_mtx(u8 num_sats, sdiff_t *sats_with_ref_first, double ref_ecef[3], double *obs_mtx)
{
  u32 obs_dim = 2 * (num_sats-1);
  u32 state_dim = 5 + num_sats;

  memset(obs_mtx, 0, obs_dim * state_dim * sizeof(double));

  double DE[(num_sats-1) * 3];
  assign_de_mtx(num_sats, sats_with_ref_first, ref_ecef, &DE[0]);

  for (u32 i=0; i+1<num_sats; i++) {
    obs_mtx[i*state_dim] = DE[i*3] / GPS_L1_LAMBDA_NO_VAC;
    obs_mtx[i*state_dim+1] = DE[i*3+1] / GPS_L1_LAMBDA_NO_VAC;
    obs_mtx[i*state_dim+2] = DE[i*3+2] / GPS_L1_LAMBDA_NO_VAC;
    // obs_mtx[i*state_dim] = DE[i*3] / 0.190293673;
    // obs_mtx[i*state_dim+1] = DE[i*3+1] / 0.190293673;
    // obs_mtx[i*state_dim+2] = DE[i*3+2] / 0.190293673;

    obs_mtx[i*state_dim+6+i] = 1;

    obs_mtx[(i+num_sats-1)*state_dim] = DE[i*3];
    obs_mtx[(i+num_sats-1)*state_dim+1] = DE[i*3+1];
    obs_mtx[(i+num_sats-1)*state_dim+2] = DE[i*3+2];
  }
}

void assign_obs_mtx_from_alms(u8 num_sats, almanac_t *alms, gps_time_t timestamp, double ref_ecef[3], double *obs_mtx)
{
  u32 obs_dim = 2 * (num_sats-1);
  u32 state_dim = 5 + num_sats;

  memset(obs_mtx, 0, obs_dim * state_dim * sizeof(double));

  double DE[(num_sats-1) * 3];
  assign_de_mtx_from_alms(num_sats, alms, timestamp, ref_ecef, &DE[0]);

  for (u32 i=0; i+1<num_sats; i++) {
    obs_mtx[i*state_dim] = DE[i*3] / GPS_L1_LAMBDA_NO_VAC;
    obs_mtx[i*state_dim+1] = DE[i*3+1] / GPS_L1_LAMBDA_NO_VAC;
    obs_mtx[i*state_dim+2] = DE[i*3+2] / GPS_L1_LAMBDA_NO_VAC;
    // obs_mtx[i*state_dim] = DE[i*3] / 0.190293673;
    // obs_mtx[i*state_dim+1] = DE[i*3+1] / 0.190293673;
    // obs_mtx[i*state_dim+2] = DE[i*3+2] / 0.190293673;

    obs_mtx[i*state_dim+6+i] = 1;

    obs_mtx[(i+num_sats-1)*state_dim] = DE[i*3];
    obs_mtx[(i+num_sats-1)*state_dim+1] = DE[i*3+1];
    obs_mtx[(i+num_sats-1)*state_dim+2] = DE[i*3+2];
  }
}

void assign_decor_obs_cov(u8 num_diffs, double phase_var, double code_var,
                          double *decor_mtx, double *decor_obs_cov)
{
  memset(decor_mtx, 0, num_diffs * num_diffs * sizeof(double)); //this is for a single set of observations (only code or only carrier)
  memset(decor_obs_cov, 0, 2 * num_diffs * sizeof(double)); // this is the whole shebang (vector b/c independence)

  for (u8 i=0; i<num_diffs; i++) {
    double i_plus_one_divisor = 1.0 / (i+1.0);

    //assign the decorrelated covariance
    decor_obs_cov[i] = phase_var + phase_var * i_plus_one_divisor;
    decor_obs_cov[i+num_diffs] = code_var + code_var * i_plus_one_divisor;

    //assign the decorrelation matrix
    decor_mtx[i*num_diffs + i] = 1;
    for (u8 j=0; j<i; j++) {
      decor_mtx[i*num_diffs + j] = - i_plus_one_divisor;
    }
  }
}

void assign_decor_obs_mtx(u8 num_sats, sdiff_t *sats_with_ref_first,
                          double ref_ecef[3], double *decor_mtx, double *obs_mtx)
{
  u32 num_diffs = num_sats-1;
  u32 state_dim = num_diffs + 6;
  u32 obs_dim = 2 * num_diffs;
  memset(obs_mtx, 0, state_dim * obs_dim * sizeof(double));

  double DE[num_diffs * 3];
  assign_de_mtx(num_sats, sats_with_ref_first, ref_ecef, &DE[0]);
  //assign L^-1 * DE into DE
  cblas_dtrmm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, // CBLAS_ORDER, CBLAS_SIDE, CBLAS_UPLO, CBLAS_TRANSPOSE, CBLAS_DIAG
              num_diffs, 3, //M, N
              1, &decor_mtx[0], num_diffs, //alpha, A, lda
              &DE[0], 3); //B, ldb

  for (u32 i=0; i<num_diffs; i++) {
    obs_mtx[i*state_dim] = DE[i*3] / GPS_L1_LAMBDA_NO_VAC;
    obs_mtx[i*state_dim + 1] = DE[i*3 + 1] / GPS_L1_LAMBDA_NO_VAC;
    obs_mtx[i*state_dim + 2] = DE[i*3 + 2] / GPS_L1_LAMBDA_NO_VAC;

    memcpy(&obs_mtx[(i+num_diffs)*state_dim], &DE[i*3], 3 * sizeof(double));
    memcpy(&obs_mtx[i*state_dim + 6], &decor_mtx[i*num_diffs], (i+1) * sizeof(double));
  }
}

void assign_decor_obs_mtx_from_alms(u8 num_sats, almanac_t *alms, gps_time_t timestamp,
                                    double ref_ecef[3], double *decor_mtx, double *obs_mtx)
{
  u32 num_diffs = num_sats-1;
  u32 state_dim = num_diffs + 6;
  u32 obs_dim = 2 * num_diffs;
  memset(obs_mtx, 0, state_dim * obs_dim * sizeof(double));

  double DE[num_diffs * 3];
  assign_de_mtx_from_alms(num_sats, alms, timestamp, ref_ecef, &DE[0]);
  cblas_dtrmm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, // CBLAS_ORDER, CBLAS_SIDE, CBLAS_UPLO, CBLAS_TRANSPOSE, CBLAS_DIAG
              num_diffs, 3, //M, N
              1, &decor_mtx[0], num_diffs, //alpha, A, lda
              &DE[0], 3); //B, ldb

  for (u32 i=0; i<num_diffs; i++) {
    obs_mtx[i*state_dim] = DE[i*3] / GPS_L1_LAMBDA_NO_VAC;
    obs_mtx[i*state_dim + 1] = DE[i*3 + 1] / GPS_L1_LAMBDA_NO_VAC;
    obs_mtx[i*state_dim + 2] = DE[i*3 + 2] / GPS_L1_LAMBDA_NO_VAC;

    memcpy(&obs_mtx[(i+num_diffs)*state_dim], &DE[i*3], 3 * sizeof(double));
    memcpy(&obs_mtx[i*state_dim + 6], &decor_mtx[i*num_diffs], (i+1) * sizeof(double));
  }
}

static void least_squares_solve(kf_t *kf, double *measurements, double *lsq_state)
{
  double decor_measurements[kf->obs_dim];
  /*VEC_PRINTF(measurements, kf->obs_dim);*/
  decorrelate(kf, measurements, decor_measurements);
  s32 obs_dim = kf->obs_dim;
  s32 state_dim = kf->state_dim;
  s32 nrhs = 1;
  double decor_obs_mtx_transpose[kf->obs_dim * kf->state_dim];
  for (u32 i=0; i<kf->obs_dim; i++) {
    for (u32 j=0; j<kf->state_dim; j++) {
      decor_obs_mtx_transpose[i + j*kf->obs_dim] = kf->decor_obs_mtx[i*kf->state_dim + j];
    }
  }
  memcpy(lsq_state, decor_measurements, kf->obs_dim * sizeof(double));
  integer ldb = (s32) MAX(kf->state_dim, kf->obs_dim);
  double s[MIN(kf->state_dim, kf->obs_dim)];
  double rcond = 1e-12;
  integer rank;
  double w[1]; //try 25 + 10*num_sats
  integer lwork = -1;
  integer info;
  dgelss_((integer *) &obs_dim, (integer *) &state_dim, (integer *) &nrhs, //M, N, NRHS
                 &decor_obs_mtx_transpose[0], (integer *) &obs_dim, //A, LDA
                 &lsq_state[0], &ldb, //B, LDB
                 &s[0], &rcond, // S, RCOND
                 &rank, //RANK
                 &w[0], &lwork, // WORK, LWORK
                 &info); //INFO
  lwork = round(w[0]);

  double work[lwork];
  dgelss_((integer *) &obs_dim, (integer *) &state_dim, (integer *) &nrhs, //M, N, NRHS
                 &decor_obs_mtx_transpose[0], (integer *) &obs_dim, //A, LDA
                 &lsq_state[0], &ldb, //B, LDB
                 &s[0], &rcond, // S, RCOND
                 &rank, //RANK
                 &work[0], &lwork, // WORK, LWORK
                 &info); //INFO

  memset(&lsq_state[3],0,3 * sizeof(double)); //should already be nearly zero, because this bit of state is independent of the obs
}

void assign_transition_cov(u32 state_dim, double pos_var, double vel_var, double int_var, double *transition_cov)
{
  memset(transition_cov, 0, state_dim * state_dim * sizeof(double));
  transition_cov[0] = pos_var;
  transition_cov[state_dim + 1] = pos_var;
  transition_cov[2 * state_dim + 2] = pos_var;
  transition_cov[3 * state_dim + 3] = vel_var;
  transition_cov[4 * state_dim + 4] = vel_var;
  transition_cov[5 * state_dim + 5] = vel_var;
  for (u32 i=6; i<state_dim; i++) {
    transition_cov[i * state_dim + i] = int_var;
  }
}

void initialize_state(kf_t *kf, double *dd_measurements,
                      double pos_init_var, double vel_init_var, double int_init_var)
{
  double lsq_solution[MAX(kf->obs_dim,kf->state_dim)];
  least_squares_solve(kf, dd_measurements, lsq_solution);
  memcpy(kf->state_mean, lsq_solution, kf->state_dim * sizeof(double));
  eye(kf->state_dim, kf->state_cov_U);
  kf->state_cov_D[0] = pos_init_var;
  kf->state_cov_D[1] = pos_init_var;
  kf->state_cov_D[2] = pos_init_var;
  kf->state_cov_D[3] = vel_init_var;
  kf->state_cov_D[4] = vel_init_var;
  kf->state_cov_D[5] = vel_init_var;
  for (u32 i=6; i<kf->state_dim; i++) {
    kf->state_cov_D[i] = int_init_var;
  }
}

void get_kf(kf_t *kf, double phase_var, double code_var, double pos_var, double vel_var, double int_var,
            double pos_init_var, double vel_init_var, double int_init_var,
            u8 num_sdiffs, sdiff_t *sats_with_ref_first, double *dd_measurements, double ref_ecef[3], double dt)
{
  u32 state_dim = num_sdiffs + 5;
  u32 num_diffs = num_sdiffs-1;
  kf->state_dim = state_dim;
  kf->obs_dim = 2*num_diffs;
  assign_transition_mtx(state_dim, dt, &kf->transition_mtx[0]);
  assign_transition_cov(state_dim, pos_var, vel_var, int_var, &kf->transition_cov[0]);
  assign_decor_obs_cov(num_diffs, phase_var, code_var, &kf->decor_mtx[0], &kf->decor_obs_cov[0]);
  assign_decor_obs_mtx(num_sdiffs, sats_with_ref_first, &ref_ecef[0], &kf->decor_mtx[0], &kf->decor_obs_mtx[0]);
  initialize_state(kf, dd_measurements,
                   pos_init_var, vel_init_var, int_init_var);
}

void reset_kf_except_state(kf_t *kf,
                           double phase_var, double code_var,
                           double pos_var, double vel_var, double int_var,
                           u8 num_sdiffs, sdiff_t *sats_with_ref_first, double ref_ecef[3], double dt)
{
  u32 state_dim = num_sdiffs + 5;
  u32 num_diffs = num_sdiffs-1;
  kf->state_dim = state_dim;
  kf->obs_dim = 2*num_diffs;
  assign_transition_mtx(state_dim, dt, &kf->transition_mtx[0]);
  assign_transition_cov(state_dim, pos_var, vel_var, int_var, &kf->transition_cov[0]);
  assign_decor_obs_cov(num_diffs, phase_var, code_var, &kf->decor_mtx[0], &kf->decor_obs_cov[0]);
  assign_decor_obs_mtx(num_sdiffs, sats_with_ref_first, &ref_ecef[0], &kf->decor_mtx[0], &kf->decor_obs_mtx[0]);
}

kf_t get_kf_from_alms(double phase_var, double code_var, double pos_var, double vel_var, double int_var,
                      u8 num_sdiffs, almanac_t *alms, gps_time_t timestamp, double ref_ecef[3], double dt)
{
  u32 state_dim = num_sdiffs + 5;
  u32 num_diffs = num_sdiffs-1;
  kf_t kf;
  kf.state_dim = state_dim;
  kf.obs_dim = 2*num_diffs;
  assign_transition_mtx(state_dim, dt, &kf.transition_mtx[0]);
  assign_transition_cov(state_dim, pos_var, vel_var, int_var, &kf.transition_cov[0]);
  assign_decor_obs_cov(num_diffs, phase_var, code_var, &kf.decor_mtx[0], &kf.decor_obs_cov[0]);
  assign_decor_obs_mtx_from_alms(num_sdiffs, alms, timestamp, &ref_ecef[0], &kf.decor_mtx[0], &kf.decor_obs_mtx[0]);
  return kf;
}

void assign_state_rebase_mtx(u8 num_sats, u8 *old_prns, u8 *new_prns, double *rebase_mtx)
{
  u32 state_dim = num_sats + 5;
  memset(rebase_mtx, 0, state_dim * state_dim * sizeof(double));
  u8 old_ref = old_prns[0];
  u8 new_ref = new_prns[0];

  for (u8 i=0; i<6; i++) {
    rebase_mtx[i*state_dim + i] = 1;
  }
  s32 index_of_new_ref_in_old = find_index_of_element_in_u8s(num_sats-1, new_ref, &old_prns[1]);
  s32 index_of_old_ref_in_new = find_index_of_element_in_u8s(num_sats-1, old_ref, &new_prns[1]);
  for (u32 i=0; i<(u32)(num_sats-1); i++) {
    rebase_mtx[(i+6)*state_dim + index_of_new_ref_in_old+6] = -1;
    if (i != (u32) index_of_old_ref_in_new) {
      s32 index_of_this_sat_in_old_basis = find_index_of_element_in_u8s(num_sats-1, new_prns[i+1], &old_prns[1]);
      rebase_mtx[(i+6)*state_dim + index_of_this_sat_in_old_basis+6] = 1;
    }
  }
  // MAT_PRINTF(rebase_mtx, state_dim, state_dim);
}

static void old_rebase_covariance_sigma(double *state_cov, u8 num_sats, u8 *old_prns, u8 *new_prns)
{
  u32 state_dim = num_sats + 5;
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

static void old_rebase_covariance_udu(double *state_cov_U, double *state_cov_D, u8 num_sats, u8 *old_prns, u8 *new_prns)
{
  u32 state_dim = num_sats + 5;
  double state_cov[state_dim * state_dim];
  reconstruct_udu(state_dim, state_cov_U, state_cov_D, state_cov);
  old_rebase_covariance_sigma(state_cov, num_sats, old_prns, new_prns);
  udu(state_dim, state_cov, state_cov_U, state_cov_D);//
}

void rebase_kf(kf_t *kf, u8 num_sats, u8 *old_prns, u8 *new_prns)
{
  rebase_mean_N(&(kf->state_mean[6]), num_sats, old_prns, new_prns);
  old_rebase_covariance_udu(kf->state_cov_U, kf->state_cov_D, num_sats, old_prns, new_prns);
}

void kalman_filter_state_projection(kf_t *kf,
                                    u8 num_old_non_ref_sats,
                                    u8 num_new_non_ref_sats,
                                    u8 *ndx_of_new_sat_in_old)
{
  u32 old_state_dim = num_old_non_ref_sats + 6;
  double old_cov[old_state_dim * old_state_dim];
  reconstruct_udu(old_state_dim, kf->state_cov_U, kf->state_cov_D, old_cov);
  /*MAT_PRINTF(old_cov, old_state_dim, old_state_dim);*/

  u32 new_state_dim = num_new_non_ref_sats + 6;
  double new_cov[new_state_dim * new_state_dim];
  double new_mean[new_state_dim];

  //copy the bits that only relate to the baseline
  memcpy(new_mean, kf->state_mean, 6 * sizeof(double));
  for (u32 i=0; i<6; i++) {
    memcpy(&new_cov[i*new_state_dim], &old_cov[i*old_state_dim],  6 * sizeof(double));
  }

  //copy the bits that relate to the dd ambiguities
  for (u32 i=0; i<num_new_non_ref_sats; i++) {
    u32 ndxi = ndx_of_new_sat_in_old[i];
    new_mean[6+i] = kf->state_mean[6+ndxi];

    //copy the bits that are baseline-ambiguity covariances
    for (u32 j=0; j<6; j++) {
      new_cov[(6+i)*new_state_dim + j] = old_cov[(6+ndxi)*old_state_dim + j];
      new_cov[j*new_state_dim + 6+i] = old_cov[j*old_state_dim + 6+ndxi];
    }
    //copy the bits that are the pure ambiguity covariances
    for (u32 j=0; j<num_new_non_ref_sats; j++) {
      u32 ndxj = ndx_of_new_sat_in_old[j];
      new_cov[(i+6)*new_state_dim + j+6] = old_cov[(6+ndxi)*old_state_dim + 6+ndxj];
    }
  }
  /*MAT_PRINTF(new_cov, new_state_dim, new_state_dim);*/

  //put it all back into the kf
  memcpy(kf->state_mean, new_mean, new_state_dim * sizeof(double));
  udu(new_state_dim, new_cov, kf->state_cov_U, kf->state_cov_D);
  //NOTE: IT DOESN'T UPDATE THE OBSERVATION OR TRANSITION MATRICES, JUST THE STATE
}

void kalman_filter_state_inclusion(kf_t *kf,
                                   u8 num_old_non_ref_sats,
                                   u8 num_new_non_ref_sats,
                                   u8 *ndx_of_old_sat_in_new,
                                   double int_init_var)
{
  u32 old_state_dim = num_old_non_ref_sats + 6;
  double old_cov[old_state_dim * old_state_dim];
  reconstruct_udu(old_state_dim, kf->state_cov_U, kf->state_cov_D, old_cov);

  u32 new_state_dim = num_new_non_ref_sats + 6;
  double new_cov[new_state_dim * new_state_dim];
  memset(new_cov, 0, new_state_dim * new_state_dim * sizeof(double));
  for (u32 i=0; i<num_new_non_ref_sats; i++) {
    new_cov[(i+6)*new_state_dim + i+6] = int_init_var;
  }
  double new_mean[new_state_dim];
  // least_squares_solve(kf, dd_measurements, new_mean); //TODO do we really want to initialize new states this way?
  memset(new_mean, 0, new_state_dim * sizeof(double));

  //copy the bits that only relate to the baseline
  memcpy(new_mean, kf->state_mean, 6 * sizeof(double));
  for (u32 i=0; i<6; i++) {
    memcpy(&new_cov[i*new_state_dim], &old_cov[i*old_state_dim],  6 * sizeof(double));
  }
  for (u32 i=0; i<num_old_non_ref_sats; i++) {
    u8 ndxi = ndx_of_old_sat_in_new[i];
    new_mean[ndxi + 6] = kf->state_mean[i+6];
    for (u32 j=0; j<6; j++) {
      new_cov[(ndxi+6)*new_state_dim + j] = old_cov[(i+6)*old_state_dim + j];
      new_cov[j*new_state_dim + ndxi + 6] = old_cov[j*old_state_dim + i + 6];
    }
    for (u32 j=0;j<num_old_non_ref_sats; j++) {
      u8 ndxj = ndx_of_old_sat_in_new[j];
      new_cov[(ndxi+6)*new_state_dim + ndxj+6] = old_cov[(i+6)*old_state_dim + j+6];
    }
  }
  udu(new_state_dim, new_cov, kf->state_cov_U, kf->state_cov_D);
  memcpy(kf->state_mean, new_mean, new_state_dim * sizeof(double));
}


