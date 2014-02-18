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
#include <linear_algebra.h>
#include "gpstime.h"
#include "ephemeris.h"
#include "float_kf.h"

s8 udu(u32 n, double *M, double *U, double *D) 
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

/** Reconsructs a UDU' decomposed matrix
 */
void reconstruct_udu(u32 n, double *U, double *D, double *M) 
{
  memset(M, 0, n * n * sizeof(double));
  // TODO: will be symmetric, only need to bother populating part of it
  for (u32 i=0; i<n; i++) {
    for (u32 k=i; k<n; k++) {
      for (u32 j=k; j<n; j++) {
        //U[i][j] is upper triangular = 0 if j < i
        //U[k][j] is upper triangular = 0 if j < k
        //U[i][j] * U[k][j] = 0 if j < k or j < i
        M[i*n + k] += U[i*n +j] * D[j] * U[k*n + j];
      }
      M[k*n + i] = M[i*n + k];
    }
  }
}

/** In place prediction of the next state mean and covariances
 */
void predict_forward(kf_t *kf, double *state_mean, double *state_cov_U, double *state_cov_D) 
{
  //TODO take advantage of sparsity in this function

  double x[kf->state_dim];
  memcpy(x, state_mean, kf->state_dim * sizeof(double));

  cblas_dgemv(CblasRowMajor, CblasNoTrans, // CBLAS_ORDER, CBLAS_TRANSPOSE
              kf->state_dim, kf->state_dim, // int M, int N,
              1, (double *) kf->transition_mtx, kf->state_dim, // double 1, double *A, int lda
              x, 1, // double *X, int incX
              0, state_mean, 1); // double beta, double *Y, int incY
  // VEC_PRINTF((double *) state_mean, kf->state_dim);

  double state_cov[kf->state_dim * kf->state_dim];
  reconstruct_udu(kf->state_dim, state_cov_U, state_cov_D, state_cov);
  // MAT_PRINTF((double *) state_cov, kf->state_dim, kf->state_dim);

  double FC[kf->state_dim * kf->state_dim];
  cblas_dsymm(CblasRowMajor, CblasRight, CblasUpper, //CBLAS_ORDER, CBLAS_SIDE, CBLAS_UPLO
              kf->state_dim, kf->state_dim, // int M, int N
              1, state_cov, kf->state_dim, // double alpha, double *A, int lda
              kf->transition_mtx, kf->state_dim, // double *B, int ldb
              0, FC, kf->state_dim); // double beta, double *C, int ldc
  // MAT_PRINTF((double *) FC, kf->state_dim, kf->state_dim);

  double FCF[kf->state_dim * kf->state_dim];
  memcpy(FCF, kf->transition_cov, kf->state_dim * kf->state_dim * sizeof(double));
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, // CBLAS_ORDER, CBLAS_TRANSPOSE transA, cBLAS_TRANSPOSE transB
              kf->state_dim, kf->state_dim, kf->state_dim, // int M, int N, int K
              1, FC, kf->state_dim, // double alpha, double *A, int lda
              kf->transition_mtx, kf->state_dim, //double *B, int ldb
              1, FCF, kf->state_dim); //beta, double *C, int ldc
  // MAT_PRINTF((double *) FCF, kf->state_dim, kf->state_dim);

  udu(kf->state_dim, FCF, state_cov_U, state_cov_D);
  // MAT_PRINTF((double *) state_cov_U, kf->state_dim, kf->state_dim);
  // VEC_PRINTF((double *) state_cov_D, kf->state_dim);
}

/** In place updating of the state mean and covariances to use the (decorrelated) observations
 */
void update_for_obs(kf_t *kf,
                    double *intermediate_mean, double *intermediate_cov_U, double *intermediate_cov_D,
                    double *decor_obs)
{

  for (u32 i=0; i<kf->obs_dim; i++) {
    double *h = &kf->decor_obs_mtx[kf->state_dim * i]; //vector of length kf->state_dim
    double R = kf->decor_obs_cov[i]; //scalar
    double k[kf->state_dim]; // vector of length kf->state_dim
    // printf("i=%i\n", i);
    // VEC_PRINTF(h, kf->state_dim);

    update_scalar_measurement(kf->state_dim, h, R, intermediate_cov_U, intermediate_cov_D, &k[0]); //updates cov and sets k
    // VEC_PRINTF(k, kf->state_dim);

    double predicted_obs = 0;
    for (u32 j=0; j<kf->state_dim; j++) {//TODO take advantage of sparsity of h
      predicted_obs += h[j] * intermediate_mean[j];
    }
    double obs_minus_predicted_obs = decor_obs[i] - predicted_obs;
    // printf("decor_obs = %f\n", decor_obs[i]);
    // printf("predi_obs = %f\n", predicted_obs);
    // printf(" diff_obs = %f\n", obs_minus_predicted_obs);

    for (u32 j=0; j<kf->state_dim; j++) {
      intermediate_mean[j] += k[j] * obs_minus_predicted_obs; // uses k to update mean
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

/** In place updating of the state mean and covariance. Modifies measurements.
 */
void filter_update(kf_t *kf,
                   double *state_mean, double *state_cov_U, double *state_cov_D, 
                   double *measurements)
{
  // VEC_PRINTF(measurements, kf->obs_dim);
  // MAT_PRINTF(kf->decor_obs_mtx, kf->obs_dim, kf->state_dim);
  // MAT_PRINTF(kf->obs_cov_root_inv, kf->obs_dim, kf->obs_dim);
  cblas_dtrmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
              kf->obs_dim, kf->obs_cov_root_inv, 
              kf->obs_dim, measurements, 1); // replaces raw measurements by its decorrelated version

  predict_forward(kf, state_mean, state_cov_U, state_cov_D);
  update_for_obs(kf, state_mean, state_cov_U, state_cov_D, measurements);
}

void assign_transition_mtx(u32 state_dim, double dt, double *transition_mtx)
{
  eye(state_dim, transition_mtx);
  transition_mtx[3] = dt;
  transition_mtx[state_dim + 4] = dt;
  transition_mtx[2 * state_dim + 5] = dt;
}

kf_t get_kf(u8 num_sats, u8 *sats_with_ref_first, ephemeris_t *ephemerides, double *ref_ecef, gps_time_t timestamp, double dt)
{
  u32 state_dim = num_sats + 5;
  double transition_mtx[state_dim * state_dim];
  assign_transition_mtx(state_dim, 1, transition_mtx);
  memset(sats_with_ref_first, 0, 1);
  memset(ephemerides, 0, 1);
  memset(ref_ecef, 0, 1);
  printf("%f\n", timestamp.tow);

  dt += 0;
  kf_t kf;
  return kf;
}







