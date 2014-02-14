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
void reconstruct_udu(u32 n, double *U, double *D, double *M) {
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

void predict_forward(kf_t *kf, double *state_mean, double *state_cov_U, double *state_cov_D) {
  double x[kf->state_dim];
  memcpy(x, state_mean, kf->state_dim * sizeof(double));

  cblas_dgemv(CblasRowMajor, CblasNoTrans, // CBLAS_ORDER, CBLAS_TRANSPOSE
              kf->state_dim, kf->state_dim, // int M, int N,
              1, (double *) kf->transition_mtx, kf->state_dim, // double 1, double *A, int lda
              x, 1, // double *X, int incX
              0, state_mean, 1); // double beta, double *Y, int incY

  double state_cov[kf->state_dim * kf->state_dim];
  reconstruct_udu(kf->state_dim, state_cov_U, state_cov_D, state_cov);

  // double *FC;
  // cblas_dsymm(CblasRowMajor, CblasRight, CblasUpper, //CBLAS_ORDER, CBLAS_SIDE, CBLAS_UPLO
  //             kf->state_dim, kf->state_dim, // int M, int N
  //             1, state_cov, kf->state_dim, // double alpha, double *A, int lda
  //             state_cov, kf->state_dim, // double *B, int ldb
  //             0, FC, kf->state_dim); // double beta, double *C, int ldc

  // double *FCF;
  // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, // CBLAS_ORDER, CBLAS_TRANSPOSE transA, cBLAS_TRANSPOSE transB
  //             kf->state_dim, kf->state_dim, kf->state_dim, // int M, int N, int K
  //             1, FC, kf->state_dim, // double alpha, double *A, int lda
  //             state_cov_U, kf->state_dim, //double *B, int ldb
  //             0, FCF, kf->state_dim); //beta, double *C, int ldc
}

void filter_update(kf_t *kf,
                   double *state_mean, double *state_cov_U, double *state_cov_D, 
                   double *raw_measurements) {
  cblas_dtrmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
              kf->obs_dim, kf->obs_cov_root_inv, 
              kf->obs_dim, raw_measurements, 1); // replaces raw_measurements by its decorrelated version





  memset(state_cov_U, 0, kf->state_dim * kf->state_dim * sizeof(double));
  memset(state_cov_D, 0,                kf->state_dim * sizeof(double));
  memset(state_mean,  0,                kf->state_dim * sizeof(double));
};







