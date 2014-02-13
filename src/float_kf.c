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

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

s8 udu(u32 n, double M[n][n], double U[n][n], double D[n]) 
{
  double alpha, beta;
  triu(n, M);
  eye(n, U);
  memset(D, 0, n * sizeof(double));

  for (u32 j=n; j>=2; j--) {
    D[j - 1] = M[j-1][j-1];
    if (D[j-1] > 0) {
      alpha = 1.0 / D[j-1];
    } else {
      alpha = 0.0;
    }
    for (u32 k=1; k<j; k++) {
      beta = M[k-1][j-1];
      U[k-1][j-1] = alpha * beta;
      for (u32 kk = 0; kk < k; kk++) {
        M[kk][k-1] = M[kk][k-1] - beta * U[kk][j-1];
      }
    }
    
  }
  D[0] = M[0][0];
  return 0;
}

void eye(u32 n, double M[n][n])
{
  memset(M, 0, n * n * sizeof(double));
  for (u32 i=0; i<n; i++) {
    M[i][i] = 1;
  }
}

void triu(u32 n, double M[n][n])
{
  for (u32 i=1; i<n; i++) {
    for (u32 j=0; j<i; j++) {
      M[i][j] = 0;
    }
  }
}

/** Reconsructs a UDU' decomposed matrix
 */
void reconstruct_udu(u32 n, double U[n][n], double D[n], double M[n][n]) {
  memset(M, 0, n * n * sizeof(double));

  for (u32 i=0; i<n; i++) {
    for (u32 k=0; k<n; k++) {
      for (u32 j=MAX(i,k); j<n; j++) {
        //U[i][j] is upper triangular = 0 if j < i
        //U[k][j] is upper triangular = 0 if j < k
        //U[i][j] * U[k][j] = 0 if j < k or j < i
        M[i][k] += U[i][j] * D[j] * U[k][j];
      }
    }
  }
}

// void filter_update(kf_t kf,
//                    double *state_mean, double *state_cov_U, double *state_cov_D, 
//                    double *raw_measurements) {
//   cblas_dtrmv(CblasRowMajor, L, N, N,
//               kf.obs_dim, kf.obs_cov_root_inv,
//               kf.obs_dim, raw_measurements)
// };







