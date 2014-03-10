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


