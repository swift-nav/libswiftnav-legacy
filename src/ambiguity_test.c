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
 #include <stdio.h>
 #include <string.h>
 #include "ambiguity_test.h"
 #include "common.h"
 #include "linear_algebra.h"

// void init_residual_matrices(residual_mtxs_t *res_mtxs, u8 num_dds, double *DE_mtx, double *obs_covar)
// {

// }


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

void assign_phase_obs_null_basis(u8 num_sats, double *DE_mtx, double *q)
{
  // use either GEQRF or GEQP3. GEQP3 is the version with pivoting
  // int dgeqrf_(__CLPK_integer *m, __CLPK_integer *n, __CLPK_doublereal *a, __CLPK_integer *
  //       lda, __CLPK_doublereal *tau, __CLPK_doublereal *work, __CLPK_integer *lwork, __CLPK_integer *info)
  // int dgeqp3_(__CLPK_integer *m, __CLPK_integer *n, __CLPK_doublereal *a, __CLPK_integer *
  //       lda, __CLPK_integer *jpvt, __CLPK_doublereal *tau, __CLPK_doublereal *work, __CLPK_integer *lwork,
  //        __CLPK_integer *info)

  //DE is num_sats-1 by 3, need to transpose it to column major
  double A[(num_sats-1)*(num_sats-1)];
  for (u8 i=0; i+1<num_sats; i++) {
    for (u8 j=0; j<3; j++) {
      A[j*(num_sats-1) + i] = DE_mtx[i*3 + j]; //set A = Transpose(DE_mtx)
    }
  }
  integer m = num_sats - 1;
  integer n = 3;
  double tau[3];
  QR_part1(m, n, A, tau);
  QR_part2(m, n, A, tau);
  memcpy(q, &A[3*(num_sats-1)], (num_sats-4) * (num_sats-1) * sizeof(double));
}

// void assign_r_vec(residual_mtxs_t *res_mtxs, u8 num_dds, s8 *dd_measurements, double *r_vec)
// {

// }

// void assign_r_mean(residual_mtxs_t *res_mtxs, u8 num_dds, s8 *hypothesis, double *r_mean)
// {

// }

// double assign_quadratic_term(residual_mtxs_t *res_mtxs, u8 num_dds, s8 *hypothesis, double *r_vec)
// {

// }