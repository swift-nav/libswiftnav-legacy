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

#include <stdio.h>
#include <string.h>
#include "constants.h"
#include "common.h"
#include "ambiguity_test.h"
#include "printing_utils.h"

/** Prints a matrix of doubles.
 *
 * Used primarily for printing matrices in a gdb session.
 * Assumes row-major (C-style) storage.
 *
 * \param m     The matrix to be printed.
 * \param _r    The number of rows to be printed.
 * \param _c    The number of columns in the matrix.
 */
void print_double_mtx(double *m, u32 _r, u32 _c)
{
    for (u32 _i = 0; _i < (_r); _i++) {
      printf(" [% 12lf", (m)[_i*(_c) + 0]);
      for (u32 _j = 1; _j < (_c); _j++)
        printf(" % 12lf", (m)[_i*(_c) + _j]);
      printf("]\n");
    }
}

/** Prints the matrix of Pearson correlation coefficients of a covariance matrix.
 *
 * Takes in a double valued covariance matrix and prints the pearson correlation
 * coefficients of each pair of variables:
 *
 * \f[
 *    \rho_{ij} = \frac{cov_{ij}}{\sigma_i \sigma_j}
 * \f]
 *
 * \param m     The covariance matrix to be transformed.
 * \param dim   The dimension of the covariance matrix.
 */
void print_pearson_mtx(double *m, u32 dim)
{
    for (u32 _i = 0; _i < dim; _i++) {
      printf(" [% 12lf", m[_i*dim + 0] / sqrt(m[_i*dim + _i]) / sqrt(m[0]));
      for (u32 _j = 1; _j < dim; _j++)
        printf(" % 12lf", m[_i*dim + _j] / sqrt(m[_i*dim + _i]) / sqrt(m[_j * dim + _j]));
      printf("]\n");
    }
}

/** Prints the difference between two s32 valued matrices.
 *
 * Given two s32 valued matrices of the same shape, prints their difference
 * (first matrix minus second).
 *
 * \param m     The number of rows in the matrices.
 * \param n     The number of columns in the matrices.
 * \param mat1  The matrix to be subtracted from.
 * \param mat2  The matrix to be subtracted.
 */
void print_s32_mtx_diff(u32 m, u32 n, s32 *mat1, s32 *mat2)
{
  for (u32 i=0; i < m; i++) {
    for (u32 j=0; j < n; j++) {
      printf("%"PRId32", ", mat1[i*n + j] - mat2[i*n + j]);
    }
    printf("\n");
  }
  printf("\n");
}

/** Prints a s32 valued matrix.
 *
 * \param mat  The matrix to be printed.
 * \param m     The number of rows to be printed in the matrices.
 * \param n     The number of columns in the matrix.
 */
void print_s32_mtx(s32 *mat, u32 m, u32 n)
{
  for (u32 i=0; i < m; i++) {
    for (u32 j=0; j < n; j++) {
      printf("%"PRId32", ", mat[i*n + j]);
    }
    printf("\n");
  }
  printf("\n");
}

/** Prints the result of a matrix-vector product of s32's.
 * Given a matrix M and vector v of s32's, prints the result of M*v.
 *
 * \param m   The number of rows in the matrix M.
 * \param n   The number of columns in the matrix M and elements of v.
 * \param M   The matrix M to be multiplied.
 * \param v   The vector v to be multiplied.
 */
void print_s32_gemv(u32 m, u32 n, s32 *M, s32 *v)
{
  s32 mv[m];
  memset(mv, 0, m * sizeof(s32));
  printf("[");
  for (u32 i=0; i < m; i++) {
    for (u32 j=0; j < n; j++) {
      mv[i] += M[i*n + j] * v[j];
    }
    if (i+1 == m) {
      printf("%"PRId32" ]\n", mv[i]);
    }
    else {
      printf("%"PRId32", ", mv[i]);
    }
  }
}


