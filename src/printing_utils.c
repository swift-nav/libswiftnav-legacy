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
 * \param m     The number of rows to be printed in the matrices.
 * \param n     The number of columns in the matrix.
 * \param mat1  The matrix to be printed.
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

/** Prints a s64 valued matrix.
 *
 * \param m     The number of rows to be printed in the matrices.
 * \param n     The number of columns in the matrix.
 * \param mat1  The matrix to be printed.
 */
static void print_s64_mtx(s64 *mat, u32 m, u32 n)
{
  for (u32 i=0; i < m; i++) {
    for (u32 j=0; j < n; j++) {
      printf("%"PRId64", ", mat[i*n + j]);
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


void print_hyp(void *arg, element_t *elem)
{
  u8 num_dds = *( (u8 *) arg );

  hypothesis_t *hyp = (hypothesis_t *)elem;
  printf("[");
  for (u8 i=0; i< num_dds; i++) {
    printf("%"PRId32", ", hyp->N[i]);
  }
  printf("]: %f\n", hyp->ll);
}

/* Utilities for debugging inclusion algorithm in ambiguity_test.c */
static void print_Z(s8 label, u8 full_dim, u8 new_dim, z_t * Z)
{
  printf("Z %i:\n", label);
  print_s64_mtx(Z, full_dim, new_dim);
}

void print_intersection_state(intersection_count_t *x)
{
  u8 full_dim = x->old_dim + x->new_dim;

  printf("itr lower bounds:\n");
  for (u8 i = 0; i < x->new_dim; i++) {
    printf("%"PRIi64"\n", x->itr_lower_bounds[i]);
  }
  printf("itr upper bounds:\n");
  for (u8 i = 0; i < x->new_dim; i++) {
    printf("%"PRIi64"\n", x->itr_upper_bounds[i]);
  }
  printf("box lower bounds:\n");
  for (u8 i = 0; i < full_dim; i++) {
    printf("%"PRIi64"\n", x->box_lower_bounds[i]);
  }
  printf("box upper bounds:\n");
  for (u8 i = 0; i < full_dim; i++) {
    printf("%"PRIi64"\n", x->box_upper_bounds[i]);
  }
  printf("transformation matrix:\n");
  print_Z(68, full_dim, x->new_dim, x->Z);
  printf("z1:\n");
  print_Z(0, full_dim, full_dim, x->Z1);
  printf("z2_inv:\n");
  print_Z(0, x->new_dim, x->new_dim, x->Z2_inv);
}

