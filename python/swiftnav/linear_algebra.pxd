# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *

cdef extern from "libswiftnav/linear_algebra.h":

  void dmtx_printf(double *mtx, u32 m, u32 n);
  void dmtx_printi(s32 *mtx, u32 m, u32 n);
  void submatrix(u32 new_rows, u32 new_cols, u32 old_cols, const double *old,
                 const u32 *new_row_to_old, const u32 *new_col_to_old,
                 double *new);
  void submatrix_ul(u32 new_rows, u32 new_cols, u32 old_cols, const double *old,
                    double *new);
  s32 qrdecomp_square(const double *a, u32 rows, double *qt, double *r);
  s32 qrdecomp(const double *a, u32 rows, u32 cols, double *qt, double *r);
  void qtmult(const double *qt, u32 n, const double *b, double *x);
  void rsolve(const double *r, u32 rows, u32 cols, const double *b, double *x);
  s32 qrsolve(const double *a, u32 rows, u32 cols, const double *b, double *x);

  int matrix_inverse(u32 n, const double *const a, double *b);
  void matrix_multiply(u32 n, u32 m, u32 p, const double *a,
                       const double *b, double *c);
  void matrix_multiply_i(u32 n, u32 m, u32 p, const s32 *a,
                         const s32 *b, s32 *c);
  void matrix_multiply_s64(u32 n, u32 m, u32 p, const s64 *a,
                           const s64 *b, s64 *c);
  void matrix_triu(u32 n, double *M);
  void matrix_eye(u32 n, double *M);
  void matrix_udu(u32 n, double *M, double *U, double *D);
  void matrix_reconstruct_udu(const u32 n, const double *U, const double *D, double *M);
  void matrix_add_sc(u32 n, u32 m, const double *a,
                     const double *b, double gamma, double *c);
  void matrix_transpose(u32 n, u32 m, const double *a, double *b);
  void matrix_copy(u32 n, u32 m, const double *a, double *b);

  int matrix_pseudoinverse(u32 n, u32 m, const double *a, double *b);
  int matrix_atwaiat(u32 n, u32 m, const double *a, const double *w, double *b);
  int matrix_ataiat(u32 n, u32 m, const double *a, double *b);
  int matrix_atawati(u32 n, u32 m, const double *a, const double *w, double *b);
  int matrix_ataati(u32 n, u32 m, const double *a, double *b);

  double vector_dot(u32 n, const double *a, const double *b);
  double vector_norm(u32 n, const double *a);
  double vector_mean(u32 n, const double *a);
  void vector_normalize(u32 n, double *a);
  void vector_add_sc(u32 n, const double *a, const double *b,
                     double gamma, double *c);
  void vector_add(u32 n, const double *a, const double *b, double *c);
  void vector_subtract(u32 n, const double *a,
                       const double *b, double *c);
  void vector_cross(const double a[3], const double b[3], double c[3]);
  double vector_distance(u32 n, const double *a, const double *b);
