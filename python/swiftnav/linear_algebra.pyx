# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

# cython: embedsignature=True

cimport numpy as np

class UDU_decomposition:

  def __init__(self, U, D):
    self.U = U
    self.D = D

  def reconstruct(self):
    return matrix_reconstruct_udu_(self)

def dmtx_printf_(mtx, m, n):
  raise NotImplementedError

def dmtx_printi_(mtx, m, n):
  raise NotImplementedError

def submatrix_(new_rows, new_cols, old_cols, old, new_row_to_old, new_col_to_old, new_):
  raise NotImplementedError

def submatrix_ul_(new_rows, new_cols, old_cols, old, new_):
  raise NotImplementedError

def qrdecomp_square_(a, rows, qt, r):
  raise NotImplementedError

def qrdecomp_(a, rows, cols, qt, r):
  raise NotImplementedError

def qtmult_(qt, n, b, x):
  raise NotImplementedError

def rsolve_(r, rows, cols, b, x):
  raise NotImplementedError

def qrsolve_(a, rows, cols, b, x):
  raise NotImplementedError

def matrix_inverse_(n, a, b):
  raise NotImplementedError

def matrix_multiply_(n, m, p, a, b, c):
  raise NotImplementedError

def matrix_multiply_i_(n, m, p, a, b, c):
  raise NotImplementedError

def matrix_multiply_s64_(n, m, p, a, b, c):
  raise NotImplementedError

def matrix_triu_(M):
  raise NotImplementedError

def matrix_eye_(n):
  raise NotImplementedError

def matrix_udu_(M):
  n = M.shape[0]
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] M_ = np.array(M, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] U = np.empty((n,n), dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] D = np.empty(n, dtype=np.double)
  matrix_udu(n, <double *> &M_[0,0], <double *> &U[0,0], <double *> &D[0])
  return UDU_decomposition(U, D)

def matrix_reconstruct_udu_(ud):
  n = ud.D.shape[0]
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] U = np.array(ud.U, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] D = np.array(ud.D, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] M = np.empty((n,n), dtype=np.double)
  matrix_reconstruct_udu(n, <double *> &U[0,0], <double *> &D[0], <double *> &M[0,0])
  return M

def matrix_add_sc_(n, m, a, b, double gamma, c):
  raise NotImplementedError

def matrix_transpose_(n, m, a, b):
  raise NotImplementedError

def matrix_copy_(n, m, a, b):
  raise NotImplementedError

def matrix_pseudoinverse_(n, m, a, b):
  raise NotImplementedError

def matrix_atwaiat_(n, m, a, w, b):
  raise NotImplementedError

def matrix_ataiat_(n, m, a, b):
  raise NotImplementedError

def matrix_atawati_(n, m, a, w, b):
  raise NotImplementedError

def matrix_ataati_(n, m, a, b):
  raise NotImplementedError

def vector_dot_(n, a, b):
  raise NotImplementedError

def vector_norm_(n, a):
  raise NotImplementedError

def vector_mean_(n, a):
  raise NotImplementedError

def vector_normalize_(n, a):
  raise NotImplementedError

def vector_add_sc_(n, a, b, gamma, c):
  raise NotImplementedError

def vector_add_(n, a, b, c):
  raise NotImplementedError

def vector_subtract_(n, a, b, c):
  raise NotImplementedError

def vector_cross_(a, b, c):
  raise NotImplementedError

def vector_distance_(n, a, b):
  raise NotImplementedError
