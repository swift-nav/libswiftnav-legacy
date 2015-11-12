# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

cimport numpy as np
cimport linear_algebra_c

class UDU_decomposition:
  def __init__(self, U, D):
    self.U = U
    self.D = D
  def reconstruct(self):
    return reconstruct_udu(self)

def udu(M):
  n = M.shape[0]
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] M_ = \
    np.array(M, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] U = \
    np.empty((n,n), dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] D = \
    np.empty(n, dtype=np.double)
  linear_algebra_c.matrix_udu(n, <double *> &M_[0,0], <double *> &U[0,0], <double *> &D[0])
  return UDU_decomposition(U, D)

def reconstruct_udu(ud):
  n = ud.D.shape[0]

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] U = \
    np.array(ud.U, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] D = \
    np.array(ud.D, dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] M = \
    np.empty((n,n), dtype=np.double)

  linear_algebra_c.matrix_reconstruct_udu(n, <double *> &U[0,0], <double *> &D[0], <double *> &M[0,0])
  return M