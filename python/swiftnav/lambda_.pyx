# Copyright (C) 2015 Swift Navigation Inc.
# Copyright (C) 2007-2008 by T.TAKASU, All rights reserved.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

# cython: embedsignature=True

import numpy as np
cimport numpy as np

def lambda_reduction_(Q):
  """Lambda reduction transformation: integer least-square
  estimation. reduction is performed by lambda and search by
  mlambda. See RTKLIB for more references.

  Parameters
  ----------
  Q : Numpy covariance matrix of float parameters (n x n)

  Returns
  ----------
  Z
  """
  assert len(Q.shape) == 2 and Q.shape[0] == Q.shape[1], "Q matrix must have shape (n, n)"
  n = Q.shape[0]
  cdef np.ndarray[np.double_t, ndim=2, mode="fortran"] Q_ = np.array(Q, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=2, mode="fortran"] Z = np.empty((n, n), dtype=np.double)
  assert lambda_reduction(n, &Q_[0,0], &Z[0,0]) == 0, "lambda error!"
  return Z

def lambda_solution_(x, sigma, m):
  """Performs integer least squares for mean x and covariance sigma,
  returning the m vectors with smallest residuals.

  Parameters
  ----------
  x : float parameters (n x 1)
  sigma : sigma covariance matrix of float parameters (n x n)
  m :  number of fixed solutions

  Returns
  ----------
  (fixed solutions (n x m), sum of squared residulas of fixed solutions (1 x m))

  """
  assert len(sigma.shape) == 2 and sigma.shape[0] == sigma.shape[1], "Q matrix must have shape (n, n)"
  assert len(x.shape) == 1 and x.shape[0] == sigma.shape[1], "x vector must have length to match sigma matrix"
  num_dds = x.shape[0]
  cdef np.ndarray[np.double_t, ndim=1, mode="fortran"] x_ = np.array(x, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=2, mode="fortran"] sigma_ = np.array(sigma, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=2, mode="fortran"] F_ = np.empty((num_dds,m), dtype=np.double, order='F')
  cdef np.ndarray[np.double_t, ndim=1, mode="fortran"] s_ = np.empty((m), dtype=np.double, order='F')
  assert lambda_solution(num_dds, m, &x_[0], &sigma_[0,0], &F_[0,0], &s_[0]) == 0,"lambda error!"
  return (F_.T, s_.T)
