# Copyright (C) 2014 Swift Navigation Inc.
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
cimport lambda_c

def reduction(Q):
  if len(Q.shape) != 2 or Q.shape[0] != Q.shape[1]:
    raise ValueError("Q matrix must have shape (n, n)")
  n = Q.shape[0]

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] Q_ = \
    np.array(Q, dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] Z = \
    np.empty((n, n), dtype=np.double)

  lambda_c.lambda_reduction(n, &Q_[0,0], &Z[0,0])

  return Z

