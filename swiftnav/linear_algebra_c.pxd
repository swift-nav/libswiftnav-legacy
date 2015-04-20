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
  void matrix_triu(u32 n, double *M)
  void matrix_eye(u32 n, double *M)
  void matrix_udu(u32 n, double *M, double *U, double *D)
  void matrix_reconstruct_udu(u32 n, double *U, double *D, double *M)

