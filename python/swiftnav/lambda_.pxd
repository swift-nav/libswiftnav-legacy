# Copyright (C) 2014 Swift Navigation Inc.
# Copyright (C) 2007-2008 by T.TAKASU, All rights reserved.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

cdef extern from "libswiftnav/lambda.h":
  int lambda_reduction(int n, double *Q, double *Z)
  int lambda_solution(int n, int m, double *a, double *Q, double *F, double *s)
