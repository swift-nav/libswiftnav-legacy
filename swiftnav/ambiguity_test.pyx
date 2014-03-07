# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

cimport ambiguity_test_c
import numpy as np
cimport numpy as np

def get_phase_obs_null_basis(DE):
  num_dds = DE.shape[0]
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] DE_ = \
    np.array(DE, dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] q_ = \
    np.empty((num_dds-3,num_dds), dtype=np.double)

  ambiguity_test_c.assign_phase_obs_null_basis(DE.shape[0]+1, &DE_[0,0], &q_[0,0])
  return q_