# Copyright (C) 2012 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from ephemeris_c cimport ephemeris_t
from libcpp cimport bool

cdef extern from "libswiftnav/nav_msg.h":
  ctypedef struct nav_msg_t:
    pass

  void nav_msg_init(nav_msg_t *n)
  s32 nav_msg_update(nav_msg_t *n, s32 corr_prompt_real)
  bool subframe_ready(nav_msg_t *n)
  void process_subframe(nav_msg_t *n, ephemeris_t *e)

