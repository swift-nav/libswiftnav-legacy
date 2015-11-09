# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from libcpp cimport bool

cdef extern from "libswiftnav/set.h":

  ctypedef int (*cmp_fn)

  int cmp_u8_u8(const void * a, const void * b)
  int cmp_s32_s32(const void * a, const void * b)

  bool is_set(u8 n, size_t sz, const void *set, cmp_fn cmp)
  bool is_prn_set(u8 len, const u8 *prns)

  s32 intersection_map(u32 na, size_t sa, const void *as,
                       u32 nb, size_t sb, const void *bs,
                       cmp_fn cmp, void *context,
                       void (*f)(void *context, u32 n, const void *a, const void *b))

  s32 intersection(u32 na, size_t sa, const void *as, void *a_out,
                   u32 nb, size_t sb, const void *bs, void *b_out,
                   cmp_fn cmp)

  u32 insertion_index(u32 na, size_t sa, const void *as, void *b, cmp_fn cmp)
  u32 remove_element(u32 na, size_t sa, const void *as, void *a_out, void *b, cmp_fn cmp)
  u32 insert_element(u32 na, size_t sa, const void *as, void *a_out, void *b, cmp_fn cmp)
