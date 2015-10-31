# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *

cdef extern from "libswiftnav/memory_pool.h":

  ctypedef u8 element_t

  ctypedef struct memory_pool_node_hdr_t:
    node_t *next

  ctypedef struct node_t:
    memory_pool_node_hdr_t hdr
    element_t elem[0]

  ctypedef struct memory_pool_t:
    node_t *allocated_nodes_head

