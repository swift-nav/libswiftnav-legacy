# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

# cython: embedsignature=True

"""Functional Memory Pool

Simple fixed size memory pool collection supporting functional
operations.

The functional memory pool container is both a memory pool handling
allocation of fixed size 'elements' and a container type that arranges
allocated elements in a linked list, exposing functional style
primitives such as map and fold that operate over the
container. Elements can be removed from the container and released
back to the pool with a filter operation.

Allocation and deallocation from the pool are guaranteed constant time
and map and fold are O(N).

"""

from libc.string cimport memcpy, memset

cdef class MemoryPool:

  def __cinit__(self, n_elements, element):
    element_size = sizeof(element._thisptr)
    cdef memory_pool_t *pool = memory_pool_new(n_elements, element_size)
    self._thisptr = pool

  def empty(self):
    return memory_pool_empty(self._thisptr)

  def __dealloc__(self):
    memory_pool_destroy(self._thisptr)

  def n_free(self):
    return memory_pool_n_free(self._thisptr)

  def n_allocated(self):
    return memory_pool_n_allocated(self._thisptr)

  def n_elements(self):
    return memory_pool_n_elements(self._thisptr)

  def add(self):
    memory_pool_add(self._thisptr)

  def to_array(self):
    raise NotImplementedError

  def map(self, arg, f):
    raise NotImplementedError

  def filter(self, arg, f):
    raise NotImplementedError

  def clear(self):
    raise NotImplementedError

  def fold(self, x0, f):
    raise NotImplementedError

  def dfold(self, x0, f):
    raise NotImplementedError

  def ffold(self, x0, f):
    raise NotImplementedError

  def ifold(self, x0, f):
    raise NotImplementedError

  def sort(self, arg, cmp):
    raise NotImplementedError

  def group_by(self, arg, cmp, x0, x_size, agg):
    raise NotImplementedError

  def product(self, xs, max_xs, x_size, prod):
    raise NotImplementedError

  def product_generator(self, x0, n_xs, size_x_size, init, next, prod):
    raise NotImplementedError
