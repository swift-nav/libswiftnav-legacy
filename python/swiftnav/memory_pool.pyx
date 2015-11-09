# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

cdef class MemoryPool:

  def __cinit__(self):
    raise NotImplementedError
    # memory_pool_t *memory_pool_new(u32 n_elements, size_t element_size)
    # s8 memory_pool_init(memory_pool_t *new_pool, u32 n_elements, size_t element_size, void *buff)

  def empty(self):
    #u8 memory_pool_empty(self):
    raise NotImplementedError

  def __dealloc__(self):
    raise NotImplementedError
   # void memory_pool_destroy(self)
   # s32 memory_pool_n_free(self)
   # s32 memory_pool_n_allocated(self)
   # u32 memory_pool_n_elements(self)


  def add(self):
    raise NotImplementedError

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
