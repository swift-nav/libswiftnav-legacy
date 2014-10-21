/*
 * Copyright (C) 2014 Swift Navigation Inc.
 * Contact: Scott Kovach <scott@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef LIBSWIFTNAV_ITERATOR_H
#define LIBSWIFTNAV_ITERATOR_H

#include <stdlib.h>
#include "common.h"

typedef s32 key;

typedef struct iterator_t {
  bool (*end)(const struct iterator_t*);
  void (*next)(struct iterator_t*);
  void (*reset)(struct iterator_t*);
  const void *(*current)(struct iterator_t*);
  void *state;
} iterator_t;

typedef struct {
  iterator_t *it;
  const void *arg;
  bool (*f)(const void *, const void *);
} filter_state_t;

typedef struct {
  const void *fst;
  const void *snd;
} tuple;
typedef struct {
  iterator_t *it1;
  iterator_t *it2;
  key (*key1)(const void *);
  key (*key2)(const void *);
  tuple current;
} intersection_state_t;

bool end(const iterator_t *it);
void next(iterator_t *it);
void reset(iterator_t *it);
const void *current(iterator_t *it);
bool more(const iterator_t *it);

/* Generic Iterators */
void mk_filter_itr(iterator_t *it, filter_state_t* s, bool (*f)(const void*, const void*),
                   const void *arg, iterator_t *base);
void mk_intersection_itr(iterator_t *it, intersection_state_t* s,
                         iterator_t *it1, iterator_t *it2,
                         key (*key1)(const void *),
                         key (*key2)(const void *));

/* Interator Utilities */
s8 freeze_itr(size_t elem_size, size_t max_len, void *buffer, iterator_t *it);
bool ptr_itr_equality(iterator_t *it1, iterator_t *it2);
void each(iterator_t *it, void (*f)(const void *, const void *), const void *arg);

#endif /* LIBSWIFTNAV_ITERATOR_H */
