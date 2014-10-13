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

#ifndef LIBSWIFTNAV_SET_H
#define LIBSWIFTNAV_SET_H

#include <stdlib.h>
#include "common.h"

typedef u8 prn;
typedef prn key;

typedef struct {
  u8 len;
  size_t size;
  void *ref;
  void *set;
  key (*key)(const void *);
} set_t;

typedef struct iterator_t {
  bool (*end)(const struct iterator_t*);
  void (*next)(struct iterator_t*);
  void *(*current)(struct iterator_t*);
  void *state;
} iterator_t;

bool is_set(const set_t *set);
void next(iterator_t *it);
bool end(const iterator_t *it);
bool more(const iterator_t *it);
void *current(iterator_t *it);

#endif /* LIBSWIFTNAV_SET_H */
