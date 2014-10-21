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
#include <common.h>
#include <iterator.h>

typedef u8 prn;

typedef struct {
  u8 len;
  size_t size;
  const void *ref;
  const void *set; // TODO rename this to array?
} set_t;

typedef struct {
  const void *current;
  const set_t *set;
} set_state_t;

bool is_set(const set_t *set, key (*key)(const void *));
key prn_key(const void *x);
void print_prn_tuple(const void *arg, const void *elem);
void mk_prn_set(set_t *set, int len, prn *arr);
void mk_set_itr(iterator_t *it, set_state_t *s, set_t *set);

void print_prn(const void *arg, const void *elem);

bool eq_ref(const void *arg, const void *elem);
bool not_eq_ref(const void *arg, const void *elem);

s8 set_ref_prn(set_t *set, prn prn);
void mk_without_ref_itr(iterator_t *it, filter_state_t *s, iterator_t *base, set_t *set);
#endif /* LIBSWIFTNAV_SET_H */
