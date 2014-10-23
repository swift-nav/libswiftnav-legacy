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

// TODO move sdiff set definitions to single_diff?
#include "single_diff.h"

typedef struct {
  u8 len;
  size_t size;
  const void *set;
} set_t;

typedef struct {
  set_t *set;
  const void *ref;
} ptd_set_t;

typedef struct {
  const void *current;
  const set_t *set;
} set_state_t;

void mk_pointed(set_t *set, int len, size_t size, const void *arr);
void freeze_to_set(size_t max_len, set_t *set, iterator_t *it);

bool is_set(const set_t *set, key (*key)(const void *));
key prn_key(const void *x);
key sdiff_key(const void *x);
void print_prn_tuple(const void *arg, const void *elem);
void mk_sdiff_set(set_t *set, int len, sdiff_t *sats);
void mk_prn_set(set_t *set, int len, prn *arr);
void mk_set_itr(iterator_t *it, set_state_t *s, set_t *set);

void print_prn(const void *arg, const void *elem);

bool eq_ref(const void *arg, const void *elem);
bool not_eq_ref(const void *arg, const void *elem);

void set_ref_prn(ptd_set_t *ptd, prn prn, key (*f)(const void *));
void mk_without_ref_itr(iterator_t *it, filter_state_t *s, iterator_t *base, set_t *set);

/* MACROS */
#define MK_PRN_ITR(name, len) \
  set_t name##_set; \
  mk_prn_set(&name##_set, len, name); \
  DECL_ITR(set, name) \
  mk_set_itr(&name##_itr, &name##_state, &name##_set)

#define DECL_ITR(type, name) \
  iterator_t name##_itr; \
  type##_state_t name##_state; \

#endif /* LIBSWIFTNAV_SET_H */
