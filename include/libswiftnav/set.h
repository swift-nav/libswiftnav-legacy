/*
 * Copyright (C) 2015 Swift Navigation Inc.
 * Contact: Fergus Noble <fergus@swift-nav.com>
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

#include "common.h"
#include "signal.h"

/** \addtogroup set
 * \{ */

/** Comparison function prototype.
 * Follows the standard C library comparison prototype e.g. used by qsort()
 * http://www.gnu.org/software/libc/manual/html_node/Comparison-Functions.html
 *
 * \param a Pointer to A
 * \param b Pointer to B
 * \return <0 if A < B,
 *          0 if A = B,
 *         >0 if A > B
 * */
typedef int (*cmp_fn) (const void* a, const void* b);

/** \} */

int cmp_s32_s32(const void * a, const void * b);

bool is_set(u8 n, size_t sz, const void *set, cmp_fn cmp);
bool is_prn_set(u8 len, const gnss_signal_t *prns);

s32 intersection_map(u32 na, size_t sa, const void *as,
                     u32 nb, size_t sb, const void *bs,
                     cmp_fn cmp, void *context,
                     void (*f)(void *context, u32 n,
                               const void *a, const void *b));

s32 intersection(u32 na, size_t sa, const void *as, void *a_out,
                 u32 nb, size_t sb, const void *bs, void *b_out,
                 cmp_fn cmp);

u32 insertion_index(u32 na, size_t sa, const void *as, void *b, cmp_fn cmp);
u32 remove_element(u32 na, size_t sa, const void *as, void *a_out,
                   void *b, cmp_fn cmp);
u32 insert_element(u32 na, size_t sa, const void *as, void *a_out,
                   void *b, cmp_fn cmp);
#endif /* LIBSWIFTNAV_SET_H */

