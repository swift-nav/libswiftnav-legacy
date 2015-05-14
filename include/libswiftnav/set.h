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

/** \addtogroup set
 * \{ */

/** Comparison function prototype.
 * Follows the standard C library comparison prototype e.g. used by qsort()
 * http://www.gnu.org/software/libc/manual/html_node/Comparison-Functions.html
 * */
typedef int (*cmp_fn) (const void*, const void*);

/** \} */

u32 intersection_map(u32 na, size_t sa, const void *as,
                     u32 nb, size_t sb, const void *bs,
                     cmp_fn cmp, void *context,
                     void (*f)(void *context, u32 n,
                               const void *a, const void *b));

#endif /* LIBSWIFTNAV_SET_H */

