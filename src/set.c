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

#include <string.h>
#include <assert.h>

#include "set.h"

/** \defgroup set Utility functions for dealing with sets
 * \{ */

/** Call a function for each element in the intersection of two sets.
 *
 * \param na Number of elements in set A
 * \param sa Size of each element of set A
 * \param as Array of elements in set A
 * \param nb Number of elements in set B
 * \param sb Size of each element of set B
 * \param as Array of elements in set B
 * \param cmp Pointer to a comparison function
 * \param context Pointer to an context passed directly through to `f`
 * \param f Pointer to function to map across intersection
 *
 * \return Number of elements in intersection
 */
u32 intersection_map(u32 na, size_t sa, const void *as,
                     u32 nb, size_t sb, const void *bs,
                     cmp_fn cmp, void *context,
                     void (*f)(void *context, u32 n,
                               const void *a, const void *b))
{
  u32 ia, ib, n = 0;

  for (ia=0, ib=0; ia<na && ib<nb; ia++, ib++) {
    const void *a = as + sa*ia;
    const void *b = bs + sb*ib;
    if (cmp(a, b) < 0) {
      ib--;
    } else if (cmp(a, b) > 0) {
      ia--;
    } else {
      f(context, n, a, b);
      n++;
    }
  }
  return n;
}

/** \} */

