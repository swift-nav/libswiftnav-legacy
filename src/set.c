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

/** Comparison function for s32s.
 * See `cmp_fn`. */
int cmp_s32(const void * a, const void * b)
{
   return *(s32*)a - *(s32*)b;
}

/** Comparison function for u8s.
 * See `cmp_fn`. */
int cmp_u8(const void * a, const void * b)
{
   return *(u8*)a - *(u8*)b;
}

/* Tests if an array of PRNs form a sorted set with no duplicate elements.
 *
 * \param n Length of PRN array
 * \param prns Array of PRNs
 * \return `TRUE` if the PRNs form an ordered set, else `FALSE`
 */
bool is_prn_set(u8 n, const u8 *prns)
{
  return is_set(n, sizeof(u8), prns, cmp_u8);
}

/* Tests if an array forms a sorted set with no duplicate elements.
 *
 * \param n Length of array
 * \param sz Size of array element in bytes
 * \param set Pointer to array
 * \param cmp Pointer to a comparison function
 * \return `TRUE` if the array forms an ordered set, else `FALSE`
 */
bool is_set(u8 n, size_t sz, const void *set, cmp_fn cmp)
{
  assert(sz != 0);
  assert(set != NULL);
  assert(cmp != NULL);

  if (n == 0) {
    return true;
  }

  const void *current = set;
  for (u8 i = 1; i < n; i++) {
    const void *next = set + sz*i;
    if (cmp(next, current) <= 0) {
      return false;
    }
    current = next;
  }
  return true;
}

/** Call a function for each element in the intersection of two sets.
 *
 * Set A and B should be represented as sorted arrays with no repeated
 * elements.
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
 * \return Number of elements in intersection on success,
 *         -1 if A is not a valid set,
 *         -2 if B is not a valid set
 */
s32 intersection_map(u32 na, size_t sa, const void *as,
                     u32 nb, size_t sb, const void *bs,
                     cmp_fn cmp, void *context,
                     void (*f)(void *context, u32 n,
                               const void *a, const void *b))
{
  assert(sa != 0);
  assert(sb != 0);
  assert(as != NULL);
  assert(bs != NULL);
  assert(cmp != NULL);
  assert(f != NULL);

  if (!is_set(na, sa, as, cmp)) {
    return -1;

  }
  if (!is_set(nb, sb, bs, cmp)) {
    return -2;
  }

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

