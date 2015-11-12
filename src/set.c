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

#define CMP_FUNCTION(ta, tb)                        \
int cmp_##ta##_##tb(const void * a, const void * b) \
{                                                   \
  const ta *da = (const ta *) a;                    \
  const tb *db = (const tb *) b;                    \
  return (*da > *db) - (*da < *db);                 \
}

CMP_FUNCTION(s32, s32)

/* Tests if an array of PRNs form a sorted set with no duplicate elements.
 *
 * \param n Length of PRN array
 * \param prns Array of PRNs
 * \return `TRUE` if the PRNs form an ordered set, else `FALSE`
 */
bool is_sid_set(u8 n, const gnss_signal_t *sids)
{
  return is_set(n, sizeof(gnss_signal_t), sids, cmp_sid_sid);
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

/** Call a function for each element in the intersection of two sorted sets.
 *
 * Set A and B should be represented as sorted arrays with no repeated
 * elements.
 *
 * \param na Number of elements in set A
 * \param sa Size of each element of set A
 * \param as Array of elements in set A
 * \param nb Number of elements in set B
 * \param sb Size of each element of set B
 * \param bs Array of elements in set B
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

struct intersection_context {
  void *a_out;
  size_t sa;
  void *b_out;
  size_t sb;
};

static void intersection_function(void *context, u32 n,
                                  const void *a, const void *b)
{
  (void)n;
  struct intersection_context *ctxt = (struct intersection_context *)context;

  if (ctxt->a_out) {
    memcpy(ctxt->a_out, a, ctxt->sa);
    ctxt->a_out += ctxt->sa;
  }
  if (ctxt->b_out) {
    memcpy(ctxt->b_out, b, ctxt->sb);
    ctxt->b_out += ctxt->sb;
  }
}

/** Take the intersection of two sets.
 *
 * Takes two arrays each representing a set and outputs two arrays containing
 * the elemnts from the input arrays that are equal under a comparison function
 * `cmp`.
 *
 * Set A and B should be represented as sorted arrays with no repeated
 * elements.
 *
 * \param na    Number of elements in set A
 * \param sa    Size of each element of set A
 * \param as    Array of elements in set A
 * \param a_out Output array of matching elements from A
 * \param nb    Number of elements in set B
 * \param sb    Size of each element of set B
 * \param bs    Array of elements in set B
 * \param b_out Output array of matching elements from B
 * \param cmp    Pointer to a comparison function
 *
 * \return Number of elements in intersection on success,
 *         -1 if A is not a valid set,
 *         -2 if B is not a valid set
 */
s32 intersection(u32 na, size_t sa, const void *as, void *a_out,
                 u32 nb, size_t sb, const void *bs, void *b_out,
                 cmp_fn cmp)
{
  struct intersection_context ctxt = {
    .a_out = a_out,
    .sa = sa,
    .b_out = b_out,
    .sb = sb
  };

  return intersection_map(na, sa, as, nb, sb, bs,
                          cmp, &ctxt, intersection_function);
}

/** Given set and element, returns index where the element should be inserted.
 *
 * \param na    Number of elements in set A
 * \param sa    Size of each element of set A
 * \param as    Array of elements in set A
 * \param b     Element of interest
 * \param cmp   Pointer to a comparison function
 *
 * \return      Index where b belongs
 *              (points past end of as if b is largest.)
 */
u32 insertion_index(u32 na, size_t sa, const void *as, void *b, cmp_fn cmp)
{
  u32 index;
  for (index = 0; index < na && cmp(as + index * sa, b) < 0; index++);
  return index;
}

/** Removes first element that isn't less than b.
 *
 * \param na    Number of elements in set A
 * \param sa    Size of each element of set A
 * \param as    Array of elements in set A
 * \param a_out Output array
 * \param b     Element of interest
 * \param cmp   Pointer to a comparison function
 *
 * \return      Old index of removed element
 */
u32 remove_element(u32 na, size_t sa, const void *as, void *a_out,
                   void *b, cmp_fn cmp)
{
  /* Index of first element that isn't less than b.
   * This element will be removed. 
   * If b is larger than all of as, the last element of as will be removed. */
  u32 index = insertion_index(na, sa, as, b, cmp);

  memcpy(a_out, as, index * sa);
  memcpy(a_out + index * sa, as + (index + 1) * sa, (na - index - 1) * sa);

  return index;
}

/** Inserts element into set.
 *
 * \param na    Number of elements in set A
 * \param sa    Size of each element of set A
 * \param as    Array of elements in set A
 * \param a_out Output array
 * \param b     Element of interest
 * \param cmp   Pointer to a comparison function
 *
 * \return      Index of inserted element in as
 *              (may point past end of as)
 */
u32 insert_element(u32 na, size_t sa, const void *as, void *a_out,
                   void *b, cmp_fn cmp)
{
  u32 index = insertion_index(na, sa, as, b, cmp);

  memcpy(a_out, as, index * sa);
  memcpy(a_out + index * sa, b, sa);
  memcpy(a_out + (index + 1) * sa, as + index * sa, (na - index) * sa);

  return index;
}

/** \} */

