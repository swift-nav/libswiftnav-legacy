/*
 * Copyright (C) 2014 Swift Navigation Inc.
 * Contact: Fergus Noble <fergus@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#include <stdlib.h>
#include <string.h>

#include "memory_pool.h"

inline static size_t calc_node_size(size_t element_size)
{
  return element_size + sizeof(memory_pool_node_hdr_t);
}

inline static node_t *get_node_n(memory_pool_t *pool, node_t *head, u32 n)
{
  return (node_t *)((u8 *)head + calc_node_size(pool->element_size) * n);
}

/** \defgroup memory_pool Functional Memory Pool
 * Simple fixed size memory pool collection supporting functional operations.
 *
 * The functional memory pool container is both a memory pool handling
 * allocation of fixed size 'elements' and a container type that arranges
 * allocated elements in a linked list, exposing functional style primitives
 * such as map and fold that operate over the container. Elements can be
 * removed from the container and released back to the pool with a filter
 * operation.
 *
 * Allocation and deallocation from the pool are guaranteed constant time and
 * map and fold are O(N).
 *
 * \{ */

/** Create a new memory pool.
 * Creates a new memory pool containing a maximum of `n_elements` elements of
 * size `element_size`. This function calles malloc() to reserve space for a
 * ::memory_pool_t struct and for the memory pool itself. Each element stored
 * in the memory pool has an overhead of one pointer size so the total space
 * used for the pool will be:
 *
 * ~~~
 * n_elements * (element_size + sizeof(void *))
 * ~~~
 *
 * Remember to free the pool with memory_pool_destroy(), and for embedded
 * systems it is recommended that all pools be initialized once at start-up to
 * prevent all the usual caveats associated with dynamic memory allocation.
 *
 * \param n_elements Number of elements that the pool can hold
 * \param element_size Size in bytes of the user payload elements
 * \returns Pointer to a new ::memory_pool_t or NULL upon a malloc() failure
 */
memory_pool_t *memory_pool_new(u32 n_elements, size_t element_size)
{
  memory_pool_t *new_pool = malloc(sizeof(memory_pool_t));
  if (!new_pool) {
    return NULL;
  }

  /* Allocate memory pool */
  size_t node_size = calc_node_size(element_size);
  node_t *buff = (node_t *)malloc(node_size * n_elements);
  if (!buff) {
    free(new_pool);
    return NULL;
  }

  memory_pool_init(new_pool, n_elements, element_size, buff);

  return new_pool;
}

/** Initialise a new memory pool.
 * Initialises a new memory pool containing a maximum of `n_elements` elements
 * of size `element_size`. This function does not allocate memory and must be
 * passed a buffer of a suitable size to hold the memory pool elements. Each
 * element stored in the memory pool has an overhead of one pointer size so the
 * total space used for the pool will be:
 *
 * ~~~
 * n_elements * (element_size + sizeof(void *))
 * ~~~
 *
 * Remember to free the pool with memory_pool_destroy(), and for embedded
 * systems it is recommended that all pools be initialized once at start-up to
 * prevent all the usual caveats associated with dynamic memory allocation.
 *
 * \param pool Pointer to a memory pool to initialise
 * \param n_elements Number of elements that the pool can hold
 * \param element_size Size in bytes of the user payload elements
 * \param buff Pointer to a buffer to use as the memory pool working area
 * \returns `0` on success, `<0` on failure.
 */
s8 memory_pool_init(memory_pool_t *new_pool, u32 n_elements,
                    size_t element_size, void *buff)
{
  if (!new_pool) {
    return -1;
  }

  new_pool->n_elements = n_elements;
  new_pool->element_size = element_size;

  /* Setup memory pool buffer area */
  new_pool->pool = (node_t *)buff;
  if (!new_pool->pool) {
    return -2;
  }

  /* Create linked list out of all the nodes in the pool, adding them to the
   * list of free nodes so they are ready to be allocated. */
  new_pool->free_nodes_head = new_pool->pool;
  node_t *current = NULL;
  node_t *next_free = NULL;
  /* Make linked list from last element to first,
   * very last element points to NULL. */
  for (s32 i=n_elements-1; i>=0; i--) {
    current = get_node_n(new_pool, new_pool->free_nodes_head, i);
    current->hdr.next = next_free;
    next_free = current;
  }

  /* No nodes currently allocated. */
  new_pool->allocated_nodes_head = NULL;

  return 0;
}

/** Destroy a memory pool.
 * Cleans up and frees the memory associated with the pool. This must only be
 * called on memory pools allocated with memory_pool_new().
 *
 * \param pool Pointer to the memory pool to destroy.
 */
void memory_pool_destroy(memory_pool_t *pool)
{
  free(pool->pool);
  free(pool);
}

/** Calculates the number of free (unallocated) elements remaining in the
 * collection.
 * This operation is O(N) in the number of free elements.
 *
 * \param pool Pointer to a memory pool
 * \returns Number of free elements or `< 0` on an error.
 */
s32 memory_pool_n_free(memory_pool_t *pool)
{
  u32 count = 0;

  node_t *p = pool->free_nodes_head;
  while (p && count <= pool->n_elements) {
    p = p->hdr.next;
    count++;
  }

  if (count == pool->n_elements && p)
    /* The list of free elements is larger than the pool,
     * something has gone horribly wrong. */
    return -1;

  return count;
}

/** Calculates the number of elements already allocated in the collection.
 * This operation is O(N) in the number of allocated elements.
 *
 * \param pool Pointer to a memory pool
 * \returns Number of allocated elements or `< 0` on an error.
 */
s32 memory_pool_n_allocated(memory_pool_t *pool)
{
  u32 count = 0;

  node_t *p = pool->allocated_nodes_head;
  while (p && count <= pool->n_elements) {
    p = p->hdr.next;
    count++;
  }

  if (count == pool->n_elements && p)
    /* The list of free elements is larger than the pool,
     * something has gone horribly wrong. */
    return -1;

  return count;
}

/** Basic getter method for memory_pool_t's n_elements field.
 */
u32 memory_pool_n_elements(memory_pool_t *pool)
{
  return pool->n_elements;
}
/** Write all of the elements of a collection to an array.
 * To determine how much space is needed in the destination array you must call
 * memory_pool_n_allocated(). The required space is:
 *
 * ~~~
 * memory_pool_n_allocated() * element_size
 * ~~~
 *
 * \param pool Pointer to a memory pool
 * \param array Array to which the elements will be written
 * \return Number of elements written to the array or `< 0` on an error.
 */
s32 memory_pool_to_array(memory_pool_t *pool, void *array)
{
  u32 count = 0;

  node_t *p = pool->allocated_nodes_head;
  while (p && count <= pool->n_elements) {
    memcpy((u8 *)array + count*pool->element_size, p->elem, pool->element_size);
    p = p->hdr.next;
    count++;
  }

  if (count == pool->n_elements && p)
    /* The list of free elements is larger than the pool,
     * something has gone horribly wrong. */
    return -1;

  return count;
}

/** Adds an element to a collection.
 * Allocates and element from the pool and adds it to the collection of
 * elements, then returns a pointer to the new element.
 *
 * \param pool Pointer to a memory pool
 * \return A pointer to the new element or NULL if the pool is full.
 */
element_t *memory_pool_add(memory_pool_t *pool)
{
  /* Take the head of the list of free nodes, insert it as the head of the
   * allocated nodes and return a pointer to the node's element. */

  if (!pool->free_nodes_head) {
    /* free_nodes_head is NULL, no free nodes available, pool is full. */
    return NULL;
  }

  node_t *new_node = pool->free_nodes_head;
  pool->free_nodes_head = new_node->hdr.next;

  new_node->hdr.next = pool->allocated_nodes_head;
  pool->allocated_nodes_head = new_node;

  return new_node->elem;
}

/** Map a function across all elements allocated in the collection.
 *
 * \param pool Pointer to a memory pool
 * \param f Pointer to a function that does an in-place update of an element.
 * \return Number of elements mapped across or `< 0` on an error.
 */
s32 memory_pool_map(memory_pool_t *pool, void *arg, void (*f)(void *arg, element_t *elem))
{
  u32 count = 0;

  node_t *p = pool->allocated_nodes_head;
  while (p && count <= pool->n_elements) {
    (*f)(arg, p->elem);
    p = p->hdr.next;
    count++;
  }

  if (count == pool->n_elements && p)
    /* The list of elements is larger than the pool,
     * something has gone horribly wrong. */
    return -1;

  return count;
}

/** Calculate a fold reduction on the collection, optionally applying a map at
 * the same time.
 *
 * \param pool Pointer to a memory pool
 * \param x0 Pointer to an initial accumulator state.
 * \param f Pointer to a function that does an in-place update of an
 *          accumulator state given an element and optionally updates that
 *          element in-place.
 * \return Number of elements folded or `< 0` on an error.
 */
s32 memory_pool_fold(memory_pool_t *pool, void *x0,
                     void (*f)(void *x, element_t *elem))
{
  u32 count = 0;

  node_t *p = pool->allocated_nodes_head;
  while (p && count <= pool->n_elements) {
    (*f)(x0, p->elem);
    p = p->hdr.next;
    count++;
  }

  if (count == pool->n_elements && p)
    /* The list of elements is larger than the pool,
     * something has gone horribly wrong. */
    return -1;

  return count;
}

/** Calculate a double valued fold reduction on the collection, optionally
 * applying a map at the same time.
 *
 * \param pool Pointer to a memory pool
 * \param x0 Initial accumulator state.
 * \param f Pointer to a function that returns a new accumulator value given a
 *          current accumulator value and an element, optionally updating the
 *          element in-place.
 * \return Result of the fold operation, i.e. final accumulator value.
 */
double memory_pool_dfold(memory_pool_t *pool, double x0,
                         double (*f)(double x, element_t *elem))
{
  u32 count = 0;
  double x = x0;

  node_t *p = pool->allocated_nodes_head;
  while (p && count <= pool->n_elements) {
    x = (*f)(x, p->elem);
    p = p->hdr.next;
    count++;
  }

  return x;
}

/** Calculate a float valued fold reduction on the collection, optionally
 * applying a map at the same time.
 *
 * \param pool Pointer to a memory pool
 * \param x0 Initial accumulator state.
 * \param f Pointer to a function that returns a new accumulator value given a
 *          current accumulator value and an element, optionally updating the
 *          element in-place.
 * \return Result of the fold operation, i.e. final accumulator value.
 */
float memory_pool_ffold(memory_pool_t *pool, float x0,
                        float (*f)(float x, element_t *elem))
{
  u32 count = 0;
  float x = x0;

  node_t *p = pool->allocated_nodes_head;
  while (p && count <= pool->n_elements) {
    x = (*f)(x, p->elem);
    p = p->hdr.next;
    count++;
  }

  return x;
}

/** Calculate a s32 valued fold reduction on the collection, optionally
 * applying a map at the same time.
 *
 * \param pool Pointer to a memory pool
 * \param x0 Initial accumulator state.
 * \param f Pointer to a function that returns a new accumulator value given a
 *          current accumulator value and an element, optionally updating the
 *          element in-place.
 * \return Result of the fold operation, i.e. final accumulator value.
 */
s32 memory_pool_ifold(memory_pool_t *pool, s32 x0,
                      s32 (*f)(s32 x, element_t *elem))
{
  u32 count = 0;
  s32 x = x0;

  node_t *p = pool->allocated_nodes_head;
  while (p && count <= pool->n_elements) {
    x = (*f)(x, p->elem);
    p = p->hdr.next;
    count++;
  }

  return x;
}

/** Filter elements in the collection, returning filtered out elements back to
 * the pool.
 *
 * \param pool Pointer to a memory pool
 * \param f Pointer to a function that takes an element and returns `0` to
 *          discard that element or `!=0` to keep that element.
 * \return Number of elements in the filtered collection or `< 0` on an error.
 */
s32 memory_pool_filter(memory_pool_t *pool, void *arg, s8 (*f)(void *arg, element_t *elem))
{
  u32 count = 0;

  /* Construct a fake 'previous' node for the head of the list, this eliminates
   * special cases for the beginning of the list. */
  node_t fake_head_node_prev = {
    .hdr = {
      .next = pool->allocated_nodes_head,
    },
  };

  node_t *p_prev = &fake_head_node_prev;
  node_t *p = pool->allocated_nodes_head;

  while (p && count <= pool->n_elements) {
    if ((*f)(arg, p->elem)) {
      /* Keep element, move along.. */
      p_prev = p;
      p = p->hdr.next;
      count++;
    } else {
      /* Drop this element from the list */
      p_prev->hdr.next = p->hdr.next;
      /* and return its node to the pool. */
      p->hdr.next = pool->free_nodes_head;
      pool->free_nodes_head = p;
      p = p_prev->hdr.next;
    }
  }

  /* Use our fake previous node to update the head pointer. */
  pool->allocated_nodes_head = fake_head_node_prev.hdr.next;

  if (count == pool->n_elements && p)
    /* The list of elements is larger than the pool,
     * something has gone horribly wrong. */
    return -1;

  return count;
}

/** Sort the elements in a collection.
 * This is implemented as a merge sort on a linked list and has O(N log N) time
 * complexity and O(1) space complexity. The implementation is stable and has
 * no pathological cases. The worst-case running time is still O(N log N).
 *
 * Ordering is defined by a comparison function `cmp` which takes two elements `a` and `b` and returns `>0` if `a` should be placed before `b`, `<0` if `b` should be placed before `a` ans `0` if they are 'equal' and their current ordering should be preserved.
 *
 * The comparison function also takes an arbitrary argument which is passed
 * through from memory_pool_sort(), this can be used to pass a key or index to
 * sort on to a general comparison function, for example.
 *
 * This implementation is based on the excellent implementation by Simon Tatham:
 *
 * http://www.chiark.greenend.org.uk/~sgtatham/algorithms/listsort.html
 *
 * > Based on code copyright 2001 Simon Tatham.
 * >
 * > Permission is hereby granted, free of charge, to any person
 * > obtaining a copy of this software and associated documentation
 * > files (the "Software"), to deal in the Software without
 * > restriction, including without limitation the rights to use,
 * > copy, modify, merge, publish, distribute, sublicense, and/or
 * > sell copies of the Software, and to permit persons to whom the
 * > Software is furnished to do so, subject to the following
 * > conditions:
 * >
 * > The above copyright notice and this permission notice shall be
 * > included in all copies or substantial portions of the Software.
 * >
 * > THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * > EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * > OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * > NONINFRINGEMENT.  IN NO EVENT SHALL SIMON TATHAM BE LIABLE FOR
 * > ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * > CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * > CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * > SOFTWARE.
 *
 * \param pool Pointer to a memory pool
 * \param arg Arbitrary argument passed through to the comparison function
 * \param cmp Comparison function used to define the sort ordering
 */
void memory_pool_sort(memory_pool_t *pool, void *arg,
                      s32 (*cmp)(void *arg, element_t *a, element_t *b))
{
  /* If collection is empty, return immediately. */
  if (!pool->allocated_nodes_head)
    return;

  u32 insize = 1;

  while (1) {
    node_t *p = pool->allocated_nodes_head;
    pool->allocated_nodes_head = NULL;
    node_t *tail = NULL;

    u32 nmerges = 0;  /* count number of merges we do in this pass */

    while (p) {
      nmerges++;  /* there exists a merge to be done */
      /* step `insize' places along from p */
      node_t *q = p;
      u32 psize = 0;
      for (u32 i = 0; i < insize; i++) {
        psize++;
        q = q->hdr.next;
        if (!q) break;
      }

      /* if q hasn't fallen off end, we have two lists to merge */
      u32 qsize = insize;

      /* now we have two lists; merge them */
      while (psize > 0 || (qsize > 0 && q)) {

        node_t *e;

        /* decide whether next element of merge comes from p or q */
        if (psize == 0) {
          /* p is empty; e must come from q. */
          e = q; q = q->hdr.next; qsize--;
        } else if (qsize == 0 || !q) {
          /* q is empty; e must come from p. */
          e = p; p = p->hdr.next; psize--;
        } else if (cmp(arg, p->elem, q->elem) <= 0) {
          /* First element of p is lower (or same);
           * e must come from p. */
          e = p; p = p->hdr.next; psize--;
        } else {
          /* First element of q is lower; e must come from q. */
          e = q; q = q->hdr.next; qsize--;
        }

        /* add the next element to the merged list */
        if (tail) {
          tail->hdr.next = e;
        } else {
          pool->allocated_nodes_head = e;
        }
        tail = e;
      }

      /* now p has stepped `insize' places along, and q has too */
      p = q;
    }
    tail->hdr.next = NULL;

    /* If we have done only one merge, we're finished. */
    if (nmerges <= 1)   /* allow for nmerges==0, the empty list case */
      return;

    /* Otherwise repeat, merging lists twice the size */
    insize *= 2;
  }
}

/** Perform a groupby type reduction on a collection.
 * A groupby reduction consists of two steps:
 *
 * 1. First the collection is grouped according to some comparison function,
 *    i.e. grouped into sets where the comparison function evaluates to `0` for
 *    any pair of elements in that set.
 *
 *    The definition of the comparison function is the same as for
 *    memory_pool_sort() and is required to define a global ordering on the
 *    whole collection.
 *
 * 2. Next each group is iterated over by an aggregation function which reduces
 *    the whole group to exactly one new element that will be in the collection
 *    after reduction. The aggregation function is similar to a fold function.
 *
 *    The aggregation function is called once for each element of the group and
 *    passed a pointer to the new element that is the aggregate of all elements
 *    in the group. The new element is initialized as being a copy of the first
 *    element of the group but should be updated in-place by the aggregation
 *    function. The aggregation function is also passed a variable `n` which
 *    indicated which element of the group is being passed and the current
 *    group element is passed as parameter `elem`.
 *
 *    The aggregation function also takes an arbitrary argument `x` which is
 *    reset to `x0` at the start of processing each group by copying `x0` into
 *    a working area each time.
 *
 * \param pool Pointer to a memory pool
 * \param arg Arbitrary argument passed through to the comparison function
 * \param cmp Comparison function used to define the grouping
 * \param x0 Arbitrary argument passed to the aggregation function, reset to
 *           this value on each new group.
 * \param x_size The size in bytes of the `x0` argument
 * \param agg The aggregation function
 */
void memory_pool_group_by(memory_pool_t *pool, void *arg,
                          s32 (*cmp)(void *arg, element_t *a, element_t *b),
                          void *x0, size_t x_size,
                          void (*agg)(element_t *new, void *x, u32 n, element_t *elem))
{
  /* If collection is empty, return immediately. */
  if (!pool->allocated_nodes_head)
    return;

  /* First sort the existing list using the compare function. */
  memory_pool_sort(pool, arg, cmp);

  /* Save the head of the unaggregated list and reset the pool head where the
   * aggregated data will be added. */
  node_t *old_head = pool->allocated_nodes_head;
  pool->allocated_nodes_head = NULL;

  /* Allocate working area for the fold function. */
  u8 x_work[x_size];

  u32 count = 0;

  node_t *p = old_head;

  while (p && count <= pool->n_elements) {
    u32 group_count = 0;

    /* Keep a pointer to the head of the group. */
    node_t *group_head = p;

    /* Re-initialize the working area. */
    if (x_size)
      memcpy(x_work, x0, x_size);

    /* Create a new element to hold the aggregate of this group. */
    element_t *new_elem = memory_pool_add(pool);

    /* Initialize the new element to the first element in the group. */
    memcpy(new_elem, p->elem, pool->element_size);

    /* Aggregate this group. */
    do {
      agg(new_elem, (void *)x_work, group_count, p->elem);
      group_count++;

      /* Store pointer to next node to process. */
      node_t *next_p = p->hdr.next;

      /* Return current node to the pool. */
      p->hdr.next = pool->free_nodes_head;
      pool->free_nodes_head = p;

      p = next_p;
    } while(p && cmp(arg, group_head->elem, p->elem) == 0);

    count++;
  }
}

/** Cartesian product of a memory pool collection with an array.
 * For each pair of an element in the original collection and an item in the
 * array `xs`, a new element is created in the updated collection formed by the
 * function `prod` acting on that element together with that array item.
 *
 * The product function takes a pointer to the new element which should be
 * updated in-place, a pointer to an item `x` from the input array `xs` and a
 * pointer to the old element `elem`.
 *
 * As a convenience the product function is also supplied with the total numer
 * of items in the array `n_xs` and the current item number from the array `n`.
 *
 * \param pool Pointer to a memory pool
 * \param xs Array to take the Cartesian product with
 * \param n_xs Number of items in array `xs`
 * \param x_size The size in bytes of each item in `xs`
 * \param prod The product function
 */
s32 memory_pool_product(memory_pool_t *pool, void *xs, u32 n_xs, size_t x_size,
                        void (*prod)(element_t *new, void *x, u32 n_xs, u32 n, element_t *elem))
{
  /* Save the head of the original list and reset the pool head where the
   * product data will be added. */
  node_t *old_head = pool->allocated_nodes_head;
  pool->allocated_nodes_head = NULL;

  u32 count = 0;

  node_t *p = old_head;
  while (p && count <= pool->n_elements) {
    /* Add one element to the new list for each pair of
     * the current element and an item in xs. */
    for (u32 i=0; i<n_xs; i++) {
      element_t *new = memory_pool_add(pool);
      if (!new) {
        /* Pool is full. */
        return -2;
      }
      /* Initialize the element to the same as the original element. */
      memcpy(new, p->elem, pool->element_size);
      prod(new, ((u8 *)xs + i*x_size), n_xs, i, p->elem);
      count++;
    }

    /* Store pointer to next node to process. */
    node_t *next_p = p->hdr.next;

    /* Return current node to the pool. */
    p->hdr.next = pool->free_nodes_head;
    pool->free_nodes_head = p;

    p = next_p;
  }

  if (count == pool->n_elements && p)
    /* The list of elements is larger than the pool,
     * something has gone horribly wrong. */
    return -1;

  return count;
}

s32 memory_pool_product_generator(memory_pool_t *pool, void *x0, u32 max_xs, size_t x_size,
                                  s8 (*next)(void *x, u32 n),
                                  void (*prod)(element_t *new, void *x, u32 n, element_t *elem))
{
  /* Save the head of the original list and reset the pool head where the
   * product data will be added. */
  node_t *old_head = pool->allocated_nodes_head;
  pool->allocated_nodes_head = NULL;

  u32 count = 0;

  node_t *p = old_head;
  while (p && count <= pool->n_elements) {
    /* Iterate through our generator adding a new element for each pair
     * of a generated element and the current element from the old list. */
    u8 x_work[x_size];
    memcpy(x_work, x0, x_size);
    u32 x_count = 0;
    do {
      if (x_count > max_xs) {
        /* Exceded maximum number of generator iterations. */
        return -3;
      }
      element_t *new = memory_pool_add(pool);
      if (!new) {
        /* Pool is full. */
        return -2;
      }
      /* Initialize the element to the same as the original element. */
      memcpy(new, p->elem, pool->element_size);
      prod(new, x_work, x_count, p->elem);
      x_count++;
      count++;
    } while (next(x_work, x_count));

    /* Store pointer to next node to process. */
    node_t *next_p = p->hdr.next;

    /* Return current node to the pool. */
    p->hdr.next = pool->free_nodes_head;
    pool->free_nodes_head = p;

    p = next_p;
  }

  if (count == pool->n_elements && p)
    /* The list of elements is larger than the pool,
     * something has gone horribly wrong. */
    return -1;

  return count;
}

/** \} */


