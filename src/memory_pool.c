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

struct node;

typedef struct {
  struct node *next;
} memory_pool_node_hdr_t;

typedef struct node {
  memory_pool_node_hdr_t hdr;
  /* C99 "flexible member array" (see C99 std. ch. 6.7.2.1),
   * Allows us to get a pointer to the top of the unknown size element. */
  element_t elem[];
} node_t;

struct _memory_pool {
  u32 n_elements;
  size_t element_size;
  node_t *pool;
  node_t *free_nodes_head;
  node_t *allocated_nodes_head;
};

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
    return 0;
  }

  new_pool->n_elements = n_elements;
  new_pool->element_size = element_size;

  /* Allocate memory pool */
  size_t node_size = calc_node_size(element_size);
  new_pool->pool = (node_t *)malloc(node_size * n_elements);
  if (!new_pool->pool) {
    free(new_pool);
    return 0;
  }

  /* Create linked list out of all the nodes in the pool, adding them to the
   * list of free nodes so they are ready to be allocated. */
  new_pool->free_nodes_head = new_pool->pool;
  node_t *current = 0;
  node_t *next_free = 0;
  /* Make linked list from last element to first,
   * very last element points to NULL. */
  for (s32 i=n_elements-1; i>=0; i--) {
    current = get_node_n(new_pool, new_pool->free_nodes_head, i);
    current->hdr.next = next_free;
    next_free = current;
  }

  /* No nodes currently allocated. */
  new_pool->allocated_nodes_head = 0;

  return new_pool;
}

/** Destroy a memory pool.
 * Cleans up and frees the memory associated with the pool.
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
    return 0;
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
s32 memory_pool_map(memory_pool_t *pool, void (*f)(element_t *elem))
{
  u32 count = 0;

  node_t *p = pool->allocated_nodes_head;
  while (p && count <= pool->n_elements) {
    (*f)(p->elem);
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
s32 memory_pool_filter(memory_pool_t *pool, s8 (*f)(element_t *elem))
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
    if ((*f)(p->elem)) {
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


/** \} */


