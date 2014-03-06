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

/** \defgroup memory_pool Memory Pool
 * Simple fixed size memory pool supporting map and fold operations.
 * \{ */

inline static size_t calc_node_size(size_t element_size)
{
  return element_size + sizeof(memory_pool_node_hdr_t);
}

inline static node_t *get_node_n(memory_pool_t *pool, node_t *head, u32 n)
{
  return (node_t *)((u8 *)head + calc_node_size(pool->element_size) * n);
}

/* The pool consists of n_elements nodes, each of which contains a node header
 * for internal bookkeeping and an 'element' which is the user defined payload.
 * */

memory_pool_t *new_memory_pool(u32 n_elements, size_t element_size)
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

void destroy_memory_pool(memory_pool_t *pool)
{
  free(pool->pool);
  free(pool);
}

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

element_t *memory_pool_append(memory_pool_t *pool)
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


