#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "iterator.h"

bool end(const  iterator_t *it       ) { return it->end(it    ) ; } 
void next(      iterator_t *it       ) { it->next(it          ) ; } 
void reset(     iterator_t *it       ) { it->reset(it         ) ; } 
const void *current(  iterator_t *it ) { return it->current(it) ; } 
bool more(const iterator_t *it       ) { return !end(it       ) ; } 

s8 freeze_itr(size_t elem_size, size_t max_len, void *buffer, iterator_t *it)
{
  size_t i = 0;
  for (reset(it); more(it); next(it), i++) {
    if (i >= max_len) {
      return -1;
    }
    memcpy(buffer + i*elem_size, current(it), elem_size);
  }
  // Return length of iterator
  return i;
}

/* Filter */
filter_state_t *filter_state(iterator_t *it)
{
  return it->state;
}
bool filter_end(const iterator_t *it)
{
  return end((const iterator_t *)filter_state((iterator_t *)it)->it);
}
void filter_next0(iterator_t *it)
{
  filter_state_t *s = filter_state(it);
  while(more(s->it) && !s->f(s->arg, current(s->it))) {
    next(s->it);
  }
}
void filter_next1(iterator_t *it)
{
  filter_state_t *s = filter_state(it);
  next(s->it);
  filter_next0(it);
}
void filter_reset(iterator_t *it)
{
  reset(filter_state(it)->it);
  filter_next0(it);
}
const void *filter_current(iterator_t *it)
{
  return current(filter_state(it)->it);
}
void mk_filter_itr(iterator_t *it, filter_state_t* s, bool (*f)(const void*, const void*),
                   const void *arg, iterator_t *base)
{
  s->it = base;
  s->arg = arg;
  s->f = f;

  it->state = s;
  it->end = &filter_end;
  it->next = &filter_next1;
  it->reset = &filter_reset;
  it->current = &filter_current;

  reset(it);
}

/* Intersection */
intersection_state_t *intersection_state(iterator_t *it)
{
  return it->state;
}
bool intersection_end(const iterator_t *it)
{
  const intersection_state_t *s = intersection_state((iterator_t *)it);
  return end((const iterator_t *)s->it1) || end((const iterator_t *)s->it2);
}
void intersection_next0(iterator_t *it)
{
  intersection_state_t *s = intersection_state(it);
  key key1, key2;
  while (more(s->it1) && more(s->it2)) {
    key1 = s->key1(current(s->it1));
    key2 = s->key2(current(s->it2));

    if (key1 == key2) {
      /* Update current ptr */
      s->map(s->arg, current(s->it1), current(s->it2), s->current);
      return;
    } else if (key1 < key2) {
      next(s->it1);
    } else if (key1 > key2) {
      next(s->it2);
    }
  }
}
void intersection_next1(iterator_t *it)
{
  intersection_state_t *s = intersection_state(it);
  next(s->it1);
  next(s->it2);
  intersection_next0(it);
}
void intersection_reset(iterator_t *it)
{
  intersection_state_t *s = intersection_state(it);
  reset(s->it1);
  reset(s->it2);
  intersection_next0(it);
}
const void *intersection_current(iterator_t *it)
{
  intersection_state_t *s = intersection_state(it);
  return s->current;
}
void mk_intersection_itr(iterator_t *it, intersection_state_t* s, void *current,
                         iterator_t *it1, iterator_t *it2,
                         key (*key1)(const void *),
                         key (*key2)(const void *),
                         void (*map)(const void *, const void *, const void *, void *),
                         const void *arg)
{
  s->it1 = it1;
  s->it2 = it2;
  s->key1 = key1;
  s->key2 = key2;
  s->map = map;
  s->arg = arg;
  s->current = current;
  
  it->state = s;
  it->end = &intersection_end;
  it->next = &intersection_next1;
  it->reset = &intersection_reset;
  it->current = &intersection_current;

  reset(it);
}

/* Checks to see if iterators iterate over the same memory */
bool ptr_itr_equality(iterator_t *it1, iterator_t *it2)
{
  // TODO remove reset?
  for(reset(it1), reset(it2); more(it1) && more(it2); next(it1), next(it2)) {
    if (current(it1) != current(it2)) {
      return false;
    }
  }
  if (more(it1) || more(it2)) {
    return false;
  }
  return true;
}

/* Checks that the keys of it1 are a subset of the keys of it2 */
bool is_subset(iterator_t *it1, iterator_t *it2,
               key (*key1)(const void *),
               key (*key2)(const void *))
{
  reset(it1), reset(it2);
  key k1, k2;
  while(more(it1) && more(it2)) {
    k1 = key1(current(it1));
    k2 = key2(current(it2));
    if (k1 == k2) {
      next(it1), next(it2);
    } else if (k2 < k1) {
      next(it2);
    } else {
      /* k1 is less than k2, so key1 cannot appear in it2 */
      return false;
    }
  }
  if (more(it1)) {
    return false;
  }
  return true;
}

void each(iterator_t *it, void (*f)(const void *, const void *), const void *arg)
{
  for(reset(it); more(it); next(it)) {
    f(arg, current(it));
  }
}

void fold(void (*f)(const void *, void *, const void *),
          const void *arg, void *init, iterator_t *it)
{
  for(reset(it); more(it); next(it)) {
    f(arg, init, current(it));
  }
}

void fst(const void *arg, const void *first, const void *second, void *current)
{
  (void) second;
  memcpy(current, first, *((size_t *) arg));
}

void snd(const void *arg, const void *first, const void *second, void *current)
{
  (void) first;
  memcpy(current, second, *((size_t *) arg));
}
