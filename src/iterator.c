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
  return 0;
}

// TODO move intersection methods
filter_state_t *filter_state(iterator_t *it)
{
  return it->state;
}
intersection_state_t *intersection_state(iterator_t *it)
{
  return it->state;
}
bool filter_end(const iterator_t *it)
{
  return end((const iterator_t *)filter_state((iterator_t *)it)->it);
}
bool intersection_end(const iterator_t *it)
{
  const intersection_state_t *s = intersection_state((iterator_t *)it);
  return end((const iterator_t *)s->it1) || end((const iterator_t *)s->it2);
}
void filter_next0(iterator_t *it)
{
  filter_state_t *s = filter_state(it);
  while(more(s->it) && !s->f(s->arg, current(s->it))) {
    next(s->it);
  }
}
void intersection_next0(iterator_t *it)
{
  intersection_state_t *s = intersection_state(it);
  key key1, key2;
  while (more(s->it1) && more(s->it2)) {
    key1 = s->key1(current(s->it1));
    key2 = s->key2(current(s->it2));

    if (key1 == key2) {
      // Update current ptrs
      s->current.fst = current(s->it1);
      s->current.snd = current(s->it2);
      return;
    } else if (key1 < key2) {
      next(s->it1);
    } else if (key1 > key2) {
      next(s->it2);
    }
  }
}
void filter_next1(iterator_t *it)
{
  filter_state_t *s = filter_state(it);
  next(s->it);
  filter_next0(it);
}
void intersection_next1(iterator_t *it)
{
  intersection_state_t *s = intersection_state(it);
  next(s->it1);
  next(s->it2);
  intersection_next0(it);
}
void filter_reset(iterator_t *it)
{
  reset(filter_state(it)->it);
  filter_next0(it);
}
void intersection_reset(iterator_t *it)
{
  intersection_state_t *s = intersection_state(it);
  reset(s->it1);
  reset(s->it2);
  intersection_next0(it);
}
const void *filter_current(iterator_t *it)
{
  return current(filter_state(it)->it);
}
const void *intersection_current(iterator_t *it)
{
  intersection_state_t *s = intersection_state(it);
  return &s->current;
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
void mk_intersection_itr(iterator_t *it, intersection_state_t* s,
                         iterator_t *it1, iterator_t *it2,
                         key (*key1)(const void *),
                         key (*key2)(const void *))
{
  s->it1 = it1;
  s->it2 = it2;
  s->key1 = key1;
  s->key2 = key2;
  
  it->state = s;
  it->end = &intersection_end;
  it->next = &intersection_next1;
  it->reset = &intersection_reset;
  it->current = &intersection_current;

  reset(it);
}

// Checks to see if iterators iterate over the same memory
bool ptr_itr_equality(iterator_t *it1, iterator_t *it2)
{
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
void each(iterator_t *it, void (*f)(const void *, const void *), const void *arg)
{
  for(reset(it); more(it); next(it)) {
    f(arg, current(it));
  }
}

/* Filter utils */
