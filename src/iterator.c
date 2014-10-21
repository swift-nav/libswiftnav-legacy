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
