#define NDEBUG

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "single_diff.h"

#include "set.h"
#include "iterator.h"

// Iterator: full set
set_state_t *set_state(iterator_t *it)
{
  return it->state;
}
void set_next(iterator_t *it)
{
  set_state_t *s = set_state(it);
  s->current += s->set->size;
}
void set_reset(iterator_t *it)
{
  set_state_t *s = set_state(it);
  s->current = s->set->set;
}
const void *set_current(iterator_t *it)
{
  return set_state(it)->current;
}
const void *set_ending(const set_t *set)
{
  return set->set + set->len * set->size;
}
bool set_end(const iterator_t *it)
{
  const set_t *set = set_state((iterator_t *) it)->set;
  return set_ending(set) == set_current((iterator_t *) it);
}
void mk_set_itr(iterator_t *it, set_state_t *s, set_t *set)
{
  s->current  = set->set;
  s->set      = set;

  it->state   = s;
  it->end     = &set_end;
  it->next    = &set_next;
  it->reset   = &set_reset;
  it->current = &set_current;
}

// DONE:
// - freeze, filter, intersect, is_subset, fold
// - set_ref_prn, without_ref_itr
// - print
//
// TODO:
//
// general:
//   map
// needs key:
//   contains?, 


// Iterator: set intersection 
// ...

// Unit tests
// integrate one function at a time
//
//
// Set operations?
// Delete prn
// Add prn

bool is_empty(const set_t *set) {
  return set->len == 0;
}

bool valid_ref(const ptd_set_t *set)
{
  const void *ref = set->ref;

  if (is_empty(set)) {
    return ref == NULL;
  }
  // Nonempty set must either have no reference, or must contain reference
  if (ref != NULL && (ref < set->set ||
                      ref > set_ending(set))) {
    return false;
  }
  return true;
}

// Tests for valid set structure
bool is_set(const set_t *set, key (*to_key)(const void *))
{
  key prev, current;
  const void *ref = set->ref;

  // Empty set is valid as long as ref is null
  if (is_empty(set)) {
    return ref == NULL;
  }

  if (!valid_ref(set)) {
    return false;
  }

  // Set must have strictly increasing sequence of keys
  prev = to_key(set->set);
  for (u8 i = 1; i < set->len; i++) {
    current = to_key(set->set + i * set->size);
    if (current <= prev) {
      return false;
    }
    prev = current;
  }
  return true;
}

void mk_sdiff_set(set_t *set, int len, sdiff_t *sats)
{
  mk_pointed(set, len, sizeof(sdiff_t), sats);
  assert(is_set(set, &sdiff_key));
}
void mk_prn_set(set_t *set, int len, prn *arr)
{
  mk_pointed(set, len, sizeof(prn), arr);
  assert(is_set(set, &prn_key));
}
void mk_pointed(set_t *set, int len, size_t size, const void *arr)
{
  set->len = len;
  set->size = size;
  set->set = arr;
  set->ref = NULL;
}

void freeze_to_set(size_t max_len, set_t *set, iterator_t *it)
{
  s8 len = freeze_itr(set->size, max_len, set->set, it);
  set->len = len;
  assert(len >= 0);
}
void print_prn(const void *arg, const void *elem)
{
  (void)arg;
  printf("prn: %i\n", *(prn *)elem);
}
void print_prn_tuple(const void *arg, const void *elem)
{
  (void)arg;
  tuple *t = (tuple *)elem;
  printf("prn1: %i\n", *((prn *)t->fst));
  printf("prn2: %i\n", *((prn *)t->snd));
}

bool eq_ref(const void *arg, const void *elem)
{
  return *(prn *)elem == *(prn *)arg;
}
bool not_eq_ref(const void *arg, const void *elem)
{
  return !(*((prn *)elem) == *((prn *)arg));
}
bool not_eq_ptr(const void *p1, const void *p2)
{
  return p1 != p2;
}

/* PRN Utilities */
key prn_key(const void *x)
{
  return (key)*((prn *) x);
}
key sdiff_key(const void *x)
{
  return (key)((sdiff_t *) x)->prn;
}
const void *find_key(iterator_t *it, key key, key (*f)(const void*))
{
  for(reset(it); more(it); next(it)) {
    if (f(current(it)) == key) {
      return current(it);
    }
  }
  return NULL;
}
void set_ref_prn(ptd_set_t *ptd, prn prn, key (*f)(const void *))
{
  iterator_t it;
  set_state_t s;
  mk_set_itr(&it, &s, ptd->set);
  const void *ref = find_key(&it, prn);
  assert(ref != NULL)
  set->ref = ref;
}
s8 set_
void mk_without_ref_itr(iterator_t *it, filter_state_t *s, iterator_t *base, set_t *set)
{
  mk_filter_itr(it, s, &not_eq_ptr, set->ref, base);
}
