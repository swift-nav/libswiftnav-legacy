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

key current_key(iterator_t* it)
{
  return set_state(it)->set->key(set_current(it));
}

// TODO
// convert array into itr, itr -> array
//
// general:
//   fold, map, filter - done: freeze, 
// needs key:
//   intersect, subset?, contains?, 
// needs pointed-set:
//   update ref, non-refs (filter)
// needs specific structure:
//   print,


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

bool valid_ref(const set_t *set)
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
bool is_set(const set_t *set)
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
  prev = set->key(set->set);
  for (u8 i = 1; i < set->len; i++) {
    current = set->key(set->set + i * set->size);
    if (current <= prev) {
      return false;
    }
    prev = current;
  }
  return true;
}

key prn_key(const void *x) {
  return *((u8 *) x);
}
void mk_prn_set(set_t *set, int len, prn *arr)
{
  set->len = len;
  set->size = sizeof(prn);
  set->set = arr;
  set->key = &prn_key;
  set->ref = NULL;

  assert(is_set(set));
}

void print_prn(const void *arg, const void *elem)
{
  (void)arg;
  printf("prn: %i\n", *(prn *)elem);
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

const void *find_key(iterator_t *it, key key)
{
  for(reset(it); more(it); next(it)) {
    if (current_key(it) == key) {
      return current(it);
    }
  }
  return NULL;
}

s8 set_ref_prn(set_t *set, prn prn)
{
  iterator_t it;
  set_state_t s;
  mk_set_itr(&it, &s, set);
  const void *ref = find_key(&it, prn);
  if (ref == NULL) {
    return -1;
  }
  set->ref = ref;
  // TODO add assert()
  return 0;
}
void mk_without_ref_itr(iterator_t *it, filter_state_t *s, iterator_t *base, set_t *set)
{
  mk_filter_itr(it, s, &not_eq_ptr, set->ref, base);
}

// TODO:
//void mk_sdiff_set(set_t *set, int len, sdiff_t *arr)
//{
//}
