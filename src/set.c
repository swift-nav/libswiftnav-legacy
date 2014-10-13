#define NDEBUG

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "single_diff.h"
#include "linear_algebra.h"
#include "set.h"

void next(iterator_t *it)
{
  it->next(it);
}
bool end(const iterator_t *it)
{
  return it->end(it);
}
bool more(const iterator_t *it)
{
  return !end(it);
}
void *current(iterator_t *it)
{
  return it->current(it);
}

typedef struct {
  void *current;
  set_t *set;
} set_state_t;

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
void *set_current(iterator_t *it)
{
  return set_state(it)->current;
}
void *set_ending(const set_t *set)
{
  return set->set + set->len * set->size;
}
bool set_end(const iterator_t *it)
{
  set_t *set = set_state((iterator_t *) it)->set;
  return set_ending(set) == set_current((iterator_t *) it);
}
void mk_set_itr(iterator_t *it, set_state_t *s, set_t *set)
{
  s->current  = set->set;
  s->set      = set;

  it->state   = s;
  it->end     = &set_end;
  it->next    = &set_next;
  it->current = &set_current;
}

key set_key(iterator_t* it)
{
  return set_state(it)->set->key(set_current(it));
}


// TODO
// Iterator: set without reference
// ...

// Iterator: set intersection 
// ...

// Unit tests
// integrate one function at a time
//
//
// Set operations?
// Delete prn
// Add prn

// Tests for valid set structure
bool is_set(const set_t *set)
{
  key prev, current;
  void *ref = set->ref;
  // Empty set is valid as long as ref is null
  if (set->len == 0) {
    return ref == NULL;
  }
  // Nonempty set must either have no reference, or must contain reference
  if (ref != NULL && (ref < set->set ||
                      ref > set_ending(set))) {
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

// TODO:
void set_ref_prn(set_t *set, prn prn)
{
}
void mk_sdiff_set(set_t *set, int len, sdiff_t *arr)
{
}
