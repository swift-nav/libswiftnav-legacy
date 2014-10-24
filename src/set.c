#define NDEBUG

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "single_diff.h"

#include "set.h"
#include "iterator.h"

/* Set Utils */
const void *set_index(const set_t *set, u16 index)
{
  return set->arr + set->size * index;
}
const void *set_ending(const set_t *set)
{
  return set_index(set, set->len);
}
/* Set Iterator */
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
  s->current = s->set->arr;
}
const void *set_current(iterator_t *it)
{
  return set_state(it)->current;
}
bool set_end(const iterator_t *it)
{
  const set_t *set = set_state((iterator_t *) it)->set;
  return set_ending(set) == set_current((iterator_t *) it);
}
void mk_set_itr(iterator_t *it, set_state_t *s, set_t *set)
{
  s->current  = set->arr;
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

bool valid_ref(const ptd_set_t *ptd)
{
  const set_t *set = &ptd->set;
  const void *ref = ptd->ref;

  if (is_empty(set)) {
    return ref == NULL;
  }
  // Nonempty set must either have no reference, or must contain reference
  if (ref != NULL && (ref < set->arr ||
                      ref > set_ending(set))) {
    return false;
  }
  return true;
}

/* Tests for valid set structure */
bool is_set(const set_t *set)
{
  key prev, next;
  if (is_empty(set)) {
    return true;
  }
  if (set->len < 0) {
    return false;
  }
  /* Set must have strictly increasing sequence of keys. */
  prev = set->keyfn(set->arr);
  for (u8 i = 1; i < set->len; i++) {
    next = set->keyfn(set_index(set, i));
    if (next <= prev) {
      return false;
    }
    prev = next;
  }
  return true;
}
/* Tests for valid pointed set structure */
bool is_ptd_set(const ptd_set_t *ptd)
{
  const void *ref = ptd->ref;
  const set_t *set = &ptd->set;

  // Empty set is valid as long as ref is null
  if (is_empty(set)) {
    return ref == NULL;
  }

  if (!valid_ref(ptd)) {
    return false;
  }
  return is_set(&ptd->set);
}

void mk_set(set_t *set, u16 len, size_t size, void *arr, key_fn_t *keyfn)
{
  set->len = len;
  set->size = size;
  set->arr = arr;
  set->keyfn = keyfn;
  assert(is_set(set));
}
void mk_ptd_set(ptd_set_t *ptd, u16 len, size_t size, void *arr, key ref, key_fn_t *keyfn)
{
  mk_set(&ptd->set, len, size, arr, keyfn);
  set_ref(ptd, ref);
}
void mk_sdiff_set(set_t *set, int len, sdiff_t *arr)
{
  mk_set(set, len, sizeof(sdiff_t), arr, &sdiff_key);
}
void mk_prn_set(set_t *set, int len, prn *arr)
{
  mk_set(set, len, sizeof(prn), arr, &prn_key);
}
void freeze_set(set_t *set, iterator_t *it, u16 max_len, size_t size, key_fn_t *keyfn)
{
  s8 len = freeze_arr(size, max_len, set->arr, it);
  mk_set(set, len, size, set->arr, keyfn);
}
void freeze_ptd(ptd_set_t *ptd, iterator_t *it, u16 max_len, size_t size, key ref, key_fn_t *keyfn)
{
  freeze_set(&ptd->set, it, max_len, size, keyfn);
  set_ref(ptd, ref);
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
const void *find_key(iterator_t *it, key key, key_fn_t *keyfn)
{
  for(reset(it); more(it); next(it)) {
    if (keyfn(current(it)) == key) {
      return current(it);
    }
  }
  return NULL;
}
void set_ref(ptd_set_t *ptd, key key)
{
  iterator_t it;
  set_state_t s;
  mk_set_itr(&it, &s, &ptd->set);
  const void *ref = find_key(&it, key, ptd->set.keyfn);
  assert(ref != NULL);
  ptd->ref = ref;
}
key get_ref(ptd_set_t *ptd)
{
  return ptd->set.keyfn(ptd->ref);
}
void mk_without_ref_itr(iterator_t *it, filter_state_t *s, iterator_t *base, ptd_set_t *ptd)
{
  mk_filter_itr(it, s, &not_eq_ptr, ptd->ref, base);
}

void ref_fst_to_ptd(ptd_set_t *ptd, u16 len, size_t size, void *arr)
{
  key ref = ptd->set.keyfn(&arr[0]);
  u8 i = 1, j = 0;
  bool set_ref = false;
  while (i < len && j < len) {
    if (ptd->set.keyfn(&arr[i]) > ref && !set_ref) {
      memcpy(ptd->set.arr + j * size, arr, size);
      set_ref = true;
      j++;
    } else {
      memcpy(ptd->set.arr + j * size, arr + i * size, size);
      j++, i++;
    }
  }
  if (j < len) {
    memcpy(ptd->set.arr + j * size, arr, size);
  }

  mk_ptd_set(ptd, len, size, ptd->set.arr, ref, ptd->set.keyfn);
}
void ptd_to_ref_fst(ptd_set_t *ptd, void *arr)
{
  /* Copy ref */
  size_t size = ptd->set.size;
  u16 len = ptd->set.len;
  memcpy(arr, ptd->ref, size);

  /* Freeze rest */
  iterator_t set_itr;
  set_state_t ss;
  iterator_t fitr;
  filter_state_t fs;
  mk_set_itr(&set_itr, &ss, &ptd->set);
  mk_without_ref_itr(&fitr, &fs, &set_itr, ptd);
  freeze_arr(size, len-1, arr+size, &fitr);
}
