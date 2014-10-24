/*
 * Copyright (C) 2014 Swift Navigation Inc.
 * Contact: Ian Horn <ian@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "single_diff.h"
#include "sats_management.h"
#include "linear_algebra.h"

#include "set.h"
#include "iterator.h"
#include "iterator_utils.h"

#define DEBUG_SATS_MAN 0

/* ORDER OF OPERATIONS
 *
 * old
 * new
 * new-box
 * test (use original name)
 *
 */


/* choose_reference_sat */
u8 old_choose_reference_sat(u8 num_sats, const sdiff_t *sats)
{
  double best_snr=sats[0].snr;
  u8 best_prn=sats[0].prn;
  for (u8 i=1; i<num_sats; i++) {
    if (sats[i].snr > best_snr) {
      best_snr = sats[i].snr;
      best_prn = sats[i].prn;
    }
  }
  return best_prn;
}
typedef struct {
  prn prn;
  double snr;
} prn_snr;
void max_snr(const void *arg, void *max, const void *elem)
{
  (void) arg;
  prn_snr *t = max;
  const sdiff_t *sdiff = elem;
  if (sdiff->snr > t->snr) {
    t->snr = sdiff->snr;
    t->prn = sdiff->prn;
  }
}
u8 new_choose_reference_sat(iterator_t *sdiffs)
{
  prn_snr best;
  best.snr = -999;
  best.prn = -1;
  fold(&max_snr, NULL, &best, sdiffs);

  return best.prn;
}
u8 box_choose_reference_sat(u8 num_sats, const sdiff_t *sats)
{
  set_t set;
  mk_sdiff_set(&set, num_sats, (sdiff_t *)sats);
  iterator_t it;
  set_state_t s;
  mk_set_itr(&it, &s, &set);

  return new_choose_reference_sat(&it);
}
u8 choose_reference_sat(u8 num_sats, const sdiff_t *sats)
{
  u8 r1 = old_choose_reference_sat(num_sats, sats);
  u8 r2 = box_choose_reference_sat(num_sats, sats);
  assert(r1 == r2);
  return r1;
}
/* END choose_reference_sat */

u8 intersect_sats(u8 num_sats1, u8 num_sdiffs, u8 *sats1, sdiff_t *sdiffs,
                  sdiff_t *intersection_sats)
{
  u8 i, j, n = 0;
  if(DEBUG_SATS_MAN) {
    printf("<INTERSECT_SATS>\n");
    printf("sdiff prns= {");
    for (u8 k=0; k< num_sdiffs; k++) {
      printf("%u, ", sdiffs[k].prn);
    }
    printf("}\n");
    printf("sats= {");
    for (u8 k=0; k< num_sats1; k++) {
      printf("%u, ", sats1[k]);
    }
    printf("}\n");
    printf("(n, i, prn1[i],\t j, prn2[j])\n");
  }
  /* Loop over sats1 and sdffs and check if a PRN is present in both. */
  for (i=0, j=0; i<num_sats1 && j<num_sdiffs; i++, j++) {
    if (sats1[i] < sdiffs[j].prn) {
      if (DEBUG_SATS_MAN) {
        printf("(%u, %u, prn1=%u,\t %u, prn2=%u)\t\t prn1 < prn2  i++\n", n, i, sats1[i], j, sdiffs[j].prn);
      }
      j--;
    }
    else if (sats1[i] > sdiffs[j].prn) {
      if (DEBUG_SATS_MAN) {
        printf("(%u, %u, prn1=%u,\t %u, prn2=%u)\t\t prn1 > prn2  j++\n", n, i, sats1[i], j, sdiffs[j].prn);
      }
      i--;
    }
    else {
      if (DEBUG_SATS_MAN) {
        printf("(%u, %u, prn1=%u,\t %u, prn2=%u)\t\t prn1 = prn2  i,j,n++\n", n, i, sats1[i], j, sdiffs[j].prn);
      }
      memcpy(&intersection_sats[n], &sdiffs[j], sizeof(sdiff_t));
      n++;
    }
  }
  if (DEBUG_SATS_MAN) {
    printf("intersection_sats= {");
    for (u8 k=0; k< n; k++) {
      printf("%u, ", intersection_sats[k].prn);
    }
    printf("}\n</INTERSECT_SATS>\n");
  }
  return n;
}

// TODO replace with pointer change
/** Puts sdiffs into sdiffs_with_ref_first with the sdiff for ref_prn first
 */
void set_reference_sat_of_prns(u8 ref_prn, u8 num_sats, u8 *prns)
{
  u8 old_ref = prns[0];
  u8 j;
  if (old_ref != ref_prn) {
    j = 1;
    u8 old_prns[num_sats];
    memcpy(old_prns, prns, num_sats * sizeof(u8));
    u8 set_old_yet = 0;
    prns[0] = ref_prn;
    for (u8 i=1; i<num_sats; i++) {
      if (old_prns[i] != ref_prn) {
        if (old_prns[i]>old_ref && set_old_yet == 0) {
          prns[j] = old_ref;
          j++;
          i--;
          set_old_yet = 1;
        }
        else {
          prns[j] = old_prns[i];
          j++;
        }
      }
      else if (i == num_sats-1) {
        prns[j] = old_ref;
      }
    }
  }
}


/** Puts sdiffs into sdiffs_with_ref_first with the sdiff for ref_prn first, while updating sats_management
 */
// TODO split function
//   - sats_management operation -> update ref
//   - sdiffs -> update ref
void set_reference_sat(u8 ref_prn, sats_management_t *sats_management,
                          u8 num_sdiffs, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{  
  u8 old_ref = sats_management->prns[0];
  if (DEBUG_SATS_MAN) {
    printf("<SET_REFERENCE_SAT>\n");
    printf("ref_prn = %u\n", ref_prn);
    printf("old_ref = %u\n", old_ref);
  }
  u8 j;
  if (old_ref != ref_prn) {
    j = 1;
    u8 old_prns[sats_management->num_sats];
    memcpy(old_prns, sats_management->prns, sats_management->num_sats * sizeof(u8));
    u8 set_old_yet = 0;
    sats_management->prns[0] = ref_prn;
    for (u8 i=1; i<sats_management->num_sats; i++) {
      if (old_prns[i] != ref_prn) {
        if (old_prns[i]>old_ref && set_old_yet == 0) {
          sats_management->prns[j] = old_ref;
          j++;
          i--;
          set_old_yet = 1;
        }
        else {
          sats_management->prns[j] = old_prns[i];
          j++;
        }
      }
      else if (i == sats_management->num_sats-1) {
        sats_management->prns[j] = old_ref;
      }
    }
  }
  j=1;
  for (u8 i=0; i<num_sdiffs; i++) {
    if (sdiffs[i].prn != ref_prn) {
      if (DEBUG_SATS_MAN) {
        printf("prn[%u] = %u\n", j, sdiffs[i].prn);
      }
      memcpy(&sdiffs_with_ref_first[j], &sdiffs[i], sizeof(sdiff_t));
      j++;
    } else {
      if (DEBUG_SATS_MAN) {
        printf("prn[0] = %u\n", sdiffs[i].prn);
      }
      memcpy(&sdiffs_with_ref_first[0], &sdiffs[i], sizeof(sdiff_t));
    }
  }
  if (DEBUG_SATS_MAN) {
    printf("</SET_REFERENCE_SAT>\n");
  }
}

/* set_reference_sat_and_prns */
// TODO is above function needed?
void old_set_reference_sat_and_prns(u8 ref_prn, sats_management_t *sats_management,
                                    u8 num_sdiffs, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  sats_management->num_sats = num_sdiffs;
  sats_management->prns[0] = ref_prn;
  u8 j=1;
  for (u8 i=0; i<num_sdiffs; i++) {
    if (sdiffs[i].prn != ref_prn) {
      sats_management->prns[j] = sdiffs[i].prn;
      if (sdiffs_with_ref_first) {
        memcpy(&sdiffs_with_ref_first[j], &sdiffs[i], sizeof(sdiff_t));
      }
      j++;
    } else {
      if (sdiffs_with_ref_first) {
        memcpy(&sdiffs_with_ref_first[0], &sdiffs[i], sizeof(sdiff_t));
      }
    }
  }
}

/* each(some_itr) {
 *   // foo bar
 * }
 * */

void new_set_reference_sat_and_prns(u8 ref_prn, iterator_t *sdiffs,
                                    ptd_set_t *sats_ptd, ptd_set_t *sdiffs_ptd)
{
  iterator_t map_itr;
  map_state_t map_state;
  prn current;
  mk_map_itr(&map_itr, &map_state, &current, sdiffs, &map_sdiff_prn, NULL);
  // TODO store max size of buffer in set_t
  freeze_ptd(sats_ptd, &map_itr, MAX_CHANNELS, sizeof(prn), (key)ref_prn, &prn_key);
  freeze_ptd(sdiffs_ptd, sdiffs, MAX_CHANNELS, sizeof(sdiff_t), (key)ref_prn, &sdiff_key);
}

void box_set_reference_sat_and_prns(u8 ref_prn, sats_management_t *sats_management,
                                    u8 num_sdiffs, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  set_t sdiffs_set;
  set_state_t sdiffs_state;
  iterator_t sdiffs_itr;
  mk_set(&sdiffs_set, num_sdiffs, sizeof(sdiff_t), sdiffs, &sdiff_key);
  mk_set_itr(&sdiffs_itr, &sdiffs_state, &sdiffs_set);

  u16 num_sats = sats_management->num_sats;
  prn prn_buff[num_sats];
  ptd_set_t sats_ptd;
  sats_ptd.set.arr = prn_buff;

  ptd_set_t sdiffs_ptd;
  sdiff_t sdiff_buff[num_sdiffs];
  sdiffs_ptd.set.arr = sdiff_buff;

  new_set_reference_sat_and_prns(ref_prn, &sdiffs_itr, &sats_ptd, &sdiffs_ptd);

  sats_management->num_sats = sats_ptd.set.len;
  ptd_to_ref_fst(&sats_ptd, sats_management->prns);
  ptd_to_ref_fst(&sdiffs_ptd, sdiffs_with_ref_first);

  //memcpy(prn_buff, sats_management->prns, num_sats * sizeof(prn));
  //ref_fst_to_ptd(&sats_ptd, num_sats, sizeof(prn),
  //               prn_buff, prn_buff[0]);
  //ptd.set.len = 

  //new_set_reference_sat_and_prns(ref_prn, &sats, &sdiffs_itr, &sdiffs_ptd);

  //ptd_to_ref_fst(&sdiffs_ptd, num_sdiffs, sizeof(sdiff_t), 
  //    sdiffs_with_ref_first, &sdiffs_ptd);
}

void set_reference_sat_and_prns(u8 ref_prn, sats_management_t *sats_management,
                                u8 num_sdiffs, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  sats_management_t man1, man2;
  sdiff_t ref_first1[MAX_CHANNELS];
  sdiff_t ref_first2[MAX_CHANNELS];
  old_set_reference_sat_and_prns(ref_prn, &man1, num_sdiffs, sdiffs, ref_first1);
  box_set_reference_sat_and_prns(ref_prn, &man2, num_sdiffs, sdiffs, ref_first2);
  assert(man1.num_sats == man2.num_sats);
  assert(man1.num_sats == num_sdiffs);
  assert(0 == memcmp(man1.prns, man2.prns, man1.num_sats * sizeof(prn)));
  assert(0 == memcmp(ref_first1, ref_first2, num_sdiffs * sizeof(sdiff_t)));

  memcpy(sats_management, &man1, sizeof(sats_management_t));
  memcpy(sdiffs_with_ref_first, ref_first1, num_sdiffs * sizeof(sdiff_t));
}
/* END set_reference_sat_and_prns */

/* init_sats_management */
// TODO get rid of 0 check?
void old_init_sats_management(sats_management_t *sats_management,
                              u8 num_sdiffs, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  if (DEBUG_SATS_MAN) {
    printf("<INIT_SATS_MANAGEMENT>\n");
  }
  if (num_sdiffs == 0) {
    sats_management->num_sats = 0;
    if (DEBUG_SATS_MAN) {
      printf("</INIT_SATS_MANAGEMENT>\n");
    }
    return;
  }

  u8 ref_prn = choose_reference_sat(num_sdiffs, sdiffs);
  set_reference_sat_and_prns(ref_prn, sats_management,
                             num_sdiffs, sdiffs, sdiffs_with_ref_first);
  if (DEBUG_SATS_MAN) {
      printf("</INIT_SATS_MANAGEMENT>\n");
    }
}
void new_init_sats_management(iterator_t *sdiffs, ptd_set_t *sats_ptd, 
                              ptd_set_t *sdiffs_ptd)
{
  if(length(sdiffs) == 0) {
    mk_empty(sats_ptd);
    mk_empty(sdiffs_ptd);
  } else {
    u8 ref_prn = new_choose_reference_sat(sdiffs);
    new_set_reference_sat_and_prns(ref_prn, sdiffs, sats_ptd, sdiffs_ptd);
  }
}
void box_init_sats_management(sats_management_t *sats_management,
                          u8 num_sdiffs, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  set_t sdiffs_set;
  set_state_t sdiffs_state;
  iterator_t sdiffs_itr;
  mk_set(&sdiffs_set, num_sdiffs, sizeof(sdiff_t), sdiffs, &sdiff_key);
  mk_set_itr(&sdiffs_itr, &sdiffs_state, &sdiffs_set);

  u16 num_sats = sats_management->num_sats;
  prn prn_buff[num_sats];
  ptd_set_t sats_ptd;
  sats_ptd.set.arr = prn_buff;

  ptd_set_t sdiffs_ptd;
  sdiff_t sdiff_buff[num_sdiffs];
  sdiffs_ptd.set.arr = sdiff_buff;

  new_init_sats_management(&sdiffs_itr, &sats_ptd, &sdiffs_ptd);

  sats_management->num_sats = sats_ptd.set.len;
  ptd_to_ref_fst(&sats_ptd, sats_management->prns);
  ptd_to_ref_fst(&sdiffs_ptd, sdiffs_with_ref_first);

}
void init_sats_management(sats_management_t *sats_management,
                          u8 num_sdiffs, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  old_init_sats_management(sats_management, num_sdiffs, sdiffs, sdiffs_with_ref_first);
  box_init_sats_management(sats_management, num_sdiffs, sdiffs, sdiffs_with_ref_first);

  sats_management_t man1, man2;
  sdiff_t ref_first1[MAX_CHANNELS];
  sdiff_t ref_first2[MAX_CHANNELS];
  old_init_sats_management(&man1, num_sdiffs, sdiffs, ref_first1);
  box_init_sats_management(&man2, num_sdiffs, sdiffs, ref_first2);
  assert(man1.num_sats == man2.num_sats);
  assert(man1.num_sats == num_sdiffs);
  assert(0 == memcmp(man1.prns, man2.prns, man1.num_sats * sizeof(prn)));
  assert(0 == memcmp(ref_first1, ref_first2, num_sdiffs * sizeof(sdiff_t)));

  memcpy(sats_management, &man1, sizeof(sats_management_t));
  memcpy(sdiffs_with_ref_first, ref_first1, num_sdiffs * sizeof(sdiff_t));
}
/* END init_sats_management */

// TODO iterator printing
// TODO prn print macro?
void print_sats_management(sats_management_t *sats_management)
{
  printf("sats_management->num_sats=%u\n", sats_management->num_sats);
  for (u8 i=0; i<sats_management->num_sats; i++) {
    printf("sats_management->prns[%u]= %u\n", i, sats_management->prns[i]);
  }
}

/** Updates sats to the new measurements' sat set
 */
// TODO intersect, contains?, call others
s8 rebase_sats_management(sats_management_t *sats_management,
                          u8 num_sdiffs, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  s8 return_code;
  u8 ref_prn;
  if (DEBUG_SATS_MAN) {
        printf("<REBASE_SATS_MANAGEMENT>\n");
  }

  if (sats_management->num_sats <= 1) {
    /* Need to init first. */
    // TODO test this change
    /* sdiffs_with_ref_first is overwritten later in this function. */
    init_sats_management(sats_management, num_sdiffs, sdiffs, sdiffs_with_ref_first);
  }
  

  // Check if old reference is in sdiffs
  if (bsearch(&(sats_management->prns[0]), sdiffs, num_sdiffs, sizeof(sdiff_t), &sdiff_search_prn)) {
    ref_prn = sats_management->prns[0];
    return_code = OLD_REF;
  }
  else {
    sdiff_t intersection_sats[num_sdiffs];
    u8 num_intersection = intersect_sats(sats_management->num_sats, num_sdiffs,
                                         &(sats_management->prns[1]), sdiffs, intersection_sats);
    if (num_intersection < INTERSECTION_SATS_THRESHOLD_SIZE) {
      if (DEBUG_SATS_MAN) {
        printf("</REBASE_SATS_MANAGEMENT>\n");
      }
      return NEW_REF_START_OVER;
    }
    else {
      if (DEBUG_SATS_MAN) {
        printf("sdiff prns= {");
        for (u8 yo_mama=0; yo_mama< num_sdiffs; yo_mama++) {
          printf("%u, ", sdiffs[yo_mama].prn);
        }
        printf("}\n");
        printf("sats_man_prns= {");
        for (u8 so_fetch=0; so_fetch < sats_management->num_sats; so_fetch++) {
          printf("%u, ", sats_management->prns[so_fetch]);
        }
        printf("}\n");
        
        printf("num intersect_sats= %u\nintersection= {", num_intersection);
        for (u8 bork=0; bork<num_intersection; bork++) {
          printf("%u, ", intersection_sats[bork].prn);
        }
        printf("}\n");
      }
      ref_prn = choose_reference_sat(num_intersection, intersection_sats);
      return_code = NEW_REF;
    }
  }
  set_reference_sat(ref_prn, sats_management,
                    num_sdiffs, sdiffs, sdiffs_with_ref_first);
  if (DEBUG_SATS_MAN) {
    printf("</REBASE_SATS_MANAGEMENT>\n");
  }
  return return_code;
}

// TODO map, freeze
void update_sats_sats_management(sats_management_t *sats_management, u8 num_non_ref_sdiffs, sdiff_t *non_ref_sdiffs)
{
  sats_management->num_sats = num_non_ref_sdiffs + 1;
  for (u8 i=1; i<num_non_ref_sdiffs+1; i++) {
    sats_management->prns[i] = non_ref_sdiffs[i-1].prn;
  }
}

// TODO intersect, subset? return code is not checked
s8 match_sdiffs_to_sats_man(sats_management_t *sats, u8 num_sdiffs, 
                            sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  u8 j = 1;
  u8 ref_prn = sats->prns[0];
  for (u8 i=0; i<num_sdiffs && j<sats->num_sats; i++) {
    if (sdiffs[i].prn == ref_prn) {
      memcpy(sdiffs_with_ref_first, &sdiffs[i], sizeof(sdiff_t));
    }
    else if (sdiffs[i].prn == sats->prns[j]) {
      memcpy(&sdiffs_with_ref_first[j], &sdiffs[i], sizeof(sdiff_t));
      j++;
    }
    else if (sdiffs[i].prn > sats->prns[j]) {
      return -1;
    }
  }
  return 0;
}
