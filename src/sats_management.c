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

#define DEBUG_SATS_MAN 0

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
u8 choose_reference_sat2(u8 num_sats, sdiff_t *sats)
{
  set_t set;
  mk_sdiff_set(&set, num_sats, sats);
  iterator_t it;
  set_state_t s;
  mk_set_itr(&it, &s, &set);

  prn_snr best;
  best.snr = -999;
  best.prn = -1;
  fold(&max_snr, NULL, &best, &it);

  return best.prn;
}
u8 choose_reference_sat(u8 num_sats, sdiff_t *sats)
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

// TODO intersection
//assumes both sets are ordered
u8 intersect_sats2(u8 num_sats, u8 num_sdiffs, u8 *sats, sdiff_t *sdiffs,
                   sdiff_t *intersection_sats)
{
  set_t sats_set, sdiffs_set;
  mk_prn_set(&sats_set, num_sats, sats);
  mk_sdiff_set(&sdiffs_set, num_sdiffs, sdiffs);
  iterator_t sats_itr, sdiffs_itr;
  set_state_t sats_state, sdiffs_state;
  mk_set_itr(&sats_itr,   &sats_state,   &sats_set);
  mk_set_itr(&sdiffs_itr, &sdiffs_state, &sdiffs_set);

  iterator_t intersection;
  intersection_state_t is;
  sdiff_t current;
  size_t sdiff_t_size = sizeof(sdiff_t);
  mk_intersection_itr(&intersection, &is, &current,
                      &sats_itr, &sdiffs_itr,
                      &prn_key, &sdiff_key,
                      &snd,
                      &sdiff_t_size);

  // TODO map snd
  int count = freeze_itr(sizeof(sdiff_t), MIN(num_sats, num_sdiffs),
                         intersection_sats, &intersection);
  assert(count != -1);
  return count;
}
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


// TODO is above function needed?
void set_reference_sat_and_prns(u8 ref_prn, sats_management_t *sats_management,
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

// TODO get rid of 0 check?
// no change
void init_sats_management(sats_management_t *sats_management,
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
    // Need to init first.
    init_sats_management(sats_management, num_sdiffs, sdiffs, 0);
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
