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

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "logging.h"
#include "observation.h"
#include "sats_management.h"
#include "linear_algebra.h"

gnss_signal_t choose_reference_sat(const u8 num_sats, const sdiff_t *sats)
{
  double best_snr=sats[0].snr;
  gnss_signal_t best_sid = sats[0].sid;
  for (u8 i=1; i<num_sats; i++) {
    if (sats[i].snr > best_snr) {
      best_snr = sats[i].snr;
      best_sid = sats[i].sid;
    }
  }
  return best_sid;
}

//assumes both sets are ordered
static u8 intersect_sats(const u8 num_sats1, const u8 num_sdiffs, const gnss_signal_t *sats1,
                         const sdiff_t *sdiffs, sdiff_t *intersection_sats)
{
  DEBUG_ENTRY();

  u8 i, j, n = 0;
  if(DEBUG) {
    printf("sdiff prns= {");
    for (u8 k=0; k< num_sdiffs; k++) {
      printf("%u, ", sdiffs[k].sid.sat);
    }
    printf("}\n");
    printf("sats= {");
    for (u8 k=0; k< num_sats1; k++) {
      printf("%u, ", sats1[k].sat);
    }
    printf("}\n");
    printf("(n, i, prn1[i],\t j, prn2[j])\n");
  }
  /* Loop over sats1 and sdffs and check if a PRN is present in both. */
  for (i=0, j=0; i<num_sats1 && j<num_sdiffs; i++, j++) {
    if (sid_compare(sats1[i], sdiffs[j].sid) < 0) {
      log_debug("(%u, %u, prn1=%u,\t %u, prn2=%u)\t\t prn1 < prn2  i++", n, i, sats1[i], j, sdiffs[j].sid.sat);
      j--;
    }
    else if (sid_compare(sats1[i], sdiffs[j].sid) > 0) {
      log_debug("(%u, %u, prn1=%u,\t %u, prn2=%u)\t\t prn1 > prn2  j++", n, i, sats1[i], j, sdiffs[j].sid.sat);
      i--;
    }
    else {
      log_debug("(%u, %u, prn1=%u,\t %u, prn2=%u)\t\t prn1 = prn2  i,j,n++", n, i, sats1[i], j, sdiffs[j].sid.sat);
      memcpy(&intersection_sats[n], &sdiffs[j], sizeof(sdiff_t));
      n++;
    }
  }
  if (DEBUG) {
    printf("intersection_sats= {");
    for (u8 k=0; k< n; k++) {
      printf("%u, ", intersection_sats[k].sid.sat);
    }
    printf("}\n");
  }

  DEBUG_EXIT();
  return n;
}

/** Updates the ref prn for a sorted prn array.
 *  Inserts the old ref prn so the tail of the array is sorted.
 */
/* TODO use the set abstraction fnoble is working on. */
void set_reference_sat_of_sids(const gnss_signal_t ref_sid, const u8 num_sats, gnss_signal_t *sids)
{
  gnss_signal_t old_ref = sids[0];
  u8 j;
  if (!sid_is_equal(old_ref, ref_sid)) {
    j = 1;
    gnss_signal_t old_sids[num_sats];
    memcpy(old_sids, sids, num_sats * sizeof(gnss_signal_t));
    u8 set_old_yet = 0;
    sids[0] = ref_sid;
    for (u8 i=1; i<num_sats; i++) {
      if (!sid_is_equal(old_sids[i], ref_sid)) {
        if ((sid_compare(old_sids[i], old_ref) > 0) &&
            set_old_yet == 0) {
          sids[j] = old_ref;
          j++;
          i--;
          set_old_yet++;
        }
        else {
          sids[j] = old_sids[i];
          j++;
        }
      }
    }
    if (set_old_yet == 0) {
      sids[j] = old_ref;
      set_old_yet++;
    }
    assert(set_old_yet == 1);
  }
}


/** Puts sdiffs into sdiffs_with_ref_first with the sdiff for ref_prn first, while updating sats_management
 *  Inserts the old ref prn so the tail of the sats_management array is sorted.
 */
/* TODO use the set abstraction fnoble is working on. */
static void set_reference_sat(const gnss_signal_t ref_sid, sats_management_t *sats_management,
                              const u8 num_sdiffs, const sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  DEBUG_ENTRY();

  gnss_signal_t old_ref = sats_management->sids[0];
  log_debug("ref_sid = %u", ref_sid);
  log_debug("old_ref = %u", old_ref);
  u8 j;
  if (!sid_is_equal(old_ref, ref_sid)) {
    j = 1;
    gnss_signal_t old_sids[sats_management->num_sats];
    memcpy(old_sids, sats_management->sids, sats_management->num_sats * sizeof(gnss_signal_t));
    u8 set_old_yet = 0;
    sats_management->sids[0] = ref_sid;
    for (u8 i=1; i<sats_management->num_sats; i++) {
      if (!sid_is_equal(old_sids[i], ref_sid)) {
        if ((sid_compare(old_sids[i], old_ref) > 0) && set_old_yet == 0) {
          sats_management->sids[j] = old_ref;
          j++;
          i--;
          set_old_yet++;
        }
        else {
          sats_management->sids[j] = old_sids[i];
          j++;
        }
      }
    }
    if (set_old_yet == 0) {
      sats_management->sids[j] = old_ref;
      set_old_yet++;
    }
    assert(set_old_yet == 1);
  }

  j=1;
  for (u8 i=0; i<num_sdiffs; i++) {
    if (!sid_is_equal(sdiffs[i].sid, ref_sid)) {
      log_debug("prn[%u] = %u", j, sdiffs[i].sid.sat);
      memcpy(&sdiffs_with_ref_first[j], &sdiffs[i], sizeof(sdiff_t));
      j++;
    } else {
      log_debug("prn[0] = %u", sdiffs[i].sid.sat);
      memcpy(&sdiffs_with_ref_first[0], &sdiffs[i], sizeof(sdiff_t));
    }
  }

  DEBUG_EXIT();
}

static void set_reference_sat_and_sids(const gnss_signal_t ref_sid, sats_management_t *sats_management,
                                       const u8 num_sdiffs, const sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  sats_management->num_sats = num_sdiffs;
  sats_management->sids[0] = ref_sid;
  u8 j=1;
  for (u8 i=0; i<num_sdiffs; i++) {
    if (!sid_is_equal(sdiffs[i].sid, ref_sid)) {
      sats_management->sids[j] = sdiffs[i].sid;
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

void init_sats_management(sats_management_t *sats_management,
                          const u8 num_sdiffs, const sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  DEBUG_ENTRY();

  if (num_sdiffs == 0) {
    sats_management->num_sats = 0;
    DEBUG_EXIT();
    return;
  }

  gnss_signal_t ref_sid = choose_reference_sat(num_sdiffs, sdiffs);
  set_reference_sat_and_sids(ref_sid, sats_management,
                             num_sdiffs, sdiffs, sdiffs_with_ref_first);
  DEBUG_EXIT();
}

/** Prints one prn per line */
void print_sats_management(sats_management_t *sats_management)
{
  printf("sats_management->num_sats=%u\n", sats_management->num_sats);
  for (u8 i=0; i<sats_management->num_sats; i++) {
    printf("sats_management->prns[%u]= %u\n", i, sats_management->sids[i].sat);
  }
}
/** Prints all prns on one line */
void print_sats_management_short(sats_management_t *sats_man) {
  printf("sats_management sats: ");
  for (u8 i=0; i<sats_man->num_sats; i++) {
    printf("%d,", sats_man->sids[i].sat);
  }
  printf("\n");
}


/** Updates sats to the new measurements' sat set
 */
s8 rebase_sats_management(sats_management_t *sats_management,
                          const u8 num_sdiffs, const sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  DEBUG_ENTRY();

  s8 return_code;
  gnss_signal_t ref_sid;

  if (sats_management->num_sats <= 1) {
    // Need to init first.
    init_sats_management(sats_management, num_sdiffs, sdiffs, 0);
  }

  // Check if old reference is in sdiffs
  if (bsearch(&(sats_management->sids[0]), sdiffs, num_sdiffs, sizeof(sdiff_t), cmp_sid_sdiff)) {
    ref_sid = sats_management->sids[0];
    return_code = OLD_REF;
  }
  else {
    sdiff_t intersection_sats[num_sdiffs];
    u8 num_intersection = intersect_sats(sats_management->num_sats-1, num_sdiffs,
                                         &(sats_management->sids[1]), sdiffs, intersection_sats);
    if (num_intersection < INTERSECTION_SATS_THRESHOLD_SIZE) {
      DEBUG_EXIT();
      return NEW_REF_START_OVER;
    }
    else {
      if (DEBUG) {
        printf("sdiff prns= {");
        for (u8 yo_mama=0; yo_mama< num_sdiffs; yo_mama++) {
          printf("%u, ", sdiffs[yo_mama].sid.sat);
        }
        printf("}\n");
        printf("sats_man_prns= {");
        for (u8 so_fetch=0; so_fetch < sats_management->num_sats; so_fetch++) {
          printf("%u, ", sats_management->sids[so_fetch].sat);
        }
        printf("}\n");
        printf("num intersect_sats= %u\nintersection= {", num_intersection);
        for (u8 bork=0; bork<num_intersection; bork++) {
          printf("%u, ", intersection_sats[bork].sid.sat);
        }
        printf("}\n");
      }
      ref_sid = choose_reference_sat(num_intersection, intersection_sats);
      return_code = NEW_REF;
    }
  }
  set_reference_sat(ref_sid, sats_management,
                    num_sdiffs, sdiffs, sdiffs_with_ref_first);

  DEBUG_EXIT();
  return return_code;
}

void update_sats_sats_management(sats_management_t *sats_management, u8 num_non_ref_sdiffs, sdiff_t *non_ref_sdiffs)
{
  sats_management->num_sats = num_non_ref_sdiffs + 1;
  for (u8 i=1; i<num_non_ref_sdiffs+1; i++) {
    sats_management->sids[i] = non_ref_sdiffs[i-1].sid;
  }
}



s8 match_sdiffs_to_sats_man(sats_management_t *sats, u8 num_sdiffs,
                            sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  u8 j = 1;
  gnss_signal_t ref_sid = sats->sids[0];
  for (u8 i=0; i<num_sdiffs && j<sats->num_sats; i++) {
    if (sid_is_equal(sdiffs[i].sid, ref_sid)) {
      memcpy(sdiffs_with_ref_first, &sdiffs[i], sizeof(sdiff_t));
    }
    else if (sid_is_equal(sdiffs[i].sid, sats->sids[j])) {
      memcpy(&sdiffs_with_ref_first[j], &sdiffs[i], sizeof(sdiff_t));
      j++;
    }
    else if (sid_compare(sdiffs[i].sid, sats->sids[j]) > 0) {
      return -1;
    }
  }
  return 0;
}

