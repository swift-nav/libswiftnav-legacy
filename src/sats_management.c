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
#include "single_diff.h"
#include "sats_management.h"
#include "linear_algebra.h"


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

// bool contains(u8 num_sats, u8 ref_prn, sdiff_t *sdiffs) 
// {
//   //todo make and use assumptions about the ordering of the prns
//   for (u8 i=0; i<num_sats; i++) {
//     if (ref_prn == sdiffs[i].prn) {
//       return true;
//     }
//   }
//   return false;
// }

//assumes both sets are ordered
u8 intersect_sats(u8 num_sats1, u8 num_sdiffs, u8 *sats1, sdiff_t *sdiffs,
                  sdiff_t *intersection_sats)
{
  u8 i, j, n = 0;

  /* Loop over sats1 and sdffs and check if a PRN is present in both. */
  for (i=0, j=0; i<num_sats1 && j<num_sdiffs; i++, j++) {
    if (sats1[i] < sdiffs[j].prn)
      j--;
    else if (sats1[i] > sdiffs[j].prn)
      i--;
    else {
      memcpy(&intersection_sats[n], &sdiffs[i], sizeof(sdiff_t));
      n++;
    }
  }
  return n;
}

// s8 check_and_rebase_reference_sat(kf_t *kf, u8 num_sats, sdiff_t *sdiffs, kf_state_t *state)
// {
//   bool still_has_ref = contains(num_sats, sdiffs);
//   if (still_has_ref == false) {
//     sdiff_t common_sats[num_sats];
//     u8 intersection_size = intersect_sats(num_sats, sdiffs, common_sats);
//     if (intersection_size>=4) {
//       u8 new_ref_prn = choose_reference_sat(intersection_size, common_sats);
//       printf("%i\n", new_ref_prn);
//       printf("%i\n", state->state_dim);
//     }
//   }
//   return 1;
// }

// u8 update_filter(kf_t *kf, u8 num_sats, sdiff_t *sdiffs, kf_state_t *state)
// {
//   //TODO add update_E logic

//   u8 control_logic = check_and_rebase_reference_sat(kf, num_sats, sdiffs, state);
// }

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


/** Puts sdiffs into sdiffs_with_ref_first with the sdiff for ref_prn first
 */
void set_reference_sat(u8 ref_prn, sats_management_t *sats_management,
                          u8 num_sdiffs, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  u8 old_ref = sats_management->prns[0];
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
      memcpy(&sdiffs_with_ref_first[j], &sdiffs[i], sizeof(sdiff_t));
      j++;
    } else {
      memcpy(&sdiffs_with_ref_first[0], &sdiffs[i], sizeof(sdiff_t));
    }
  }
}

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

void init_sats_management(sats_management_t *sats_management,
                          u8 num_sdiffs, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  if (num_sdiffs == 0) {
    sats_management->num_sats = 0;
    return;
  }

  u8 ref_prn = choose_reference_sat(num_sdiffs, sdiffs);
  set_reference_sat_and_prns(ref_prn, sats_management,
                             num_sdiffs, sdiffs, sdiffs_with_ref_first);
}

void print_sats_management(sats_management_t *sats_management)
{
  printf("sats_management->num_sats=%u\n", sats_management->num_sats);
  for (u8 i=0; i<sats_management->num_sats; i++) {
    printf("sats_management->prns[%u]= %u\n", i, sats_management->prns[i]);
  }
}

/** Updates sats to the new measurements' sat set
 */
s8 rebase_sats_management(sats_management_t *sats_management,
                          u8 num_sdiffs, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  s8 return_code;
  u8 ref_prn;

  if (sats_management->num_sats == 0) {
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
      return NEW_REF_START_OVER;
    }
    else {
      ref_prn = choose_reference_sat(num_intersection, intersection_sats);
      return_code = NEW_REF;
    }
  }
  set_reference_sat(ref_prn, sats_management,
                    num_sdiffs, sdiffs, sdiffs_with_ref_first);
  return return_code;
}

void update_sats_sats_management(sats_management_t *sats_management, u8 num_non_ref_sdiffs, sdiff_t *non_ref_sdiffs)
{
  sats_management->num_sats = num_non_ref_sdiffs + 1;
  for (u8 i=1; i<num_non_ref_sdiffs+1; i++) {
    sats_management->prns[i] = non_ref_sdiffs[i-1].prn;
  }
}







