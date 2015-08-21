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

/** Chooses a reference sat from a set of sdiffs.
 * Currently just chooses the one with highest SNR.
 * \param num_sats  The length of the sdiffs array.
 * \param sats      An sdiff list from which the reference sat will be chosen.
 * \return The PRN of the chosen reference sat.
 */
static u8 choose_reference_sat(const u8 num_sats, const sdiff_t *sats)
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

/** Finds the intersection between a list of PRNs and a list of sdiffs.
 * Both are assumed to be ordered (ascending by PRN).
 * \param num_sats1         The number of satellites in the sats1 list.
 * \param num_sdiffs        The number of satellites in the sdiffs list.
 * \param sats1             An ascending list of PRNs.
 * \param sdiffs            A list of sdiffs ascending by PRN.
 * \param intersection_sats The sdiffs that are also in sats1.
 * \return The length of the intersection_sats.
 */
static u8 intersect_sats(const u8 num_sats1, const u8 num_sdiffs, const u8 *sats1,
                         const sdiff_t *sdiffs, sdiff_t *intersection_sats)
{
  DEBUG_ENTRY();

  u8 i, j, n = 0;
  if(DEBUG) {
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
      log_debug("(%u, %u, prn1=%u,\t %u, prn2=%u)\t\t prn1 < prn2  i++\n", n, i, sats1[i], j, sdiffs[j].prn);
      j--;
    }
    else if (sats1[i] > sdiffs[j].prn) {
      log_debug("(%u, %u, prn1=%u,\t %u, prn2=%u)\t\t prn1 > prn2  j++\n", n, i, sats1[i], j, sdiffs[j].prn);
      i--;
    }
    else {
      log_debug("(%u, %u, prn1=%u,\t %u, prn2=%u)\t\t prn1 = prn2  i,j,n++\n", n, i, sats1[i], j, sdiffs[j].prn);
      memcpy(&intersection_sats[n], &sdiffs[j], sizeof(sdiff_t));
      n++;
    }
  }
  if (DEBUG) {
    printf("intersection_sats= {");
    for (u8 k=0; k< n; k++) {
      printf("%u, ", intersection_sats[k].prn);
    }
    printf("}\n");
  }

  DEBUG_EXIT();
  return n;
}

/* TODO use the set abstraction fnoble is working on. */
/** Changes the ref PRN of a PRN array.
 * Assumes that the tail of the array is sorted (and keeps it sorted).
 * Assumes the new ref_prn is also in the original array.
 * \param ref_prn   The PRN to set as reference.
 * \param num_sats  The length of the PRN array.
 * \param prns      On input, an array with ref_prn in it, whose elements are
 *                  unique and whose tail is sorted.
 *                  On exit, an array with the same elements, whose tail is
 *                  sorted, and whose head is ref_prn.
 */
void set_reference_sat_of_prns(const u8 ref_prn, const u8 num_sats, u8 *prns)
{
  u8 old_ref = prns[0];
  u8 j;
  /* Return early if we're setting the ref to what it already is. */
  if (old_ref == ref_prn) {
    return;
  }
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
        set_old_yet++;
      }
      else {
        prns[j] = old_prns[i];
        j++;
      }
    }
  }
  if (set_old_yet == 0) {
    prns[j] = old_ref;
    set_old_yet++;
  }
  assert(set_old_yet == 1);
}


/* TODO use the set abstraction fnoble is working on. */
/** Changes the reference PRN of a sdiff set (in a copy) and a sats_management_t
 * Assumes that sats_management.prns has unique elements and a sorted tail.
 * Assumes that sats_management.prns has the same elements as the PRNs of
 * sdiffs (and thus that sats_management.num_sats == num_sdiffs).
 * Assumes that sdiffs is sorted ascending by PRN.
 * Assumes that ref_prn is in the both input lists.
 * \params ref_prn                The PRN to be used as reference in
 *                                sats_management and sdiffs_with_ref_first.
 * \params sats_management        The sats_management_t to be changed.
 *                                On exit, its .prns field will have the same
 *                                elements as on entry, but will have ref_prn
 *                                as the head and a sorted tail.
 * \params num_sdiffs             The length of the sdiffs list.
 * \params sdiffs                 A list of sdiffs to be copied to
 *                                sdiffs_with_ref_first with ref_prn's sdiff
 *                                first.
 * \params sdiffs_with_ref_first  A list of sdiffs with the same elemements as
 *                                sdiffs, but with the sdiff whose PRN is
 *                                ref_prn first, and with the tail to be sorted
 *                                by PRN.
 */
static void set_reference_sat(const u8 ref_prn, sats_management_t *sats_management,
                              const u8 num_sdiffs, const sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  DEBUG_ENTRY();

  u8 old_ref = sats_management->prns[0];
  log_debug("ref_prn = %u\n", ref_prn);
  log_debug("old_ref = %u\n", old_ref);
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
          set_old_yet++;
        }
        else {
          sats_management->prns[j] = old_prns[i];
          j++;
        }
      }
    }
    if (set_old_yet == 0) {
      sats_management->prns[j] = old_ref;
      set_old_yet++;
    }
    assert(set_old_yet == 1);
  }

  j=1;
  for (u8 i=0; i<num_sdiffs; i++) {
    if (sdiffs[i].prn != ref_prn) {
      log_debug("prn[%u] = %u\n", j, sdiffs[i].prn);
      memcpy(&sdiffs_with_ref_first[j], &sdiffs[i], sizeof(sdiff_t));
      j++;
    } else {
      log_debug("prn[0] = %u\n", sdiffs[i].prn);
      memcpy(&sdiffs_with_ref_first[0], &sdiffs[i], sizeof(sdiff_t));
    }
  }

  DEBUG_EXIT();
}

/** Copies sdiffs with ref_prn as reference sat, filling in sats_management.
 * Assumes sdiffs is sorted by PRN (ascending).
 * Assumes sdiffs has an sdiff whose PRN is ref_prn.
 * \param ref_prn               The reference sat to set sats_management and
 *                              sdiffs_with_ref_first with.
 * \param sats_management       On exit, .num_sats = num_sdiffs, and .prns
 *                              is all of the prns of sdiffs, with ref_prn at
 *                              the head and a sorted tail (ascending).
 * \param num_sdiffs            The length of the sdiffs array.
 * \param sdiffs                The sdiffs used to populate sats_management and
 *                              sdiffs_with_ref_first. Sorted by PRN.
 * \param sdiffs_with_ref_first On exit, has the same elements as sdiffs, but
 *                              with the sdiff corresponding to ref_prn at the
 *                              head and the tail sorted.
 */
static void set_reference_sat_and_prns(const u8 ref_prn, sats_management_t *sats_management,
                                       const u8 num_sdiffs, const sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
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

/** Chooses a reference sat and sets up sats_management, while costructing
 * sdiffs_with_ref_first.
 * Assumes sdiffs is sorted by PRN.
 * \param sats_management       On exit, .num_sats is num_sdiffs, and .prns
 *                              has all the same PRNs as sdiffs, while its head
 *                              is a chosen reference (same ref as
 *                              sdiffs_with_ref_first) and its tail is sorted.
 *                              The order of .prns should match the order of
 *                              sdiffs_with_ref_first.
 * \param num_sdiffs            The length of the sdiffs array.
 * \param sdiffs                A set of sdiffs sorted by PRN.
 * \param sdiffs_with_ref_first On exit is a permutation of sdiffs. Shares
 *                              order with sats_management.prns. See that arg.
 */
void init_sats_management(sats_management_t *sats_management,
                          const u8 num_sdiffs, const sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  DEBUG_ENTRY();

  if (num_sdiffs == 0) {
    sats_management->num_sats = 0;
    DEBUG_EXIT();
    return;
  }

  u8 ref_prn = choose_reference_sat(num_sdiffs, sdiffs);
  set_reference_sat_and_prns(ref_prn, sats_management,
                             num_sdiffs, sdiffs, sdiffs_with_ref_first);
  DEBUG_EXIT();
}

/** Prints one prn per line */
void print_sats_management(sats_management_t *sats_management)
{
  printf("sats_management->num_sats=%u\n", sats_management->num_sats);
  for (u8 i=0; i<sats_management->num_sats; i++) {
    printf("sats_management->prns[%u]= %u\n", i, sats_management->prns[i]);
  }
}
/** Prints all prns on one line */
void print_sats_management_short(sats_management_t *sats_man) {
  printf("sats_management sats: ");
  for (u8 i=0; i<sats_man->num_sats; i++) {
    printf("%d,", sats_man->prns[i]);
  }
  printf("\n");
}


/** Rebases sats_management to have a reference sat in sdiffs, if possible.
 * Otherwise returns a reset code.
 * Initializes sats_management if empty.
 * Assumes that sdiffs is sorted by PRN (ascending).
 * \param sats_management       If, on entry, is empty (<= 1 sat), then on
 *                              exit it will have all the sats of sdiff, with
 *                              a reference chosen.
 *                              If, on entry, it's not empty, then on exit it
 *                              will be unchanged if its reference is in
 *                              sdiffs, it will be updated to a new reference
 *                              if it shares enough sats with sdiffs to do so,
 *                              otherwise unchanged.
 * \param num_sdiffs            The length of sdiffs.
 * \param sdiffs                A list of sdiffs sorted by PRN.
 * \param sdiffs_with_ref_first On exit, if sats_management on exit had a
 *                              reference in sdiffs, then this will have the
 *                              same elements as sdiffs, with that reference at
 *                              the head, and a sorted tail. Otherwise it's
 *                              unchanged.
 * \return  If, on exit, sats_management doesn't have a reference sat in the
 *          sdiffs set, then NEW_REF_START_OVER. Otherwise, if the ref is the
 *          same as the old ref, then it's OLD_REF, else it's NEW_REF.
 */
s8 rebase_sats_management(sats_management_t *sats_management,
                          const u8 num_sdiffs, const sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  DEBUG_ENTRY();

  s8 return_code;
  u8 ref_prn;

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
    u8 num_intersection = intersect_sats(sats_management->num_sats-1, num_sdiffs,
                                         &(sats_management->prns[1]), sdiffs, intersection_sats);
    if (num_intersection < INTERSECTION_SATS_THRESHOLD_SIZE) {
      DEBUG_EXIT();
      return NEW_REF_START_OVER;
    }
    else {
      if (DEBUG) {
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

  DEBUG_EXIT();
  return return_code;
}

/** Sets sats_management's non reference PRNs (and num_sats).
 * \param sats_management     On exit, has .num_sats = num_non_ref_sdiffs +1,
 *                            and the tail of .prns = non_ref_sdiffs[*].prn
 * \param num_non_ref_sdiffs  The length of non_ref_sdiffs.
 * \param non_ref_sdiffs      The sdiffs to set the tail of sats_management with.
 */
void update_sats_sats_management(sats_management_t *sats_management,
                                 u8 num_non_ref_sdiffs,
                                 const sdiff_t *non_ref_sdiffs)
{
  sats_management->num_sats = num_non_ref_sdiffs + 1;
  for (u8 i=1; i<num_non_ref_sdiffs+1; i++) {
    sats_management->prns[i] = non_ref_sdiffs[i-1].prn;
  }
}


/** Creates an sdiff set from sdiffs whose order matches a sats_management_t.
 * \param sats                  A PRN set to give the desired sat set and ref.
 * \param num_sdiffs            The length of sdiffs.
 * \param sdiffs                A set of sdiffs from which to take the PRNs
 *                              in sats.
 * \param sdiffs_with_ref_first On exit, if sats.prns is a subset of sdiffs,
 *                              then it will be the sdiffs corresponding to the
 *                              elements of sats.prns in sats.prns's order.
 * \return If sats.prns is a subset of sdiffs, then 0, else -1.
 */
s8 match_sdiffs_to_sats_man(const sats_management_t *sats, u8 num_sdiffs,
                            const sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
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

