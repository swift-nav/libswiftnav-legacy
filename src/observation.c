/*
 * Copyright (C) 2014 Swift Navigation Inc.
 * Contact: Fergus Noble <fergus@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#include <string.h>
#include <assert.h>

#include "logging.h"
#include "linear_algebra.h"
#include "ephemeris.h"
#include "constants.h"
#include "sats_management.h"
#include "set.h"
#include "observation.h"

/** \defgroup single_diff Single Difference Observations
 * Functions for storing and manipulating single difference observations.
 * \{ */

/** Comparison function for `sdiff_t` by PRN.
 * See `cmp_fn`. */
int cmp_sdiff_prn(const void *a_, const void *b_)
{
  const sdiff_t *a = (const sdiff_t *)a_;
  const sdiff_t *b = (const sdiff_t *)b_;

  u16 a_prn = a->sid.prn;
  u16 b_prn = b->sid.prn;

  return cmp_u8_u8(&a_prn, &b_prn);
}

/** Create a single difference from two observations.
 * Used by single_diff() to map two `navigation_measurement_t`s
 * into an `sdiff_t`.
 *
 * SNR in the output is the lesser of the SNRs of inputs a and b.
 *
 * `sat_pos` and `sat_vel` are taken from input b.
 *
 * Called once for each pair of observations with matching PRNs.
 */
static void single_diff_(void *context, u32 n, const void *a, const void *b)
{
  const navigation_measurement_t *m_a = (const navigation_measurement_t *)a;
  const navigation_measurement_t *m_b = (const navigation_measurement_t *)b;
  sdiff_t *sds = (sdiff_t *)context;

  sds[n].sid.prn = m_a->sid.prn;
  sds[n].pseudorange = m_a->raw_pseudorange - m_b->raw_pseudorange;
  sds[n].carrier_phase = m_a->carrier_phase - m_b->carrier_phase;
  sds[n].doppler = m_a->raw_doppler - m_b->raw_doppler;
  sds[n].snr = MIN(m_a->snr, m_b->snr);
  sds[n].lock_counter = m_a->lock_counter + m_b->lock_counter;

  /* NOTE: We use the position and velocity from B (this is required by
   * make_propagated_sdiffs(). */
  memcpy(&(sds[n].sat_pos), &(m_b->sat_pos), 3*sizeof(double));
  memcpy(&(sds[n].sat_vel), &(m_b->sat_vel), 3*sizeof(double));
}

/** Calculate single differences from two sets of observations.
 * Undifferenced input observations are assumed to be both taken at the
 * same time, `t`.
 *
 * SNR in the output is the lesser of the SNRs of inputs a and b.
 *
 * `sat_pos` and `sat_vel` are taken from input b.
 *
 * \param n_a Number of measurements in set `m_a`
 * \param m_a Array of undifferenced observations, as a set sorted by PRN
 * \param n_b Number of measurements in set `m_b`
 * \param m_b Array of undifferenced observations, as a set sorted by PRN
 * \param sds Single difference observations
 *
 * \return The number of observations written to `sds` on success,
 *         -1 if `m_a` is not a valid set,
 *         -2 if `m_b` is not a valid set
 */
u8 single_diff(u8 n_a, navigation_measurement_t *m_a,
               u8 n_b, navigation_measurement_t *m_b,
               sdiff_t *sds)
{
  return intersection_map(n_a, sizeof(navigation_measurement_t), m_a,
                          n_b, sizeof(navigation_measurement_t), m_b,
                          nav_meas_cmp, sds, single_diff_);
}

typedef struct {
  sdiff_t *sds;
  double *remote_pos_ecef;
} make_propagated_sdiff_ctxt;

static void make_propagated_sdiff_(void *context, u32 n,
                                   const void *a, const void *b)
{
  const navigation_measurement_t *m_a = (const navigation_measurement_t *)a;
  const navigation_measurement_t *m_b = (const navigation_measurement_t *)b;
  make_propagated_sdiff_ctxt *ctxt = (make_propagated_sdiff_ctxt *)context;

  /* Construct sds[n].
   * NOTE: The resulting sds have sat_pos and sat_vel taken from b. */
  single_diff_(ctxt->sds, n, a, b);

  double old_dist = vector_distance(3, m_a->sat_pos, ctxt->remote_pos_ecef);
  double new_dist = vector_distance(3, m_b->sat_pos, ctxt->remote_pos_ecef);
  double dr = new_dist - old_dist;

  ctxt->sds[n].pseudorange += dr;
  ctxt->sds[n].carrier_phase -= dr / GPS_L1_LAMBDA_NO_VAC;
}

/** New more efficient version of make_propagated_sdiffs(), NOT WORKING!!!
 * WIP */
u8 make_propagated_sdiffs_wip(u8 n_local, navigation_measurement_t *m_local,
                              u8 n_remote, navigation_measurement_t *m_remote,
                              double remote_pos_ecef[3], sdiff_t *sds)
{
  make_propagated_sdiff_ctxt ctxt = {
    .sds = sds,
    .remote_pos_ecef = remote_pos_ecef
  };
  return intersection_map(n_local, sizeof(navigation_measurement_t), m_local,
                          n_remote, sizeof(navigation_measurement_t), m_remote,
                          nav_meas_cmp, &ctxt, make_propagated_sdiff_);
}

int sdiff_search_prn(const void *a, const void *b)
{
  return (*(u8*)a - ((sdiff_t *)b)->sid.prn);
}

/** Propagates remote measurements to a local time and makes sdiffs.
 * When we get two sets of observations that aren't time matched to each
 * other (but are internally time matched within each set), we need to
 * adjust one set of measurements to be our best guess of what it would have
 * been had we measured it at the other set's time. This function does that
 * and differences those measurements from sats present in both sets.
 *
 * It returns the number of sats common in both.
 *
 * \remark This is actually using the sat positions at the time the receiver
 *         got the signal. You can backtrack to the sat position at the time
 *         the signal was sent. At the very least, if you backtrack assuming
 *         sat constant velocity, the difference is miniscule.
 *
 * \todo  Integrate this with single_diff via a higher order function.
 *
 * \param n_local           The number of measurements taken locally.
 * \param m_local           The measurements taken locally (sorted by prn).
 * \param n_remote          The number of measurements taken remotely.
 * \param m_remote          THe measurements taken remotely (sorted by prn).
 * \param remote_dists      The distances from the remote receiver to each
 *                           satellite at the time the remote measurements
 *                           were taken (i-th element of this list must
 *                           correspond to the i-th element of m_remote).
 * \param remote_pos_ecef   The position of the remote receiver (presumed
 *                           constant in ecef).
 * \param es
 * \param t
 * \param sds               The single differenced propagated measurements.
 * \return The number of sats common in both local and remote sdiffs.
 */
u8 make_propagated_sdiffs(u8 n_local, navigation_measurement_t *m_local,
                          u8 n_remote, navigation_measurement_t *m_remote,
                          double *remote_dists, double remote_pos_ecef[3],
                          ephemeris_kepler_t *es, gps_time_t t,
                          sdiff_t *sds)
{
  u8 i, j, n = 0;

  /* Loop over m_a and m_b and check if a PRN is present in both. */
  for (i=0, j=0; i<n_local && j<n_remote; i++, j++) {
    ephemeris_t eph;
    if (m_local[i].sid.constellation == GPS_CONSTELLATION) {
      eph.ephemeris_kep = &es[m_local[i].sid.prn];
      eph.ephemeris_xyz = NULL;
    } else {
      //TODO
      continue;
    }
    if (m_local[i].sid.prn < m_remote[j].sid.prn)
      j--;
    else if (m_local[i].sid.prn > m_remote[j].sid.prn)
      i--;
    else if (ephemeris_good(&eph, m_local[i].sid, t)) {
      double clock_err;
      double clock_rate_err;
      double local_sat_pos[3];
      double local_sat_vel[3];
      legacy_calc_sat_state(&es[m_local[i].sid.prn], t, local_sat_pos, local_sat_vel,
                     &clock_err, &clock_rate_err);
      sds[n].sid.prn = m_local[i].sid.prn;
      sds[n].sid.constellation = m_local[i].sid.constellation;
      sds[n].sid.band = m_local[i].sid.band;
      double dx = local_sat_pos[0] - remote_pos_ecef[0];
      double dy = local_sat_pos[1] - remote_pos_ecef[1];
      double dz = local_sat_pos[2] - remote_pos_ecef[2];
      double new_dist = sqrt( dx * dx + dy * dy + dz * dz);
      double dist_diff = new_dist - remote_dists[j];
      /* Explanation:
       * pseudorange = dist + c
       * To update a pseudorange in time:
       *  new_pseudorange = new_dist + c
       *                  = old_dist + c + (new_dist - old_dist)
       *                  = old_pseudorange + (new_dist - old_dist)
       *
       * So to get the single differenced pseudorange:
       *  local_pseudorange - new_remote_pseudorange
       *    = local_pseudorange - (old_remote_pseudorange + new_dist - old_dist)
       *
       * For carrier phase, it's the same thing, but the update has opposite sign. */
      sds[n].pseudorange = m_local[i].raw_pseudorange
                         - (m_remote[j].raw_pseudorange
                            + dist_diff);
      sds[n].carrier_phase = m_local[i].carrier_phase
                           - (m_remote[j].carrier_phase
                              - dist_diff / GPS_L1_LAMBDA);

      /* Doppler is not propagated.
       * sds[n].doppler = m_local[i].raw_doppler - m_remote[j].raw_doppler; */
      sds[n].snr = MIN(m_local[i].snr, m_remote[j].snr);
      memcpy(&(sds[n].sat_pos), &(local_sat_pos[0]), 3*sizeof(double));
      memcpy(&(sds[n].sat_vel), &(local_sat_vel[0]), 3*sizeof(double));

      n++;
    }
  }

  return n;
}

/* Checks to see if any satellites have had their lock counter values have
 * changed.
 *
 * If the lock counter changes it indicates that the satellite should be
 * reinitialized in the filter. This function checks the current lock counter
 * values in a set of single difference observations with a stored state
 * (`lock_counters`).
 *
 * It outputs a the number of satellites which have had their lock counter
 * change and a list of the corresponding PRNs.
 *
 * The lock counter values are then updated to reflect their changed values.
 *
 * \param n_sds        Number of single difference observations passed in
 * \param sds          Array of single difference observations
 * \param lock_counter Array of lock counter values, indexed by PRN
 * \param sats_to_drop Output array of PRNs for which the lock counter value
 *                     has changed
 * \return Number of sats with changed lock counter
 */
u8 check_lock_counters(u8 n_sds, const sdiff_t *sds, u16 *lock_counters,
                       u8 *sats_to_drop)
{
  assert(sds != NULL);
  assert(lock_counters != NULL);
  assert(sats_to_drop != NULL);

  u8 num_sats_to_drop = 0;
  for (u8 i = 0; i<n_sds; i++) {
    u16 prn = sds[i].sid.prn;
    u16 new_count = sds[i].lock_counter;
    if (new_count != lock_counters[prn]) {
      sats_to_drop[num_sats_to_drop++] = prn;
      lock_counters[prn] = new_count;
    }
  }
  return num_sats_to_drop;
}

/* Constructs the double differenced measurements and sdiffs needed for IAR.
 * This requires that all IAR prns are in the sdiffs used (sdiffs is a superset
 * of IAR's PRNs), that the sdiffs are ordered by prn, and that the IAR
 * non_ref_prns are also ordered (both ascending).
 *
 * The dd_meas output will be ordered like the non_ref_prns
 * with the whole vector of carrier phase elements coming before the
 * pseudoranges. The amb_sdiffs will have the ref sat first, then the rest in
 * ascending order like the non_ref_prns.
 *
 * \param ref_prn                    The current reference PRN for the IAR.
 * \param non_ref_prns               The rest of the current PRNs for the IAR.
 * \param num_dds                    The number of dds used in the IAR
 *                                   (length of non_ref_prns).
 * \param num_sdiffs                 The number of sdiffs being passed in.
 * \param sdiffs                     The sdiffs to pull measurements out of.
 * \param dd_meas  The output vector of DD measurements
 *                                   to be used to update the IAR.
 * \param sdiffs_out                 The sdiffs that correspond to the IAR PRNs.
 *                                   Should have length = num_dds+1.
 * \return 0 if the input sdiffs are superset of the IAR sats,
 *        -1 if they are not,
 *        -2 if non_ref_prns is not an ordered set.
 */
s8 make_dd_measurements_and_sdiffs(u8 ref_prn, const u8 *non_ref_prns, u8 num_dds,
                                   u8 num_sdiffs, const sdiff_t *sdiffs_in,
                                   double *dd_meas, sdiff_t *sdiffs_out)
{
  DEBUG_ENTRY();

  /* Can't be a subset if num_dds+1 > num_sdiffs. */
  if (num_dds >= num_sdiffs) {
    return -1;
  }

  if (DEBUG) {
    printf("ref_prn = %u\nnon_ref_prns = {", ref_prn);
    for (u8 i=0; i < num_dds; i++) {
      printf("%d, ", non_ref_prns[i]);
    }
    printf("}\nnum_dds = %u\nnum_sdiffs = %u\nsdiffs[*].prn = {", num_dds, num_sdiffs);
    for (u8 i=0; i < num_sdiffs; i++) {
      printf("%d, ", sdiffs_in[i].sid.prn);
    }
    printf("}\n");
  }

  if (!is_prn_set(num_dds, non_ref_prns)) {
    log_error("There is disorder in the amb_test sats.");
    printf("amb_test sat prns = {%u, ", ref_prn);
    for (u8 k=0; k < num_dds; k++) {
      printf("%u, ", non_ref_prns[k]);
    }
    printf("}\n");
    DEBUG_EXIT();
    return -2;
  }

  double ref_phase = 0;
  double ref_pseudorange = 0;
  u8 i=0;
  u8 j=0;
  u8 found_ref = 0;
  /* Go through the sdiffs, pulling out the measurements of the non-ref amb sats
   * and the reference sat. */
  while (i < num_dds) {
    if (non_ref_prns[i] == sdiffs_in[j].sid.prn) {
      /* When we find a non-ref sat, we fill in the next measurement. */
      memcpy(&sdiffs_out[i+1], &sdiffs_in[j], sizeof(sdiff_t));
      dd_meas[i] = sdiffs_in[j].carrier_phase;
      dd_meas[i+num_dds] = sdiffs_in[j].pseudorange;
      i++;
      j++;
    } else if (ref_prn == sdiffs_in[j].sid.prn) {
      /* when we find the ref sat, we copy it over and raise the FOUND flag */
      memcpy(&sdiffs_out[0], &sdiffs_in[j], sizeof(sdiff_t));
      ref_phase =  sdiffs_in[j].carrier_phase;
      ref_pseudorange = sdiffs_in[j].pseudorange;
      j++;
      found_ref = 1;
    }
    else if (non_ref_prns[i] > sdiffs_in[j].sid.prn) {
      /* If both sets are ordered, and we increase j (and possibly i), and the
       * i prn is higher than the j one, it means that the i one might be in the
       * j set for higher j, and that the current j prn isn't in the i set. */
      j++;
    } else {
      /* if both sets are ordered, and we increase j (and possibly i), and the
       * j prn is higher than the i one, it means that the j one might be in the
       * i set for higher i, and that the current i prn isn't in the j set.
       * This means a sat in the IAR's sdiffs isn't in the sdiffs.
       * */
      DEBUG_EXIT();
      return -1;
    }
  }
  /* This awkward case deals with the situation when sdiffs and sats have the
   * same satellites only the ref of amb_test.sats is the last PRN in sdiffs.
   * This case is never checked for j = num_dds as i only runs to num_dds-1. */
  /* TODO: This function could be refactored to be a lot clearer. */
  while (!found_ref && j < num_sdiffs ) {
    if (ref_prn == sdiffs_in[j].sid.prn) {
      memcpy(&sdiffs_out[0], &sdiffs_in[j], sizeof(sdiff_t));
      ref_phase =  sdiffs_in[j].carrier_phase;
      ref_pseudorange = sdiffs_in[j].pseudorange;
      found_ref = 1;
    }
    j++;
  }

  if (found_ref == 0) {
    DEBUG_EXIT();
    return -1;
  }
  for (i=0; i < num_dds; i++) {
    dd_meas[i] -= ref_phase;
    dd_meas[i+num_dds] -= ref_pseudorange;
  }
  if (DEBUG) {
    printf("amb_sdiff_prns = {");
    for (i = 0; i < num_dds+1; i++) {
      printf("%u, ", sdiffs_out[i].sid.prn);
    }
    printf("}\ndd_measurements = {");
    for (i=0; i < 2 * num_dds; i++) {
      printf("%f, \t", dd_meas[i]);
    }
    printf("}\n");
  }

  DEBUG_EXIT();
  return 0;
}

s8 copy_sdiffs_put_ref_first(const u8 ref_prn, const u8 num_sdiffs, const sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  s8 not_found = -1;
  u8 j = 1;
  for (u8 i=0; i<num_sdiffs; i++) {
    if (sdiffs[i].sid.prn == ref_prn) {
      memcpy(sdiffs_with_ref_first, &sdiffs[i], sizeof(sdiff_t));
      not_found = 0;
    }
    else {
      if (j == num_sdiffs) {
        return not_found; //note: not_found should always be -1 if it reaches this point. if anyone thinks "return -1" is more clear, go ahead and change it.
      }
      memcpy(&sdiffs_with_ref_first[j], &sdiffs[i], sizeof(sdiff_t));
      j++;
    }
  }
  return not_found;
}

static bool contains_prn(u8 len, u8 *prns, u8 prn)
{
  for (u8 i = 0; i < len; i++) {
    if (prns[i] == prn) {
      return true;
    }
  }
  return false;
}

u8 filter_sdiffs(u8 num_sdiffs, sdiff_t *sdiffs, u8 num_sats_to_drop, u8 *sats_to_drop)
{
  u8 new_num_sdiffs = 0;
  for (u8 i = 0; i < num_sdiffs; i++) {
    if (!contains_prn(num_sats_to_drop, sats_to_drop, sdiffs[i].sid.prn)) {
      if (new_num_sdiffs != i) {
        memcpy(&sdiffs[new_num_sdiffs], &sdiffs[i], sizeof(sdiff_t));
      }
      new_num_sdiffs++;
    }
  }
  return new_num_sdiffs;
}

/** Prints an sdiff_t
 * \param sd    the sdiff_t to print.
 */
void debug_sdiff(sdiff_t sd)
{
  log_debug("sdiff_t:"
    "\tprn = %u\n"
    "\tsnr = %f\n"
    "\tpseudorange   = %f\n"
    "\tcarrier_phase = %f\n"
    "\tdoppler       = %f\n"
    "\tsat_pos = [%f, %f, %f]\n"
    "\tsat_vel = [%f, %f, %f]\n",
    sd.sid.prn, sd.snr,
    sd.pseudorange, sd.carrier_phase, sd.doppler,
    sd.sat_pos[0], sd.sat_pos[1], sd.sat_pos[2],
    sd.sat_vel[0], sd.sat_vel[1], sd.sat_vel[2]);
}

/** Prints an array of sdiffs
 * \param n     The number of sdiffs to print.
 * \param sds   A pointer to the head of the array of sdiffs to print.
 */
void debug_sdiffs(u8 n, sdiff_t *sds)
{
  log_debug("[");
  for (u8 i=0; i<n; i++) {
    debug_sdiff(sds[i]);
  }
  log_debug("]");
}

/** \} */

