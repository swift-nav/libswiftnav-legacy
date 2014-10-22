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

#include <stdio.h>
#include <string.h>

#include "linear_algebra.h"
#include "single_diff.h"
#include "ephemeris.h"
#include "constants.h"
#include "sats_management.h"

/** \defgroup single_diff Single Difference Observations
 * Functions for storing and manipulating single difference observations.
 * \{ */

#define GPS_L1_LAMBDA (GPS_C / GPS_L1_HZ)

u8 propagate(u8 n, double ref_ecef[3],
             navigation_measurement_t *m_in_base, gps_time_t *t_base,
             navigation_measurement_t *m_in_rover, gps_time_t *t_rover,
             navigation_measurement_t *m_out_base)
{
  double dt = gpsdifftime(*t_rover, *t_base);
  (void)dt;

  double dr_[n];

  for (u8 i=0; i<n; i++) {
    m_out_base[i].prn = m_in_base[i].prn;
    m_out_base[i].snr = m_in_base[i].snr;
    m_out_base[i].lock_time = m_in_base[i].lock_time;

    /* Calculate delta range. */
    double dr[3];
    vector_subtract(3, m_in_rover[i].sat_pos, m_in_base[i].sat_pos, dr);

    /* Subtract linear term (initial satellite velocity * dt),
     * we are going to add back in the linear term derived from the Doppler
     * instead. */
    /*vector_add_sc(3, dr, m_in_base[i].sat_vel, -dt, dr);*/

    /* Make unit vector to satellite, e. */
    double e[3];
    vector_subtract(3, m_in_rover[i].sat_pos, ref_ecef, e);
    vector_normalize(3, e);

    /* Project onto the line of sight vector. */
    dr_[i] = vector_dot(3, dr, e);

    /*printf("# ddr_ = %f\n", dr_[i]);*/

    /* Add back in linear term now using Doppler. */
    /*dr_[i] -= m_in_base[i].raw_doppler * dt * GPS_L1_LAMBDA;*/

    /*printf("# dr_dopp = %f\n", -m_in_base[i].raw_doppler * dt * GPS_L1_LAMBDA);*/

    /*printf("# raw dopp = %f\n", m_in_rover[i].raw_doppler);*/
    /*printf("# my dopp = %f\n", -vector_dot(3, e, m_in_rover[i].sat_vel) / GPS_L1_LAMBDA);*/
    /*printf("# ddopp = %f\n", m_in_rover[i].raw_doppler - vector_dot(3, e, m_in_rover[i].sat_vel) / GPS_L1_LAMBDA);*/
    /*printf("[%f, %f],", m_in_base[i].raw_doppler, vector_dot(3, e, m_in_rover[i].sat_vel) / GPS_L1_LAMBDA);*/

    /*printf("# dr_ = %f\n", dr_[i]);*/

    m_out_base[i].raw_pseudorange = m_in_base[i].raw_pseudorange + dr_[i];
    m_out_base[i].pseudorange = m_in_base[i].pseudorange;
    m_out_base[i].carrier_phase = m_in_base[i].carrier_phase - dr_[i] / GPS_L1_LAMBDA;
    m_out_base[i].raw_doppler = m_in_base[i].raw_doppler;
    m_out_base[i].doppler = m_in_base[i].doppler;
    /*m_in_base[i].carrier_phase -= dr_[i] / GPS_L1_LAMBDA;*/
  }
  return 0;
}

/** Calculate single difference observations.
 * Undifferenced input observations are assumed to be both taken at the
 * same time, `t`.
 *
 * SNR in the output is the lesser of the SNRs of inputs a and b.
 *
 * `sat_pos` and `sat_vel` are taken from input a.
 *
 * \param n_a Number of measurements in `m_a`
 * \oaram m_a Array of undifferenced observations, sorted by PRN
 * \param n_b Number of measurements in `m_b`
 * \oaram m_new Array of b navigation measurements, sorted by PRN
 * \param sds Single difference observations
 * \return The number of observations written to `sds`
 */
u8 single_diff(u8 n_a, navigation_measurement_t *m_a,
               u8 n_b, navigation_measurement_t *m_b,
               sdiff_t *sds)
{
  u8 i, j, n = 0;

  /* Loop over m_a and m_b and check if a PRN is present in both. */
  for (i=0, j=0; i<n_a && j<n_b; i++, j++) {
    if (m_a[i].prn < m_b[j].prn)
      j--;
    else if (m_a[i].prn > m_b[j].prn)
      i--;
    else {
      sds[n].prn = m_a[i].prn;
      sds[n].pseudorange = m_a[i].raw_pseudorange - m_b[j].raw_pseudorange;
      sds[n].carrier_phase = m_a[i].carrier_phase - m_b[j].carrier_phase;
      sds[n].doppler = m_a[i].raw_doppler - m_b[j].raw_doppler;
      sds[n].snr = MIN(m_a[i].snr, m_b[j].snr);

      memcpy(&(sds[n].sat_pos), &(m_a[i].sat_pos), 3*sizeof(double));
      memcpy(&(sds[n].sat_vel), &(m_a[i].sat_vel), 3*sizeof(double));

      n++;
    }
  }

  return n;
}

int sdiff_search_prn(const void *a, const void *b)
{
  return (*(u8*)a - ((sdiff_t *)b)->prn);
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
 * \param sds               The single differenced propagated measurements.
 * \return The number of sats common in both local and remote sdiffs.
 */
u8 make_propagated_sdiffs(u8 n_local, navigation_measurement_t *m_local,
                          u8 n_remote, navigation_measurement_t *m_remote,
                          double *remote_dists, double remote_pos_ecef[3],
                          ephemeris_t *es, gps_time_t t,
                          sdiff_t *sds)
{
  u8 i, j, n = 0;

  /* Loop over m_a and m_b and check if a PRN is present in both. */
  for (i=0, j=0; i<n_local && j<n_remote; i++, j++) {
    if (m_local[i].prn < m_remote[j].prn)
      j--;
    else if (m_local[i].prn > m_remote[j].prn)
      i--;
    else if (ephemeris_good(es[m_local[i].prn], t)) {
      double clock_err;
      double clock_rate_err;
      double local_sat_pos[3];
      double local_sat_vel[3];
      calc_sat_pos(&local_sat_pos[0],
                   &local_sat_vel[0],
                   &clock_err, &clock_rate_err, &es[m_local[i].prn], t);
      sds[n].prn = m_local[i].prn;
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


/** Convert a list of almanacs to a list of single differences.
 * This only fills the position, velocity and prn.
 *
 * It's only useful for using functions that need a sdiff_t when you have almanac_t.
 */
void almanacs_to_single_diffs(u8 n, almanac_t *alms, gps_time_t timestamp, sdiff_t *sdiffs)
{
  for (u8 i=0; i<n; i++) {
    double p[3];
    double v[3];
    calc_sat_state_almanac(&alms[i], timestamp.tow, timestamp.wn, p, v);
    memcpy(sdiffs[i].sat_pos, &p[0], 3 * sizeof(double));
    memcpy(sdiffs[i].sat_vel, &v[0], 3 * sizeof(double));
    sdiffs[i].prn = alms[i].prn;
    if (i==0) {
      sdiffs[i].snr = 1;
      // printf("ref_prn=%d\n", sdiffs[i].prn);
    }
    else {
      sdiffs[i].snr = 0;
    }
    // if (sdiffs[i].prn == 16) {
    //   sdiffs[i].snr = 2;
    // }
  }
}

void double_diff(u8 n, sdiff_t *sds, sdiff_t *dds, u8 ref_idx)
{
  for (u8 i=0; i<n; i++) {
    dds[i].prn = sds[i].prn;
    dds[i].pseudorange = sds[i].pseudorange - sds[ref_idx].pseudorange;
    dds[i].carrier_phase = sds[i].carrier_phase - sds[ref_idx].carrier_phase;
    dds[i].doppler = sds[i].doppler - sds[ref_idx].doppler;
    dds[i].snr = MIN(sds[i].snr, sds[ref_idx].snr);

    memcpy(&(dds[i].sat_pos), &(sds[i].sat_pos), 3*sizeof(double));
    memcpy(&(dds[i].sat_vel), &(sds[i].sat_vel), 3*sizeof(double));
  }
}

s8 copy_sdiffs_put_ref_first(u8 ref_prn, u8 num_sdiffs, sdiff_t *sdiffs, sdiff_t *sdiffs_with_ref_first)
{
  s8 not_found = -1;
  u8 j = 1;
  for (u8 i=0; i<num_sdiffs; i++) {
    if (sdiffs[i].prn == ref_prn) {
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

bool _contains_prn(u8 len, u8 *prns, u8 prn)
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
    if (!_contains_prn(num_sats_to_drop, sats_to_drop, sdiffs[i].prn)) {
      if (new_num_sdiffs != i) {
        memcpy(&sdiffs[new_num_sdiffs], &sdiffs[i], sizeof(sdiff_t));
      }
      new_num_sdiffs++;
    }
  }
  return new_num_sdiffs;
}

/** \} */

