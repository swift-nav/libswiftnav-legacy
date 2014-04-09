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
#include <stdio.h>

#include "linear_algebra.h"
#include "single_diff.h"
#include "constants.h"

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

/** \} */

