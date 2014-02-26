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

#include "single_diff.h"

/** \defgroup single_diff Single Difference Observations
 * Functions for storing and manipulating single difference observations.
 * \{ */

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
  }
}

/** \} */

