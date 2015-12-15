/*
 * Copyright (c) 2015 Swift Navigation Inc.
 * Contact: Jacob McNamee <jacob@swiftnav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#include <math.h>
#include <assert.h>
#include <string.h>

#include <libswiftnav/bit_sync.h>

/* Approx number of nav bit edges needed to accept bit sync for a
   strong signal (sync will take longer on a weak signal) */
#define BITSYNC_THRES 22

/* Bit lengths for different constellations. Bounded by BIT_LENGTH_MAX */
#define BIT_LENGTH_GPS_L1 20
#define BIT_LENGTH_SBAS_L1 2

static void histogram_update(bit_sync_t *b, s32 corr_prompt_real);

/** \defgroup bit_sync Bit Sync
 * Functions and calculations related to data bit synchronization
 *
 * \{ */

/** Initialize a bit sync structure
 *
 * \param b    Pointer to a bit sync structure
 * \param sid  Signal identifier
 */
void bit_sync_init(bit_sync_t *b, gnss_signal_t sid)
{
  memset(b, 0, sizeof(bit_sync_t));
  b->bit_phase_ref = BITSYNC_UNSYNCED;

  assert(sid.band == BAND_L1);
  switch (sid.constellation) {
  case CONSTELLATION_GPS:
    b->bit_length = BIT_LENGTH_GPS_L1;
    break;
  case CONSTELLATION_SBAS:
    b->bit_length = BIT_LENGTH_SBAS_L1;
    break;
  default:
    assert("unsupported constellation");
  }
}

/** Update bit sync and get bit integration output
 *
 * \param b                 Pointer to a bit sync structure
 * \param corr_prompt_real  Real part of the prompt correlation
 * \param ms                Integration time (ms) of the correlation
 * \param bit_integrate     Pointer to output bit integration (if valid)
 *
 * \return  True if *bit_integrate contains a valid bit integration,
 *          False otherwise
 */
bool bit_sync_update(bit_sync_t *b, s32 corr_prompt_real, u32 ms,
                     s32 *bit_integrate)
{
  b->bit_phase += ms;
  b->bit_phase %= b->bit_length;
  b->bit_integrate += corr_prompt_real;

  /* Search for bit phase if not yet locked. */
  if (b->bit_phase_ref == BITSYNC_UNSYNCED)
    histogram_update(b, corr_prompt_real);

  /* Return the integration at the end of the bit period */
  if (b->bit_phase == b->bit_phase_ref) {
    *bit_integrate = b->bit_integrate;
    b->bit_integrate = 0;
    return true;
  }

  return false;
}

/* TODO: Bit synchronization that can operate with multi-ms integration times
   e.g. http://www.thinkmind.org/download.php?articleid=spacomm_2013_2_30_30070
 */
static void histogram_update(bit_sync_t *b, s32 corr_prompt_real)
{
  /* On 20th call:
     bit_phase = 0
     bitsync_count = 20
     bit_integrate holds sum of first 20 correlations
     bitsync_histogram is all zeros
     bitsync_prev_corr[0]=0, others hold previous correlations ([1] = first)
     In this function:
       bitsync_prev_corr[0] <= corr_prompt_real (20th correlation)
       bitsync_histogram[0] <= sum of corrs [0..19]
     On 21st call:
     bit_phase = 1
     bit_integrate holds sum of corrs [0..20]
       bit_integrate -= bitsync_prev_corr[1]
         bit_integrate now holds sum of corrs [1..20]
       bitsync_histogram[1] <= bit_integrate
     ...
     On 39th call:
     bit_phase = 19
00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40
____bit_integrate is sum of these after 20th call__________
                                                         __________________after 39th call__________________________

       after subtraction, bit_integrate holds sum of corrs [19..38]
       bitsync_histogram[19] <= bit_integrate. Now fully populated.
       Suppose first correlation happened to be first ms of a nav bit
        then max_i = bit_phase_ref = 0.
  */

  /* Maintain a rolling sum of the 20 most recent correlations in
     bit_integrate */
  b->bit_integrate -= b->bitsync_prev_corr[b->bit_phase];
  b->bitsync_prev_corr[b->bit_phase] = corr_prompt_real;
  if (b->bitsync_count < b->bit_length) {
    b->bitsync_count++;
    return;  /* That rolling accumulator is not valid yet */
  }

  /* Add the accumulator to the histogram for the relevant phase */
  b->bitsync_histogram[(b->bit_phase) % b->bit_length] += ABS(b->bit_integrate);

  if (b->bit_phase == b->bit_length - 1) {
    /* Histogram is valid.  Find the two highest values. */
    u32 max = 0, next_best = 0;
    u32 max_prev_corr = 0;
    u8 max_i = 0;
    for (u8 i = 0; i < b->bit_length; i++) {
      u32 v = b->bitsync_histogram[i];
      if (v > max) {
        next_best = max;
        max = v;
        max_i = i;
      } else if (v > next_best) {
        next_best = v;
      }
      /* Also find the highest value from the last 20 correlations.
         We'll use this to normalize the threshold score. */
      v = ABS(b->bitsync_prev_corr[i]);
      if (v > max_prev_corr)
        max_prev_corr = v;
    }
    /* Form score from difference between the best and the second-best */
    if (max - next_best > BITSYNC_THRES * 2 * max_prev_corr) {
      /* We are synchronized! */
      b->bit_phase_ref = max_i;
      /* TODO: Subtract necessary older prev_corrs from bit_integrate to
         ensure it will be correct for the upcoming first dump */
    }
  }
}

/** \} */

