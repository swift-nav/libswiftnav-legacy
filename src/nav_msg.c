/*
 * Copyright (C) 2010 Swift Navigation Inc.
 * Contact: Henry Hallam <henry@swift-nav.com>
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
#include <stdlib.h>
#include <math.h>

#include "logging.h"
#include "constants.h"
#include "bits.h"
#include "nav_msg.h"

/* Approx number of nav bit edges needed to accept bit sync for a
   strong signal (sync will take longer on a weak signal) */
#define BITSYNC_THRES 22

void nav_msg_init(nav_msg_t *n)
{
  /* Initialize the necessary parts of the nav message state structure. */
  memset(n, 0, sizeof(nav_msg_t));
  n->bit_phase_ref = BITSYNC_UNSYNCED;
  n->next_subframe_id = 1;
  n->bit_polarity = BIT_POLARITY_UNKNOWN;
}

static u32 extract_word(nav_msg_t *n, u16 bit_index, u8 n_bits, u8 invert)
{
  /* Extract a word of n_bits length (n_bits <= 32) at position bit_index into
   * the subframe. Takes account of the offset stored in n, and the circular
   * nature of the n->subframe_bits buffer. */

  /* Offset for the start of the subframe in the buffer. */
  if (n->subframe_start_index) {
    if (n->subframe_start_index > 0) {
      bit_index += n->subframe_start_index; /* Standard. */
    } else {
      bit_index -= n->subframe_start_index; /* Bits are inverse! */
      invert = !invert;
    }

    bit_index--;
  }

  /* Wrap if necessary. */
  if (bit_index > NAV_MSG_SUBFRAME_BITS_LEN * 32) {
    bit_index -= NAV_MSG_SUBFRAME_BITS_LEN * 32;
  }

  u8 bix_hi = bit_index >> 5;
  u8 bix_lo = bit_index & 0x1F;
  u32 word = n->subframe_bits[bix_hi] << bix_lo;

  if (bix_lo) {
    bix_hi++;
    if (bix_hi == NAV_MSG_SUBFRAME_BITS_LEN) {
      bix_hi = 0;
    }
    word |=  n->subframe_bits[bix_hi] >> (32 - bix_lo);
  }

  if (invert) {
    word = ~word;
  }

  return word >> (32 - n_bits);
}

/* TODO: Bit synchronization that can operate with multi-ms integration times
   e.g. http://www.thinkmind.org/download.php?articleid=spacomm_2013_2_30_30070
 */
static void update_bit_sync(nav_msg_t *n, s32 corr_prompt_real)
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
  n->bit_integrate -= n->bitsync_prev_corr[n->bit_phase];
  n->bitsync_prev_corr[n->bit_phase] = corr_prompt_real;
  if (n->bitsync_count < 20) {
    n->bitsync_count++;
    return;  /* That rolling accumulator is not valid yet */
  }

  /* Add the accumulator to the histogram for the relevant phase */
  n->bitsync_histogram[(n->bit_phase) % 20] += abs(n->bit_integrate);

  if (n->bit_phase == 20 - 1) {
    /* Histogram is valid.  Find the two highest values. */
    u32 max = 0, next_best = 0;
    u32 max_prev_corr = 0;
    u8 max_i = 0;
    for (u8 i = 0; i < 20; i++) {
      u32 v = n->bitsync_histogram[i];
      if (v > max) {
        next_best = max;
        max = v;
        max_i = i;
      } else if (v > next_best) {
        next_best = v;
      }
      /* Also find the highest value from the last 20 correlations.
         We'll use this to normalize the threshold score. */
      v = abs(n->bitsync_prev_corr[i]);
      if (v > max_prev_corr) {
        max_prev_corr = v;
      }
    }
    /* Form score from difference between the best and the second-best */
    if (max - next_best > BITSYNC_THRES * 2 * max_prev_corr) {
      /* We are synchronized! */
      n->bit_phase_ref = max_i;
      /* TODO: Subtract necessary older prev_corrs from bit_integrate to
         ensure it will be correct for the upcoming first dump */
    }
  }
}

/** Navigation message decoding update.
 * Called once per tracking loop update. Performs the necessary steps to
 * recover the nav bit clock, store the nav bits and decode them.
 *
 * Also extracts and returns the GPS time of week each time a new subframe is
 * received.
 *
 * \param n Nav message decode state struct
 * \param corr_prompt_real In-phase prompt correlation from tracking loop
 * \param ms Number of milliseconds integration performed in the correlation
 *
 * \return The GPS time of week in milliseconds of the current code phase
 *         rollover, or `TOW_INVALID` (-1) if unknown
 */
s32 nav_msg_update(nav_msg_t *n, s32 corr_prompt_real, u8 ms)
{
  s32 TOW_ms = TOW_INVALID;

  n->bit_phase += ms;
  n->bit_phase %= 20;
  n->bit_integrate += corr_prompt_real;
  /* Do we have bit phase lock yet? (Do we know which of the 20 possible PRN
   * offsets corresponds to the nav bit edges?) */
  if (n->bit_phase_ref == BITSYNC_UNSYNCED) {
    update_bit_sync(n, corr_prompt_real);
  }

  if (n->bit_phase != n->bit_phase_ref) {
    /* Either we don't have bit phase lock, or this particular
       integration is not aligned to a nav bit boundary. */
    return TOW_INVALID;
  }

  /* Dump the nav bit, i.e. determine the sign of the correlation over the
   * nav bit period. */
  bool bit_val = n->bit_integrate > 0;
  /* Zero the accumulator for the next nav bit. */
  n->bit_integrate = 0;

  /* The following is offset by 27 to allow the handover word to be
   * overwritten.  This is not a problem as it's handled below and
   * isn't needed by the ephemeris decoder.
   */
  u16 last_subframe_bit_index = ABS(n->subframe_start_index) + 27;
  last_subframe_bit_index %= NAV_MSG_SUBFRAME_BITS_LEN * 32;
  if (n->subframe_start_index &&
      (n->subframe_bit_index == last_subframe_bit_index)) {
    /* Subframe buffer is full: the nav message decoder has missed it's
     * deadline.  Clobbering the buffer can result in invalid nav data
     * being used.
     */
    n->overrun = true;
    return -2;
  }

  if (n->subframe_bit_index >= NAV_MSG_SUBFRAME_BITS_LEN * 32) {
    log_error("subframe bit index gone wild %d", (int)n->subframe_bit_index);
    return -22;
  }

  if (bit_val) {
    n->subframe_bits[n->subframe_bit_index >> 5] |= \
      1 << (31 - (n->subframe_bit_index & 0x1F));
  } else {
    /* Integrated correlation is negative, so bit is 0. */
    n->subframe_bits[n->subframe_bit_index >> 5] &= \
      ~(1 << (31 - (n->subframe_bit_index & 0x1F)));
  }

  n->subframe_bit_index++;
  if (n->subframe_bit_index == NAV_MSG_SUBFRAME_BITS_LEN * 32) {
    n->subframe_bit_index = 0;
  }

  /* Yo dawg, are we still looking for the preamble? */
  if (!n->subframe_start_index) {
    /* We're going to look for the preamble at a time 360 nav bits ago,
     * then again 60 nav bits ago. */
    #define SUBFRAME_START_BUFFER_OFFSET (NAV_MSG_SUBFRAME_BITS_LEN * 32 - 360)

    /* Check whether there's a preamble at the start of the circular
     * subframe_bits buffer. */
    u8 preamble_candidate = extract_word(n,
                                         n->subframe_bit_index + SUBFRAME_START_BUFFER_OFFSET, 8,
                                         0);

    if (preamble_candidate == 0x8B) {
      n->subframe_start_index = n->subframe_bit_index +
                                SUBFRAME_START_BUFFER_OFFSET + 1;
    } else if (preamble_candidate == 0x74) {
      n->subframe_start_index =
        -(n->subframe_bit_index + SUBFRAME_START_BUFFER_OFFSET + 1);
    }

    if (n->subframe_start_index) {
      /* Looks like we found a preamble, but let's confirm. */
      if (extract_word(n, 300, 8, 0) == 0x8B) {
        /* There's another preamble in the following subframe.  Looks good so far. */
        /* Extract the TOW: */
        unsigned int TOW_trunc =
          extract_word(n, 30, 17, extract_word(n, 29, 1, 0));
        /* (bit 29 is D30* for the second word, where the TOW resides) */
        if (TOW_trunc < 7 * 24 * 60 * 10) {
          /* TOW in valid range */
          TOW_trunc++;                         /* Increment it, to see what we expect at the start of the next subframe */
          if (TOW_trunc == 7 * 24 * 60 * 10) { /* Handle end of week rollover */
            TOW_trunc = 0;
          }

          if (TOW_trunc ==
              extract_word(n, 330, 17, extract_word(n, 329, 1, 0))) {
            /* We got two appropriately spaced preambles, and two matching TOW counts.  Pretty certain now. */
            /* TODO: should still check parity? */
            /* The TOW in the message is for the start of the NEXT subframe. */
            /* That is, 240 nav bits' time from now, since we are 60 nav bits into the second subframe that we recorded. */
            if (TOW_trunc == 0) {
              /* end-of-week special case */
              TOW_ms = 7 * 24 * 60 * 60 * 1000 - (300 - 60) * 20;
            } else {
              TOW_ms = TOW_trunc * 6000 - (300 - 60) * 20;
            }
          }
        }
      }
      /* If we didn't find a matching pair of preambles + TOWs, this offset can't be right. Move on. */
      if (TOW_ms < 0) {
        n->subframe_start_index = 0;
      }
    }
  }
  return TOW_ms;
}

/* Tests the parity of a L1 C/A NAV message word.
 * Inverts the data bits if necessary, and checks the parity.
 * Expects a word where MSB = D29*, bit 30 = D30*, bit 29 = D1, ... LSB = D30.
 *
 * \note This function may modify the value of `word`.
 *
 * References:
 *   -# ICD-GPS-200E Table 20-XIV
 *
 * \param word Pointer to word to check. Note, if D30* is set then the data
 *             bits in this word will be inverted in place.
 * \return 0 if the parity is correct,
 *         otherwise returns the number of the first incorrect parity bit.
 */
static u8 nav_parity(u32 *word)
{
  if (*word & 1 << 30) { /* Inspect D30* */
    *word ^= 0x3FFFFFC0; /* D30* = 1, invert all the data bits! */
  }

  /* Check D25 */
  if (parity(*word & 0xBB1F34A0 /* 0b10111011000111110011010010100000 */)) {
    return 25;
  }
  /* Check D26 */
  if (parity(*word & 0x5D8F9A50 /* 0b01011101100011111001101001010000 */)) {
    return 26;
  }
  /* Check D27 */
  if (parity(*word & 0xAEC7CD08 /* 0b10101110110001111100110100001000 */)) {
    return 27;
  }
  /* Check D28 */
  if (parity(*word & 0x5763E684 /* 0b01010111011000111110011010000100 */)) {
    return 28;
  }
  /* Check D29 */
  if (parity(*word & 0x6BB1F342 /* 0b01101011101100011111001101000010 */)) {
    return 29;
  }
  /* Check D30 */
  if (parity(*word & 0x8B7A89C1 /* 0b10001011011110101000100111000001 */)) {
    return 30;
  }

  return 0;
}

bool subframe_ready(nav_msg_t *n)
{
  return n->subframe_start_index != 0;
}

s8 process_subframe(nav_msg_t *n, ephemeris_t *e)
{
  /* Check parity and parse out the ephemeris from the most recently received subframe */

  if (!e) {
    log_error("process_subframe: CALLED WITH e = NULL!");
    n->subframe_start_index = 0;  /* Mark the subframe as processed */
    n->next_subframe_id = 1;      /* Make sure we start again next time */
    return -1;
  }

  /* First things first - check the parity, and invert bits if necessary. */
  /* process the data, skipping the first word, TLM, and starting with HOW */

  /* Detect half cycle slip. */
  s8 prev_bit_polarity = n->bit_polarity;
  n->bit_polarity = (n->subframe_start_index > 0) ? BIT_POLARITY_NORMAL :
                    BIT_POLARITY_INVERTED;
  if ((prev_bit_polarity != BIT_POLARITY_UNKNOWN)
      && (prev_bit_polarity != n->bit_polarity)) {
    log_error("PRN %02d Nav phase flip - half cycle slip detected, "
              "but not corrected", e->prn + 1);
    /* TODO: declare phase ambiguity to IAR */
  }

  /* Complain if buffer overrun */
  if (n->overrun) {
    log_error("PRN %02d nav_msg subframe buffer overrun!", e->prn + 1);
    n->overrun = false;
  }

  /* Extract word 2, and the last two parity bits of word 1 */
  u32 sf_word2 = extract_word(n, 28, 32, 0);
  if (nav_parity(&sf_word2)) {
    log_info("PRN %02d subframe parity mismatch (word 2)", e->prn + 1);
    n->subframe_start_index = 0;  /* Mark the subframe as processed */
    n->next_subframe_id = 1;      /* Make sure we start again next time */
    return -2;
  }

  u8 sf_id = sf_word2 >> 8 & 0x07;                                             /* Which of 5 possible subframes is it? */

  if (sf_id <= 3 && sf_id == n->next_subframe_id) {                            /* Is it the one that we want next? */

    for (int w = 0; w < 8; w++) {                                              /* For words 3..10 */
      n->frame_words[sf_id - 1][w] = extract_word(n, 30 * (w + 2) - 2, 32, 0); /* Get the bits */
      /* MSBs are D29* and D30*.  LSBs are D1...D30 */
      if (nav_parity(&n->frame_words[sf_id - 1][w])) {                         /* Check parity and invert bits if D30* */
        log_info("PRN %02d subframe parity mismatch (word %d)", e->prn + 1,
                 w + 3);
        n->next_subframe_id = 1;      /* Make sure we start again next time */
        n->subframe_start_index = 0;  /* Mark the subframe as processed */
        return -3;
      }
    }
    n->subframe_start_index = 0;  /* Mark the subframe as processed */
    n->next_subframe_id++;

    if (sf_id == 3) {
      /* Got all of subframes 1 to 3 */
      n->next_subframe_id = 1;      /* Make sure we start again next time */

      /* Now let's actually decode the ephemeris... */
      decode_ephemeris(n->frame_words, e);

      return 1;

    }
  } else {                       /* didn't get the subframe that we want next */
    n->next_subframe_id = 1;     /* Make sure we start again next time */
    n->subframe_start_index = 0; /* Mark the subframe as processed */
  }

  return 0;

}

