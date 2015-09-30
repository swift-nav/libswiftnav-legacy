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
#include <assert.h>

#include "logging.h"
#include "constants.h"
#include "bits.h"
#include "nav_msg.h"
#include "edc.h"

/* Approx number of nav bit edges needed to accept bit sync for a
   strong signal (sync will take longer on a weak signal) */
#define L1_BITSYNC_THRES 22
#define SBAS_BITSYNC_THRES 200

static void l1_legacy_nav_msg_init(l1_legacy_nav_msg_t *n)
{
  /* Initialize the necessary parts of the nav message state structure. */
  memset(n, 0, sizeof(l1_legacy_nav_msg_t));
  n->bit_phase_ref = BITSYNC_UNSYNCED;
  n->next_subframe_id = 1;
  n->bit_polarity = BIT_POLARITY_UNKNOWN;
  n->bit_length = 20;
}

static void l1_sbas_nav_msg_init(l1_sbas_nav_msg_t *n)
{
  /* Initialize the necessary parts of the nav message state structure. */
  memset(n, 0, sizeof(l1_sbas_nav_msg_t));
  n->bit_phase_ref = BITSYNC_UNSYNCED;
  n->bit_polarity = BIT_POLARITY_UNKNOWN;
  n->bit_length = 2;

  n->dec_passes = 0;

  n->polarity = 0;

  n->msg_normal = 0;
  n->msg_inverse = 0;
  n->init = 1;
}

void nav_msg_init(nav_msg_t *n)
{
  switch(n->type) {
    case L1_LEGACY_NAV:
      l1_legacy_nav_msg_init(n->l1_nav_msg);
      break;
    case L1_SBAS:
      l1_sbas_nav_msg_init(n->sbas_nav_msg);
  }
}

static u32 extract_word(l1_legacy_nav_msg_t *n, u16 bit_index, u8 n_bits, u8 invert)
{
  /* Extract a word of n_bits length (n_bits <= 32) at position bit_index into
   * the subframe. Takes account of the offset stored in n, and the circular
   * nature of the n->subframe_bits buffer. */

  /* Offset for the start of the subframe in the buffer. */
  if (n->subframe_start_index) {
    if (n->subframe_start_index > 0)
      bit_index += n->subframe_start_index; /* Standard. */
    else {
      bit_index -= n->subframe_start_index; /* Bits are inverse! */
      invert = !invert;
    }

    bit_index--;
  }

  /* Wrap if necessary. */
  if (bit_index > L1_NAV_MSG_SUBFRAME_BITS_LEN*32)
    bit_index -= L1_NAV_MSG_SUBFRAME_BITS_LEN*32;

  u8 bix_hi = bit_index >> 5;
  u8 bix_lo = bit_index & 0x1F;
  u32 word = n->subframe_bits[bix_hi] << bix_lo;

  if (bix_lo) {
    bix_hi++;
    if (bix_hi == L1_NAV_MSG_SUBFRAME_BITS_LEN)
      bix_hi = 0;
    word |=  n->subframe_bits[bix_hi] >> (32 - bix_lo);
  }

  if (invert)
    word = ~word;

  return word >> (32 - n_bits);
}


/* TODO: Bit synchronization that can operate with multi-ms integration times
   e.g. http://www.thinkmind.org/download.php?articleid=spacomm_2013_2_30_30070
 */
static void update_bit_sync(nav_msg_t *n, s32 corr_prompt_real, u8 bit_len)
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
  u8 thres;
  s32 *bit_integrate;
  u8 *bit_phase, *bitsync_count;
  s8 *bit_phase_ref;
  s32 *bitsync_prev_corr;
  u32 *bitsync_histogram;
  switch (n->type) {
    case L1_LEGACY_NAV:
      thres = L1_BITSYNC_THRES;
      bit_integrate = &n->l1_nav_msg->bit_integrate;
      bitsync_count = &n->l1_nav_msg->bitsync_count;
      bit_phase = &n->l1_nav_msg->bit_phase;
      bit_phase_ref = &n->l1_nav_msg->bit_phase_ref;
      bitsync_histogram = n->l1_nav_msg->bitsync_histogram;
      bitsync_prev_corr = n->l1_nav_msg->bitsync_prev_corr;
      break;
    case L1_SBAS:
      thres = SBAS_BITSYNC_THRES;
      bit_integrate = &n->sbas_nav_msg->bit_integrate;
      bitsync_count = &n->sbas_nav_msg->bitsync_count;
      bit_phase = &n->sbas_nav_msg->bit_phase;
      bit_phase_ref = &n->sbas_nav_msg->bit_phase_ref;
      bitsync_histogram = n->sbas_nav_msg->bitsync_histogram;
      bitsync_prev_corr = n->sbas_nav_msg->bitsync_prev_corr;
      break;
    default:
      return;
  }
  *bit_integrate -= bitsync_prev_corr[*bit_phase];
  bitsync_prev_corr[*bit_phase] = corr_prompt_real;
  if (*bitsync_count < bit_len) {
    *bitsync_count = *bitsync_count + 1;
    return;  /* That rolling accumulator is not valid yet */
  }

  /* Add the accumulator to the histogram for the relevant phase */
  bitsync_histogram[*bit_phase % bit_len] += abs(*bit_integrate);

  if (*bit_phase == bit_len - 1) {
    /* Histogram is valid.  Find the two highest values. */
    u32 max = 0, next_best = 0;
    u32 max_prev_corr = 0;
    u8 max_i = 0;
    for (u8 i = 0; i < bit_len; i++) {
      u32 v = bitsync_histogram[i];
      if (v > max) {
        next_best = max;
        max = v;
        max_i = i;
      } else if (v > next_best) {
        next_best = v;
      }
      /* Also find the highest value from the last 20 correlations.
         We'll use this to normalize the threshold score. */
      v = abs(bitsync_prev_corr[i]);
      if (v > max_prev_corr)
        max_prev_corr =  v;
    }
    /* Form score from difference between the best and the second-best */
    if (max - next_best > thres * 2 * max_prev_corr) {
      /* We are synchronized! */
      *bit_phase_ref = max_i;
      /* TODO: Subtract necessary older prev_corrs from bit_integrate to
         ensure it will be correct for the upcoming first dump */
    }
  }
}

static s32 l1_sbas_nav_msg_update(l1_sbas_nav_msg_t *n, s32 corr_prompt_real, u8 ms)
{
  n->bit_phase += ms;
  n->bit_phase %= n->bit_length;
  n->bit_integrate += corr_prompt_real;
  /* Do we have bit phase lock yet? (Do we know which of the 20 possible PRN
   * offsets corresponds to the nav bit edges?) */
  if (n->bit_phase_ref == BITSYNC_UNSYNCED) {
    nav_msg_t nav = {
      .sbas_nav_msg = n,
      .type = L1_SBAS
    };
    update_bit_sync(&nav, corr_prompt_real, n->bit_length);
  }

  if (n->bit_phase != n->bit_phase_ref) {
    /* Either we don't have bit phase lock, or this particular
       integration is not aligned to a nav bit boundary. */
    n->good_bit = false;
    return TOW_INVALID;
  }

  /* Dump the nav bit, i.e. determine the sign of the correlation over the
   * nav bit period. */
  bool bit_val = n->bit_integrate > 0;
  /* Zero the accumulator for the next nav bit. */
  n->bit_integrate = 0;
  n->good_bit = true;

  if (bit_val) {
    n->last_bit = 255;
  } else {
    n->last_bit = 0;
  }

  return TOW_INVALID;
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
static s32 l1_legacy_nav_msg_update(l1_legacy_nav_msg_t *n, s32 corr_prompt_real, u8 ms)
{
  s32 TOW_ms = TOW_INVALID;

  n->bit_phase += ms;
  n->bit_phase %= n->bit_length;
  n->bit_integrate += corr_prompt_real;
  /* Do we have bit phase lock yet? (Do we know which of the 20 possible PRN
   * offsets corresponds to the nav bit edges?) */
  if (n->bit_phase_ref == BITSYNC_UNSYNCED) {
    nav_msg_t nav = {
      .l1_nav_msg = n,
      .type = L1_LEGACY_NAV
    };
    update_bit_sync(&nav, corr_prompt_real, n->bit_length);
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
  last_subframe_bit_index %= L1_NAV_MSG_SUBFRAME_BITS_LEN * 32;
  if (n->subframe_start_index &&
      (n->subframe_bit_index == last_subframe_bit_index)) {
    /* Subframe buffer is full: the nav message decoder has missed it's
     * deadline.  Clobbering the buffer can result in invalid nav data
     * being used.
     */
    n->overrun = true;
    return -2;
  }

  if (n->subframe_bit_index >= L1_NAV_MSG_SUBFRAME_BITS_LEN*32) {
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
  if (n->subframe_bit_index == L1_NAV_MSG_SUBFRAME_BITS_LEN*32)
    n->subframe_bit_index = 0;

  /* Yo dawg, are we still looking for the preamble? */
  if (!n->subframe_start_index) {
    /* We're going to look for the preamble at a time 360 nav bits ago,
     * then again 60 nav bits ago. */
    #define SUBFRAME_START_BUFFER_OFFSET (L1_NAV_MSG_SUBFRAME_BITS_LEN*32 - 360)

    /* Check whether there's a preamble at the start of the circular
     * subframe_bits buffer. */
    u8 preamble_candidate = extract_word(n, n->subframe_bit_index + SUBFRAME_START_BUFFER_OFFSET, 8, 0);

    if (preamble_candidate == 0x8B) {
       n->subframe_start_index = n->subframe_bit_index + SUBFRAME_START_BUFFER_OFFSET + 1;
    }
    else if (preamble_candidate == 0x74) {
       n->subframe_start_index = -(n->subframe_bit_index + SUBFRAME_START_BUFFER_OFFSET + 1);
    }

    if (n->subframe_start_index) {
      // Looks like we found a preamble, but let's confirm.
      if (extract_word(n, 300, 8, 0) == 0x8B) {
        // There's another preamble in the following subframe.  Looks good so far.
        // Extract the TOW:
        unsigned int TOW_trunc = extract_word(n,30,17,extract_word(n,29,1,0));
        /* (bit 29 is D30* for the second word, where the TOW resides) */
        if (TOW_trunc < 7*24*60*10) {
          /* TOW in valid range */
          TOW_trunc++;  // Increment it, to see what we expect at the start of the next subframe
          if (TOW_trunc == 7*24*60*10)  // Handle end of week rollover
            TOW_trunc = 0;

          if (TOW_trunc == extract_word(n,330,17,extract_word(n,329,1,0))) {
            // We got two appropriately spaced preambles, and two matching TOW counts.  Pretty certain now.
            /* TODO: should still check parity? */
            // The TOW in the message is for the start of the NEXT subframe.
            // That is, 240 nav bits' time from now, since we are 60 nav bits into the second subframe that we recorded.
            if (TOW_trunc == 0)
              /* end-of-week special case */
              TOW_ms = 7*24*60*60*1000 - (300-60)*20;
            else
              TOW_ms = TOW_trunc * 6000 - (300-60)*20;
          }
        }
      }
      /* If we didn't find a matching pair of preambles + TOWs, this offset can't be right. Move on. */
      if (TOW_ms < 0)
        n->subframe_start_index = 0;
    }
  }
  return TOW_ms;
}


s32 nav_msg_update(nav_msg_t *n, s32 corr_prompt_real, u8 ms)
{
  switch (n->type) {
    case L1_LEGACY_NAV:
      return l1_legacy_nav_msg_update(n->l1_nav_msg, corr_prompt_real, ms);
      break;
    case L1_SBAS:
      return l1_sbas_nav_msg_update(n->sbas_nav_msg, corr_prompt_real, ms);
  }

  return TOW_INVALID;
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
  if (*word & 1<<30) { /* Inspect D30* */
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

bool subframe_ready(l1_legacy_nav_msg_t *n) {
  return (n->subframe_start_index != 0);
}

s8 l1_legacy_process_subframe(l1_legacy_nav_msg_t *n, ephemeris_kepler_t *e)
{
  // Check parity and parse out the ephemeris from the most recently received subframe

  if (!e) {
    log_error("process_subframe: CALLED WITH NULL ephemeris!");
    n->subframe_start_index = 0;  // Mark the subframe as processed
    n->next_subframe_id = 1;      // Make sure we start again next time
    return -1;
  }

  // First things first - check the parity, and invert bits if necessary.
  // process the data, skipping the first word, TLM, and starting with HOW

  /* Detect half cycle slip. */
  s8 prev_bit_polarity = n->bit_polarity;
  n->bit_polarity = (n->subframe_start_index > 0) ? BIT_POLARITY_NORMAL :
    BIT_POLARITY_INVERTED;
  if ((prev_bit_polarity != BIT_POLARITY_UNKNOWN)
      && (prev_bit_polarity != n->bit_polarity)) {
    log_error("PRN %02d Nav phase flip - half cycle slip detected, "
              "but not corrected", e->sid.prn+1);
    /* TODO: declare phase ambiguity to IAR */
  }

  /* Complain if buffer overrun */
  if (n->overrun) {
    log_error("PRN %02d nav_msg subframe buffer overrun!", e->sid.prn + 1);
    n->overrun = false;
  }

  /* Extract word 2, and the last two parity bits of word 1 */
  u32 sf_word2 = extract_word(n, 28, 32, 0);
  if (nav_parity(&sf_word2)) {
    log_info("PRN %02d subframe parity mismatch (word 2)", e->sid.prn + 1);
    n->subframe_start_index = 0;  // Mark the subframe as processed
    n->next_subframe_id = 1;      // Make sure we start again next time
    return -2;
  }

  u8 sf_id = sf_word2 >> 8 & 0x07;    // Which of 5 possible subframes is it?

  if (sf_id <= 3 && sf_id == n->next_subframe_id) {  // Is it the one that we want next?

    for (u16 w = 0; w < 8; w++) {   // For words 3..10
      n->frame_words[sf_id-1][w] = extract_word(n, 30*(w+2) - 2, 32, 0);    // Get the bits
      // MSBs are D29* and D30*.  LSBs are D1...D30
      if (nav_parity(&n->frame_words[sf_id-1][w])) {  // Check parity and invert bits if D30*
        log_info("PRN %02d subframe parity mismatch (word %d)", e->sid.prn + 1, w+3);
        n->next_subframe_id = 1;      // Make sure we start again next time
        n->subframe_start_index = 0;  // Mark the subframe as processed
        return -3;
      }
    }
    n->subframe_start_index = 0;  // Mark the subframe as processed
    n->next_subframe_id++;

    if (sf_id == 3) {
      // Got all of subframes 1 to 3
      n->next_subframe_id = 1;      // Make sure we start again next time

      // Now let's actually decode the ephemeris...
      decode_ephemeris(n->frame_words, e);

      return 1;

    }
  } else {  // didn't get the subframe that we want next
      n->next_subframe_id = 1;      // Make sure we start again next time
      n->subframe_start_index = 0;  // Mark the subframe as processed
  }

  return 0;

}

static u8 isNormal(u8 preamble_candidate)
{
  if (preamble_candidate == 0x53 ||
      preamble_candidate == 0x9A ||
      preamble_candidate == 0xC6)
    return 1;
  return 0;
}

static u8 isInverse(u8 preamble_candidate)
{
  if (preamble_candidate == 0xAC ||
      preamble_candidate == 0x65 ||
      preamble_candidate == 0x39)
    return 1;
  return 0;
}

s16 l1_sbas_process_subframe(l1_sbas_nav_msg_t *n, ephemeris_xyz_t *e,
                            sbas_almanac_t *alm)
{
  u16 i = 0;
  u8 *buffer = n->decoded;
  u16 total_bits = sizeof(n->decoded) * 8;
  bool found = false;

  /*
   *Loops through data buffer and finds first message with
   *good preamble and CRC. If the first preamble is an inverse one then the
   *full buffer is inversed (bit inverse).
   */
  for(i = 0; i < total_bits - 8 && !found; i++) {
    u8 preamble_candidate = getbitu(buffer, i, 8);
    if (isNormal(preamble_candidate)) {
      u32 msg_crc = getbitu(buffer, i + 226, 24);
      u8 crc_buf[29] = {0};
      crc_buf[0] = getbitu(buffer, i, 2);
      for(u8 k = 0; k < 28; k++) {
        crc_buf[k + 1] = getbitu(buffer, i + 2 + k * 8, 8);
      }
      u32 computed_crc = crc24q(crc_buf, 29, 0);
      if (computed_crc == msg_crc) {
        found = true;
      }
    } else if (isInverse(preamble_candidate)) {
      for (u8 j = 0; j < total_bits/8 + 1; j++)
        buffer[j] = ~buffer[j];
      u32 msg_crc = getbitu(buffer, i + 226, 24);
      u8 crc_buf[29] = {0};
      crc_buf[0] = getbitu(buffer, i, 2);
      for(u8 k = 0; k < 28; k++) {
        crc_buf[k + 1] = getbitu(buffer, i + 2 + k * 8, 8);
      }
      u32 computed_crc = crc24q(crc_buf, 29, 0);
      if (computed_crc == msg_crc) {
        found = true;
      } else {
        for (u8 j = 0; j < total_bits/8; j++)
          buffer[j] = ~buffer[j];
      }
    }
  }

  /*
   *Align with start of first good message from the buffer.
   */
  i--;
  s16 preamble_index = -1;
  for (u16 j = i; j < total_bits - 250; j += 250) {
    /*
     *Check CRC-24Q.
     */
    u32 crc = getbitu(buffer, j + 226, 24);
    u8 type = getbitu(buffer, j + 8, 6);
    u8 crc_buf[29] = {0};
    crc_buf[0] = getbitu(buffer, j, 2);
    for(u8 i = 0; i < 28; i++) {
      crc_buf[i + 1] = getbitu(buffer, j + 2 + i * 8, 8);
    }
    u32 computed_crc = crc24q(crc_buf, 29, 0);
    if (computed_crc != crc) {
      /*
       *log_info("Msg %d has failed CRC check.", type);
       */
      continue;
    }

    /*
     *Get first preamble which is aligned with 6 second GPS subframe.
     */
    u8 preamble_candidate = getbitu(buffer, j, 8);
    if (preamble_candidate == 0x53 || preamble_candidate == 0xAC) {
        /*
         *preamble_candidate == 0x9A || preamble_candidate == 0x65 ||
         *preamble_candidate == 0xC6 || preamble_candidate == 0x39) {
         */
      preamble_index = total_bits - j;
    }
    if (n->polarity == 0) {
      n->msg_normal++;
      n->msg_normal %= 50;
    } else {
      n->msg_inverse++;
      n->msg_inverse %= 50;
    }

    /*
     *Parse for GEO Ephemeris.
     */
    if (type == 9) {
      assert(e != NULL);
      u16 prev = j + 8 + 6;

      e->iod = getbitu(buffer, prev, 8); prev += 8;
      e->toa = getbitu(buffer, prev, 13) * 16; prev += 13;
      e->ura = getbitu(buffer, prev, 4); prev += 4;

      e->pos[0] = getbits(buffer, prev, 30) * 0.08; prev += 30;
      e->pos[1] = getbits(buffer, prev, 30) * 0.08; prev += 30;
      e->pos[2] = getbits(buffer, prev, 25) * 0.4; prev += 25;

      e->rate[0] = getbits(buffer, prev, 17) * 0.000625; prev += 17;
      e->rate[1] = getbits(buffer, prev, 17) * 0.000625; prev += 17;
      e->rate[2] = getbits(buffer, prev, 18) * 0.004; prev += 18;

      e->acc[0] = getbits(buffer, prev, 10) * 0.0000125; prev += 10;
      e->acc[1] = getbits(buffer, prev, 10) * 0.0000125; prev += 10;
      e->acc[2] = getbits(buffer, prev, 10) * 0.0000625; prev += 10;

      e->a_gf0 = getbits(buffer, prev, 12) * pow(2, -31); prev += 12;
      e->a_gf1 = getbits(buffer, prev, 8) * pow(2, -40); prev += 12;

      e->valid = 1;
      e->healthy = 1;

    } else if (type == 17) { /* Parse for GEO Almanac */
      assert(alm != NULL);
      u8 base = 8 + 6;
      u8 pass = 0;

      for (u16 i = 0; i < 67*3; i+= 67) {
        u16 prev = i + j + base;
        alm[pass].data_id = getbitu(buffer, prev, 2); prev += 2;
        alm[pass].sid.prn = getbitu(buffer, prev, 8) - 1; prev += 8;
        alm[pass].sid.constellation = SBAS_CONSTELLATION;
        alm[pass].sid.band = L1_BAND;
        alm[pass].health = getbitu(buffer, prev, 8); prev += 8;
        if (alm[pass].health == 0)
          alm[pass].valid = 1;

        alm[pass].x = getbits(buffer, prev, 15) * 2600; prev += 15;
        alm[pass].y = getbits(buffer, prev, 15) * 2600; prev += 15;
        alm[pass].z = getbits(buffer, prev, 9) * 26000; prev+= 9;

        alm[pass].x_rate = getbits(buffer, prev, 3) * 10; prev += 3;
        alm[pass].y_rate = getbits(buffer, prev, 3) * 10; prev += 3;
        alm[pass].z_rate = getbits(buffer, prev, 4) * 40.96; prev += 4;

        pass++;
      }
      alm[0].t0 = alm[1].t0 = alm[2].t0 = getbitu(buffer, j + base + 67*3, 11);
    }
  }

  return preamble_index;
}

