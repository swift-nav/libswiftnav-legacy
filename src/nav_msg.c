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
#include <math.h>

#include "logging.h"
#include "constants.h"
#include "bits.h"
#include "nav_msg.h"

#define NAV_MSG_BIT_PHASE_THRES 10

void nav_msg_init(nav_msg_t *n)
{
  /* Initialize the necessary parts of the nav message state structure. */
  n->subframe_bit_index = 0;
  n->bit_phase = 0;
  n->bit_phase_ref = 0;
  n->bit_phase_count = 0;
  n->nav_bit_integrate = 0;
  n->subframe_start_index = 0;
  memset(n->subframe_bits, 0, sizeof(n->subframe_bits));
  n->next_subframe_id = 1;
}

static u32 extract_word(nav_msg_t *n, u16 bit_index, u8 n_bits, u8 invert)
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
  if (bit_index > NAV_MSG_SUBFRAME_BITS_LEN*32)
    bit_index -= NAV_MSG_SUBFRAME_BITS_LEN*32;

  u8 bix_hi = bit_index >> 5;
  u8 bix_lo = bit_index & 0x1F;
  u32 word = n->subframe_bits[bix_hi] << bix_lo;

  if (bix_lo) {
    bix_hi++;
    if (bix_hi == NAV_MSG_SUBFRAME_BITS_LEN)
      bix_hi = 0;
    word |=  n->subframe_bits[bix_hi] >> (32 - bix_lo);
  }

  if (invert)
    word = ~word;

  return word >> (32 - n_bits);
}

/* Check for a sign change in the correlation and add to the histogram.
 * After NAV_MSG_BIT_PHASE_THRES sign changes, the histogram is evaluated
 * and bit_phase_ref is updated.
 */
static void update_bit_sync(nav_msg_t *n, s32 corr_prompt_real, u8 ms)
{
  float dot = corr_prompt_real * n->prev_corr;
  n->prev_corr = corr_prompt_real;

  if (dot > 0)
    return;

  /* Sign change - Add to histogram */
  n->hist[n->bit_phase] += -dot;

  n->bit_phase_count++;
  if (n->bit_phase_count < NAV_MSG_BIT_PHASE_THRES)
    return;

  /* Find the max and run with it */
  float max = 0;
  u8 max_index = 0;
  for (u8 i = 0; i < 20; i++) {
    if (n->hist[i] > max) {
      max = n->hist[i];
      max_index = i;
    }
    n->hist[i] = 0;
  }

  n->bit_phase_ref = (max_index + 20 - ms) % 20;

  n->bit_phase_count = 0;
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

  /* Do we have bit phase lock yet? (Do we know which of the 20 possible PRN
   * offsets corresponds to the nav bit edges?) */
  n->bit_phase += ms;
  n->bit_phase %= 20;

  update_bit_sync(n, corr_prompt_real, ms);

  /* We have bit phase lock. */
  n->nav_bit_integrate += corr_prompt_real;

  if (n->bit_phase != n->bit_phase_ref) {
    return TOW_INVALID;
  }

  /* Dump the nav bit, i.e. determine the sign of the correlation over the
   * nav bit period. */

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
    n->nav_bit_integrate = 0;
    return -2;
  }

  /* Is bit 1? */
  if (n->nav_bit_integrate > 0) {
    n->subframe_bits[n->subframe_bit_index >> 5] |= \
      1 << (31 - (n->subframe_bit_index & 0x1F));
  } else {
    /* Integrated correlation is negative, so bit is 0. */
    n->subframe_bits[n->subframe_bit_index >> 5] &= \
      ~(1 << (31 - (n->subframe_bit_index & 0x1F)));
  }

  /* Zero the integrator for the next nav bit. */
  n->nav_bit_integrate = 0;

  n->subframe_bit_index++;
  if (n->subframe_bit_index == NAV_MSG_SUBFRAME_BITS_LEN*32)
    n->subframe_bit_index = 0;

  /* Yo dawg, are we still looking for the preamble? */
  if (!n->subframe_start_index) {
    /* We're going to look for the preamble at a time 360 nav bits ago,
     * then again 60 nav bits ago. */
    #define SUBFRAME_START_BUFFER_OFFSET (NAV_MSG_SUBFRAME_BITS_LEN*32 - 360)

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

        unsigned int TOW_trunc = extract_word(n,30,17,extract_word(n,29,1,0)); // bit 29 is D30* for the second word, where the TOW resides.
        TOW_trunc++;  // Increment it, to see what we expect at the start of the next subframe
        if (TOW_trunc >= 7*24*60*10)  // Handle end of week rollover
          TOW_trunc = 0;

        if (TOW_trunc == extract_word(n,330,17,extract_word(n,329,1,0))) {
          // We got two appropriately spaced preambles, and two matching TOW counts.  Pretty certain now.

          // The TOW in the message is for the start of the NEXT subframe.
          // That is, 240 nav bits' time from now, since we are 60 nav bits into the second subframe that we recorded.
          if (TOW_trunc)
            TOW_ms = TOW_trunc * 6000 - (300-60)*20;
          else  // end of week special case
            TOW_ms = 7*24*60*60*1000 - (300-60)*20;

        } else
          n->subframe_start_index = 0;  // the TOW counts didn't match - disregard.
      } else
        n->subframe_start_index = 0;    // didn't find a second preamble in the right spot - disregard.
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

bool subframe_ready(nav_msg_t *n) {
  return (n->subframe_start_index != 0);
}

s8 process_subframe(nav_msg_t *n, ephemeris_t *e) {
  // Check parity and parse out the ephemeris from the most recently received subframe

  // First things first - check the parity, and invert bits if necessary.
  // process the data, skipping the first word, TLM, and starting with HOW
  /* Complain if buffer overrun */
  if (n->overrun) {
    log_warn("nav_msg subframe buffer overrun!\n");
    n->overrun = false;
  }

  /* TODO: Check if inverted has changed and detect half cycle slip. */
  if (n->inverted != (n->subframe_start_index < 0)) {
    log_info("Nav phase flip\n");
  }
  n->inverted = (n->subframe_start_index < 0);

  if (!e) {
    log_error("process_subframe: CALLED WITH e = NULL!\n");
    n->subframe_start_index = 0;  // Mark the subframe as processed
    n->next_subframe_id = 1;      // Make sure we start again next time
    return -1;
  }
  u32 sf_word2 = extract_word(n, 28, 32, 0);
  if (nav_parity(&sf_word2)) {
      log_info("subframe parity mismatch (word 2)\n");
      n->subframe_start_index = 0;  // Mark the subframe as processed
      n->next_subframe_id = 1;      // Make sure we start again next time
      return -2;
  }

  u8 sf_id = sf_word2 >> 8 & 0x07;    // Which of 5 possible subframes is it?

  if (sf_id <= 3 && sf_id == n->next_subframe_id) {  // Is it the one that we want next?

    for (int w = 0; w < 8; w++) {   // For words 3..10
      n->frame_words[sf_id-1][w] = extract_word(n, 30*(w+2) - 2, 32, 0);    // Get the bits
      // MSBs are D29* and D30*.  LSBs are D1...D30
      if (nav_parity(&n->frame_words[sf_id-1][w])) {  // Check parity and invert bits if D30*
        log_info("subframe parity mismatch (word %d)\n", w+3);
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

