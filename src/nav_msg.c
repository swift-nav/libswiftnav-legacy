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

#include <libswiftnav/logging.h>
#include <libswiftnav/constants.h>
#include <libswiftnav/bits.h>
#include <libswiftnav/nav_msg.h>

void nav_msg_init(nav_msg_t *n)
{
  /* Initialize the necessary parts of the nav message state structure. */
  memset(n, 0, sizeof(nav_msg_t));
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

/** Navigation message decoding update.
 * Called once per nav bit interval. Performs the necessary steps to
 * store the nav bits and decode them.
 *
 * Also extracts and returns the GPS time of week each time a new subframe is
 * received.
 *
 * \param n Nav message decode state struct
 * \param bit_val State of the nav bit to process
 *
 * \return The GPS time of week in milliseconds of the current code phase
 *         rollover, or `TOW_INVALID` (-1) if unknown
 */
s32 nav_msg_update(nav_msg_t *n, bool bit_val)
{
  s32 TOW_ms = TOW_INVALID;

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

  if (n->subframe_bit_index >= NAV_MSG_SUBFRAME_BITS_LEN*32) {
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
  char buf[SID_STR_LEN_MAX];
  sid_to_string(buf, sizeof(buf), e->sid);

  // Check parity and parse out the ephemeris from the most recently received subframe

  if (!e) {
    log_error("process_subframe: CALLED WITH e = NULL!");
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
    log_error("%s Nav phase flip - half cycle slip detected, "
              "but not corrected", buf);
    /* TODO: declare phase ambiguity to IAR */
  }

  /* Complain if buffer overrun */
  if (n->overrun) {
    log_error("%s nav_msg subframe buffer overrun!", buf);
    n->overrun = false;
  }

  /* Extract word 2, and the last two parity bits of word 1 */
  u32 sf_word2 = extract_word(n, 28, 32, 0);
  if (nav_parity(&sf_word2)) {
    log_info("%s subframe parity mismatch (word 2)", buf);
    n->subframe_start_index = 0;  // Mark the subframe as processed
    n->next_subframe_id = 1;      // Make sure we start again next time
    return -2;
  }

  n->alert = sf_word2 >> 12 & 0x01; // Alert flag, bit 18
  if (n->alert) {
    char buf[SID_STR_LEN_MAX];
    sid_to_string(buf, sizeof(buf), e->sid);
    log_warn("%s alert flag set! Ignoring satellite.", buf);
  }

  u8 sf_id = sf_word2 >> 8 & 0x07;    // Which of 5 possible subframes is it? bits 20-22

  if (sf_id <= 3 && sf_id == n->next_subframe_id) {  // Is it the one that we want next?

    for (int w = 0; w < 8; w++) {   // For words 3..10
      n->frame_words[sf_id-1][w] = extract_word(n, 30*(w+2) - 2, 32, 0);    // Get the bits
      // MSBs are D29* and D30*.  LSBs are D1...D30
      if (nav_parity(&n->frame_words[sf_id-1][w])) {  // Check parity and invert bits if D30*
        log_info("%s subframe parity mismatch (word %d)", buf, w+3);
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
