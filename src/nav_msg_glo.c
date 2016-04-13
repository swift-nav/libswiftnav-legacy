/*
 * Copyright (C) 2016 Swift Navigation Inc.
 * Contact: Dmitry Tatarinov <dmitry.tatarinov@exafore.com>
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
#include <stdio.h>
#include <math.h>
#include <libswiftnav/nav_msg_glo.h>
#include <libswiftnav/time.h>

/* Word Ft (accuracy of measurements), refer to GLO ICD, Table 4.4 */
const float f_t[] = { 1.0f, 2.0f, 2.5f, 4.0f, 5.0f, 7.0f, 10.0f, 12.0f, 14.0f,
                      16.0f, 32.0f, 64.0f, 128.0f, 256.0f, 512.0f, -1.0f };

/* Word P1 (Time interval between adjacent values of tb, minutes), refer Table 4.3 */
const u8 p1[] = { 0, 30, 45, 60 }; /* min */

/** Initialize the necessary parts of the nav message state structure.
 * \param n Pointer to GLO nav message structure to be initialized
 */
void nav_msg_init_glo(nav_msg_glo_t *n)
{
  memset(n, 0, sizeof(nav_msg_glo_t));
  n->next_string_id = 1; /* start parsing from string 1 */
  n->state = SYNC_TM;
}

/** Extract a word of n_bits length (n_bits <= 32) at position bit_index into
 * the subframe. Refer to bit index to Table 4.6 and 4.11 in GLO ICD 5.1 (pg. 34)
 * \param n pointer to GLO nav message structure to be parsed
 * \param bit_index number of bit the extract process start with. Range [1..85]
 * \param n_bits how many bits should be extracted [1..32]
 * \return word extracted from navigation string
 */
u32 extract_word_glo(const nav_msg_glo_t *n, u16 bit_index, u8 n_bits)
{
  assert(n_bits <= 32 && n_bits > 0);
  assert(bit_index <= 85 && bit_index > 0);

  /* Extract a word of n_bits length (n_bits <= 32) at position bit_index into
   * the GLO string.*/
  bit_index--;
  u32 word = 0;
  u8 bix_hi = bit_index >> 5;
  u8 bix_lo = bit_index & 0x1F;
  if (bix_lo + n_bits <= 32) {
    word = n->string_bits[bix_hi] >> bix_lo;
    word &= (0xffffffff << (32 - n_bits)) >> (32 - n_bits);
  } else {
    u8 s = 32 - bix_lo;
    word = extract_word_glo(n, bit_index + 1, s)
         | extract_word_glo(n, bit_index + 1 + s, n_bits - s) << s;
  }

  return word;
}

/** Navigation message GLO decoding update.
 * Called once per nav bit interval (10 ms).
 * Performs the necessary steps to store the nav bits in buffer.
 *
 * \param n GLO Nav message decode state struct
 * \param bit_val State of the nav bit to process, 0 or 1
 *
 * \return  1 if Glo nav string ready for decoding,
 *         -1 otherwise.
 */
s8 nav_msg_update_glo(nav_msg_glo_t *n, bool bit_val)
{
  s8 ret = -1;

  switch (n->state) {
  case SYNC_TM: /* try to find time mark */
    /* put incoming bit at the tail of the buffer */
    n->string_bits[0] <<= 1; /* use one word of buffer for that purpose */
    n->string_bits[0] |= bit_val;
    /* collected bits match time mark? if not stay at this state */
    if (extract_word_glo(n, 1, GLO_TM_LEN) == GLO_TM) {
      /* time mark found, next time start collecting data bits */
      n->meander_bits_cnt = 0;
      n->manchester = 0;
      n->state = GET_DATA_BIT;
      n->string_bits[0] = 0;
    }
    break;
  case GET_DATA_BIT: /* collect data bits of string */
    n->meander_bits_cnt++;
    n->manchester <<= 1;
    n->manchester |= bit_val; /* store incoming bit */
    /* did we take 2 bits of line code?
     * if yes, remove meander and store bit in buffer,
     * if no, stay at the state */
    if (n->meander_bits_cnt == 2) {
      /* shift whole buffer by 1 bit left */
      for (u8 i = NAV_MSG_GLO_STRING_BITS_LEN - 1; i > 0; i--) {
        u32 tmp = (n->string_bits[i] << 1)
                | ((n->string_bits[i - 1] & (1 << 31)) >> 31);
        n->string_bits[i] = tmp;
      }
      n->string_bits[0] <<= 1;
      /* store bit after meander removal to buffer */
      n->string_bits[0] |= (n->manchester ^ 2) & 1;
      n->current_head_bit_index++;
      n->meander_bits_cnt = 0;
      n->manchester = 0;
      /* did we received all bits of a string?
       * if yes, notify user and start searching time mark again*/
      if (n->current_head_bit_index == 85) {
        n->current_head_bit_index = 0;
        n->state = SYNC_TM;
        ret = 1;
      }
    }
    break;
  default:
    nav_msg_init_glo(n); //TODO: probably not needed to initialize next_string_id
    break;
  }
  return ret;
}

/** The function decodes a GLO navigation string.
 *  Assume we receive signal from GLONASS-M
 * \param n Pointer to nav_msg_glo_t structure which contains GLO string
 * \param e Pointer to Ephemeris to store
 * \return 0 - decode not completed,
 *         1 -- decode completed, all ephemeris data stored
 *         <0 -- in case of an error.
 */
s8 process_string_glo(nav_msg_glo_t *n, ephemeris_t *e)
{
  /* Extract and check dummy bit from GLO string, bit 85 in GLO string */
  if (extract_word_glo(n, 85, 1) != 0) {
    log_error("GLO dummy bit is not 0.");
    return -1;
  }
  /* Extract string number */
  u32 m = extract_word_glo(n, 81, 4);
  u32 ret;
  u8 sign;
  /* is the string we need? */
  if (n->next_string_id == m) {
    switch (m) {
    case 1: /* string 1 */
      /* extract x */
      ret = extract_word_glo(n, 9, 26);
      sign = extract_word_glo(n, 9 + 26, 1);
      e->glo.pos[0] = sign ?
              -1.0 * ret * pow(2, -11) * 1000.0 : ret * pow(2, -11) * 1000.0;
      /* extract Vx */
      ret = extract_word_glo(n, 41, 23);
      sign = extract_word_glo(n, 41 + 23, 1);
      e->glo.vel[0] = sign ?
              -1.0 * ret * pow(2, -20) * 1000.0 : ret * pow(2, -20) * 1000.0;
      /* extract Ax */
      ret = extract_word_glo(n, 36, 4);
      sign = extract_word_glo(n, 36 + 4, 1);
      e->glo.acc[0] = sign ?
              -1.0 * ret * pow(2, -30) * 1000.0 : ret * pow(2, -30) * 1000.0;
      /* extract tk */
      n->hrs = (u8) extract_word_glo(n, 65, 5);
      n->min = (u8) extract_word_glo(n, 70, 6);
      n->sec = (u8) extract_word_glo(n, 76, 1);

      n->next_string_id = 2;
      break;
    case 2: /* string 2 */
      /* extract y */
      ret = extract_word_glo(n, 9, 26);
      sign = extract_word_glo(n, 9 + 26, 1);
      e->glo.pos[1] = sign ?
              -1.0 * ret * pow(2, -11) * 1000.0 : ret * pow(2, -11) * 1000.0;
      /* extract Vy */
      ret = extract_word_glo(n, 41, 23);
      sign = extract_word_glo(n, 41 + 23, 1);
      e->glo.vel[1] = sign ?
              -1.0 * ret * pow(2, -20) * 1000.0 : ret * pow(2, -20) * 1000.0;
      /* extract Ay */
      ret = extract_word_glo(n, 36, 4);
      sign = extract_word_glo(n, 36 + 4, 1);
      e->glo.acc[1] = sign ?
              -1.0 * ret * pow(2, -30) * 1000.0 : ret * pow(2, -30) * 1000.0;
      /* extract MSB of B (if the bit is clear the SV is OK ) */
      e->healthy = extract_word_glo(n, 80, 1);
      /* extract P1 */
      e->fit_interval = p1[extract_word_glo(n, 77, 2)];

      n->next_string_id = 3;
      break;
    case 3: /* string 3 */
      /* extract z */
      ret = extract_word_glo(n, 9, 26);
      sign = extract_word_glo(n, 9 + 26, 1);
      e->glo.pos[2] = sign ?
              -1.0 * ret * pow(2, -11) * 1000.0 : ret * pow(2, -11) * 1000.0;
      /* extract Vz */
      ret = extract_word_glo(n, 41, 23);
      sign = extract_word_glo(n, 41 + 23, 1);
      e->glo.vel[2] = sign ?
              -1.0 * ret * pow(2, -20) * 1000.0 : ret * pow(2, -20) * 1000.0;
      /* extract Az */
      ret = extract_word_glo(n, 36, 4);
      sign = extract_word_glo(n, 36 + 4, 1);
      e->glo.acc[2] = sign ?
              -1.0 * ret * pow(2, -30) * 1000.0 : ret * pow(2, -30) * 1000.0;
      /* extract gamma */
      ret = extract_word_glo(n, 69, 10);
      sign = extract_word_glo(n, 69 + 10, 1);
      e->glo.gamma = sign ? -1.0 * ret * pow(2, -40) : ret * pow(2, -40);
      /* extract l, if the it is clear the SV is OK, so OR it with B */
      e->healthy |= extract_word_glo(n, 65, 1);

      n->next_string_id = 4;
      break;
    case 4: /* string 4 */
      /* extract tau */
      ret = extract_word_glo(n, 59, 21);
      sign = extract_word_glo(n, 59 + 21, 1);
      e->glo.tau = sign ? -1.0 * ret * pow(2, -30) : ret * pow(2, -30);
      /* extract n */
      e->sid.sat = extract_word_glo(n, 11, 5);
      /* extract Ft (URA) */
      e->ura = f_t[extract_word_glo(n, 30, 4)];
      /*extract Nt*/
      n->nt = (u16) extract_word_glo(n, 16, 11);

      n->next_string_id = 5;
      break;
    case 5: /* string 5 */
      /* extract N4 */
      n->n4 = (u8) extract_word_glo(n, 32, 5);

      n->next_string_id = 0; /* all string parsed */
      break;
    default:
      break;
    }
  }
  /* all needed strings decoded?
   * fill ephemeris structure if we was not able to do it before*/
  if (n->next_string_id == 0) {
    e->sid.code = CODE_GLO_L1CA;
    /* convert GLO TOE to GPS TOE */
    e->toe = glo_time2gps_time(n->nt, n->n4, n->hrs, n->min, n->sec);
    e->healthy ^= 1; /* invert healthy bit */
    e->valid = e->healthy; //NOTE: probably Valid needs to be defined by other way
    n->next_string_id = 1; /* start from string 1 */
    return 1;
  }


  return 0;
}
