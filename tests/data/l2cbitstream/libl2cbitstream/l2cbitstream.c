/*
* Copyright (C) 2016 Swift Navigation Inc.
* Contact: Pasi Miettinen <pasi.miettinen@exafore.com>
*
* This source is subject to the license found in the file 'LICENSE' which must
* be be distributed together with this source. All other rights reserved.
*
* THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
* EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
*/


#include "l2cbitstream.h"
#include <libswiftnav/bits.h>
#include <libswiftnav/edc.h>

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BITS_IN_BYTE 8

#define CNAVMSG_OFFSET_BITS            4
#define CNAVMSG_PREAMBLE               0x8B
#define CNAVMSG_PREAMBLE_LEN_BITS      8
#define CNAVMSG_PRN_LEN_BITS           6
#define CNAVMSG_TYPE_LEN_BITS          6
#define CNAVMSG_CRC_LEN_BYTES          3
#define CNAVMSG_CRC_LEN_BITS           (BITS_IN_BYTE * CNAVMSG_CRC_LEN_BYTES)
#define CNAVMSG_CRC_IDX_BITS           276
#define CNAVMSG_TOW_COUNT_IDX_BITS     20
#define CNAVMSG_TOW_COUNT_INC          2 // equals 12 seconds
#define CNAVMSG_TOW_COUNT_LEN_BITS     17
#define CNAVMSG_TOW_COUNT_IGNORE_LSB   2

u8 get_l2c_message_length(void)
{
  return CNAVMSG_LEN_BYTES;
}

/** Creates L2C message
 *
 * Fills buffer with L2C message frame including preamble
 * and CRC-24q. Bits are right-aligned so four leftmost bits
 * of the 38 byte buffer are always zero.
 *
 * \param[out] au_message      Buffer for the message.
 * \param[in]  prn        GPS SV PRN number [0..31].
 * \param[in]  msg_id     Message idntifier [0..31].
 * \param[in]  tow        ToW in 6 second units [0..131071].
 * \return true if success
 *         false if error
 */
bool get_l2c_message(u8 *au_message, u8 prn, u8 msg_id, u32 tow)
{
  u32 q_idx = CNAVMSG_OFFSET_BITS;

  memset(au_message, 0, CNAVMSG_LEN_BYTES);

  /*
   * Offset is 4 -> 4 leftmost bits are always '0'.
   * Preamble.
   */
  setbitu(au_message, q_idx, CNAVMSG_PREAMBLE_LEN_BITS, CNAVMSG_PREAMBLE);
  q_idx += CNAVMSG_PREAMBLE_LEN_BITS;

  /* PRN */
  setbitu(au_message, q_idx, CNAVMSG_PRN_LEN_BITS, prn);
  q_idx += CNAVMSG_PRN_LEN_BITS;

  /* Message Type */
  setbitu(au_message, q_idx, CNAVMSG_TYPE_LEN_BITS, msg_id);
  q_idx += CNAVMSG_TYPE_LEN_BITS;

  /* 17 MSBs of the TOW count -> ignore 2 LSBs */
  setbitu(au_message, q_idx, CNAVMSG_TOW_COUNT_LEN_BITS, tow);
  q_idx += CNAVMSG_TOW_COUNT_LEN_BITS;

  /* Random payload */
  while (q_idx < CNAVMSG_CRC_IDX_BITS + CNAVMSG_OFFSET_BITS)
  {
    /* ignore sign bit (always 0) */
    setbitu(au_message, q_idx, 31, rand());
    q_idx += 31;
  }
  q_idx = CNAVMSG_CRC_IDX_BITS + CNAVMSG_OFFSET_BITS;

  setbitu(au_message, q_idx, CNAVMSG_CRC_LEN_BITS,
    crc24q(au_message, CNAVMSG_LEN_BYTES - CNAVMSG_CRC_LEN_BYTES, 0));
  q_idx += CNAVMSG_TOW_COUNT_LEN_BITS;

  return true;
}

/** Fills the buffer with 4bit offset removed
 *
 * \param au_msg_write    buffer to fill
 * 
 * \return -1 if error, otherwise the amount of fetched messages
 */
static s32 get_aligned_message(u8 *au_msg_write)
{
  u8 u_bytes = 0;
  static bool b_offset = true;
  static u8 au_msg_nxt[CNAVMSG_LEN_BYTES] = {0};
  static u8 u_idx = CNAVMSG_LEN_BYTES;
  u32 q_tow_count = 0;
  u64 t_msg_amount = 0;

  /* First call */
  if (u_idx >= CNAVMSG_LEN_BYTES) {
    if (!get_l2c_message(au_msg_nxt, rand(), rand(), q_tow_count)) {
        return -1;
    } else {
      q_tow_count += CNAVMSG_TOW_COUNT_INC;
      u_idx = 0;
    }
  }

  /* Fill 304 bits long buffer with 300 bit messages */
  while (u_bytes < CNAVMSG_LEN_BYTES) {
    if (b_offset) {
      au_msg_write[u_bytes] = au_msg_nxt[u_idx++] << CNAVMSG_OFFSET_BITS;
    } else {
      au_msg_write[u_bytes++] = au_msg_nxt[u_idx++];
    }

    if (u_idx >= CNAVMSG_LEN_BYTES) {
      /* get next 300 bit message since the buffer is empty*/
      if (!get_l2c_message(au_msg_nxt, rand(), rand(), q_tow_count)) {
        return -1;
      } else {
        q_tow_count += CNAVMSG_TOW_COUNT_INC;
        u_idx = 0;
        t_msg_amount++;
        if (b_offset) {
          /* handle offset from previous buffer */
          au_msg_write[u_bytes++] |= au_msg_nxt[u_idx++] &0x0F;
          b_offset = false;
        } else {
          b_offset = true;
          continue;
        }
      }
    }

    if (b_offset) {
        au_msg_write[u_bytes++] |= au_msg_nxt[u_idx] >> CNAVMSG_OFFSET_BITS;
    }
  }
  return t_msg_amount;
}

/** Fills file with L2C messages
 *
 * Writes L2C messages into given file. Utilizes
 * get_aligned_message function. Flushes, but doesn't close
 * the file.
 *
 * \param i_fileno                Filenumber from fileno()
 * \param t_wanted_msg_amount   Indicates how many messages
 *                                should be written. Depending
 *                                on the offset the final result
 *                                can be 1 message off.
 * \return -1 if error, otherwise the amount of fetched messages
 */
s32 write_l2c_to_file(int i_fileno, u64 t_wanted_msg_amount)
{
  u64 t_msg_amount = 0;
  u8 au_msg_write[CNAVMSG_LEN_BYTES] = {0};
  FILE *pz_file = fdopen(i_fileno, "wb");

  srand(time(NULL) * clock());

  if (NULL == pz_file) {
    printf("ERROR opening file\n");
    return -1;
  }

  t_wanted_msg_amount++;
  while (t_msg_amount < t_wanted_msg_amount) {
    s32 i_ret = get_aligned_message(au_msg_write);

    if (0 > i_ret) {
      printf("ERROR getting aligned message\n");
      fflush(pz_file);
      return t_msg_amount;
    }

    t_msg_amount += i_ret;

    if (fwrite(au_msg_write, sizeof(au_msg_write[0]),
      sizeof(au_msg_write), pz_file) != sizeof(au_msg_write)) {
      printf("ERROR writing file\n");
      fflush(pz_file);
      return t_msg_amount;
    }
    memset(au_msg_write, 0, sizeof(au_msg_write)); 
  }

  fflush(pz_file);

  return t_msg_amount;
}
