/*
 * Copyright (C) 2013 Swift Navigation Inc.
 * Contact: Fergus Noble <fergus@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#include "bits.h"

/** \defgroup bits Bit Utils
 * Bit field packing, unpacking and utility functions.
 * \{ */

/** Get bit field from buffer as an unsigned integer.
 * Unpacks `len` bits at bit position `pos` from the start of the buffer.
 * Maximum bit field length is 32 bits, i.e. `len <= 32`.
 *
 * \param pos Position in buffer of start of bit field in bits.
 * \param len Length of bit field in bits.
 * \return Bit field as an unsigned value.
 */
u32 getbitu(const u8 *buff, u32 pos, u8 len)
{
    u32 bits = 0;

    for (u32 i = pos; i < pos + len; i++) {
      bits = (bits << 1) +
             ((buff[i/8] >> (7 - i%8)) & 1u);
    }

    return bits;
}

/** Get bit field from buffer as a signed integer.
 * Unpacks `len` bits at bit position `pos` from the start of the buffer.
 * Maximum bit field length is 32 bits, i.e. `len <= 32`.
 *
 * This function sign extends the `len` bit field to a signed 32 bit integer.
 *
 * \param pos Position in buffer of start of bit field in bits.
 * \param len Length of bit field in bits.
 * \return Bit field as a signed value.
 */
s32 getbits(const u8 *buff, u32 pos, u8 len)
{
    s32 bits = (s32)getbitu(buff, pos, len);

    /* Sign extend, taken from:
     * http://graphics.stanford.edu/~seander/bithacks.html#VariableSignExtend
     */
    s32 m = 1u << (len - 1);
    return (bits ^ m) - m;
}

/** Set bit field in buffer from an unsigned integer.
 * Packs `len` bits into bit position `pos` from the start of the buffer.
 * Maximum bit field length is 32 bits, i.e. `len <= 32`.
 *
 * \param pos Position in buffer of start of bit field in bits.
 * \param len Length of bit field in bits.
 * \param data Unsigned integer to be packed into bit field.
 */
void setbitu(u8 *buff, u32 pos, u32 len, u32 data)
{
  u32 mask = 1u << (len - 1);

  if (len <= 0 || 32 < len)
    return;

  for (u32 i = pos; i < pos + len; i++, mask >>= 1) {
    if (data & mask)
      buff[i/8] |= 1u << (7 - i % 8);
    else
      buff[i/8] &= ~(1u << (7 - i % 8));
  }
}

/** Set bit field in buffer from a signed integer.
 * Packs `len` bits into bit position `pos` from the start of the buffer.
 * Maximum bit field length is 32 bits, i.e. `len <= 32`.
 *
 * \param pos Position in buffer of start of bit field in bits.
 * \param len Length of bit field in bits.
 * \param data Signed integer to be packed into bit field.
 */
void setbits(u8 *buff, u32 pos, u32 len, s32 data)
{
  setbitu(buff, pos, len, (u32)data);
}

/** \} */

