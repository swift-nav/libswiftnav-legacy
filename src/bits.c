/*
 * Copyright (C) 2013, 2016 Swift Navigation Inc.
 * Contact: Fergus Noble <fergus@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#include <libswiftnav/bits.h>

#include <limits.h>
#include <string.h>

/** \defgroup bits Bit Utils
 * Bit field packing, unpacking and utility functions.
 * \{ */

/** Computes the parity of a 32-bit word.
 *
 * References:
 *   -# https://graphics.stanford.edu/~seander/bithacks.html#ParityParallel
 *
 * \param x Word for which to compute parity
 * \return 1 if there are an odd number of bits set.
 *         0 if there are an even number of bits set.
 */
u8 parity(u32 x)
{
  x ^= x >> 16;
  x ^= x >> 8;
  x ^= x >> 4;
  x &= 0xF;
  return (0x6996 >> x) & 1;
}


/** Get bit field from buffer as an unsigned integer.
 * Unpacks `len` bits at bit position `pos` from the start of the buffer.
 * Maximum bit field length is 32 bits, i.e. `len <= 32`.
 *
 * \param buff
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
 * \param buff
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
 * \param buff
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
 * \param buff
 * \param pos Position in buffer of start of bit field in bits.
 * \param len Length of bit field in bits.
 * \param data Signed integer to be packed into bit field.
 */
void setbits(u8 *buff, u32 pos, u32 len, s32 data)
{
  setbitu(buff, pos, len, (u32)data);
}

/**
 * Shift MSB bit buffer contents to the left.
 * The method performs in-place shift operation.
 *
 * \param[in,out] buf   Pointer to buffer head.
 * \param[in]     bits  Number of bytes in the buffer.
 * \param[in]     shift Number of bits for left shift operation.
 *
 * \return None
 */
void bitshl(void *buf, u32 size, u32 shift)
{
  if (shift > size * CHAR_BIT) {
    /* Quick check: if the shift is larger, than the buffer, zero the data */
    memset(buf, 0, size);
    return;
  }

  unsigned char       *dst = buf;                /* Destination byte. */
  const unsigned char *src = dst + shift / CHAR_BIT; /* First source byte,
                                                      * possibly incomplete. */

  u32 copy_bits  = size * CHAR_BIT - shift; /* Number of bits to move */
  u32 byte_shift = copy_bits % CHAR_BIT;    /* Shift of data */
  u32 full_bytes = copy_bits / CHAR_BIT;    /* Number of bytes to move */

  if (0 == byte_shift) {
    /* When moving data in character boundaries, use built-in functions: move
     * data, and then zero the tail. */
    memmove(dst, src, full_bytes);
    memset(dst + full_bytes, 0, size - full_bytes);
  } else {
    /* Create an accumulator: it will hold a value of two consecutive bytes */
    u32      acc   = *src++;
    for (u32 i = 0; i < full_bytes; ++i) {
      acc = (acc << CHAR_BIT) | *src++;
      *dst++ = acc >> byte_shift;
    }
    *dst++ = acc << CHAR_BIT >> byte_shift;
    if (full_bytes + 1 < size) {
      memset(dst, 0, size - full_bytes - 1);
    }
  }
}

/**
 * Performs block bit copy operation.
 *
 * The function copies given number of bits from the source to destination.
 *
 * \param[in,out] dst       Pointer to destination buffer.
 * \param[in]     dst_index Destination bit index.
 * \param[in]     src       Source buffer.
 * \param[in]     src_index Source bit index.
 * \param[in]     count     Number of bits to copy.
 *
 * \return None
 *
 * \todo This function can be optimized for copying aligned data and using
 *       proper native type like long.
 */
void bitcopy(void *dst, u32 dst_index, const void *src, u32 src_index,
              u32 count)
{
  u32 limit1 = count / 32;
  u32 limit2 = count % 32;

  for (u32 idx = 0; idx < limit1; ++idx) {
    u32 tmp = getbitu(src, src_index, 32);
    setbitu(dst, dst_index, 32, tmp);
    src_index += 32;
    dst_index += 32;
  }
  if (0 != limit2) {
    u32 tmp = getbitu(src, src_index, limit2);
    setbitu(dst, dst_index, limit2, tmp);
  }
}

/** \} */
