/*
 * Copyright (C) 2016 Swift Navigation Inc.
 * Contact: Valeri Atamaniouk <valeri@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#include <libswiftnav/edc.h>
#include <libswiftnav/bits.h>
#include <libswiftnav/cnav_msg.h>

#include <limits.h>
#include <string.h>

/** \defgroup GPS_L2 GPS L2
 * GPS L2 operations
 * \{ */
/** \defgroup gps_cnav_decoder Decoder
 * GPS L2C CNAV message decoding.
 * \{ */
/*
 * Block Viterbi decoding parameters.
 */
/** Viterbi decoder inverted polynome A */
#define GPS_L2C_V27_POLY_A       (0x4F) /* 0b01001111 - inverted 0171*/
/** Viterbi decoder inverted polynome B */
#define GPS_L2C_V27_POLY_B       (0x6D) /* 0b01101101 - inverted 0133 */

/*
 * GPS L2C message constants.
 */

/** GPS L2C preamble */
#define GPS_CNAV_PREAMBLE1          (0b10001011u)
/** Inverted GPS L2C preamble */
#define GPS_CNAV_PREAMBLE2          (0b01110100u)
/** GPS L2C preamble length in bits */
#define GPS_CNAV_PREAMBLE_LENGTH    (8)
/** GPS L2C CNAV message length in bits */
#define GPS_CNAV_MSG_LENGTH         (300)
/** GPS LC2 CNAV CRC length in bits */
#define GPS_CNAV_MSG_CRC_LENGTH     (24)
/** GPS L2C CNAV message payload length in bits */
#define GPS_CNAV_MSG_DATA_LENGTH    (GPS_CNAV_MSG_LENGTH-GPS_CNAV_MSG_CRC_LENGTH)
/** GPS L2C CNAV message lock detector threshold */
#define GPS_CNAV_LOCK_MAX_CRC_FAILS (10)

// forward declarations for test access
void _cnav_bitlshift(void *buf, size_t bits, u32 shift);
u32 _cnav_compute_crc(cnav_v27_part_t *part);
u32 _cnav_extract_crc(const cnav_v27_part_t *part);
void _cnav_rescan_preamble(cnav_v27_part_t *part);
void _cnav_add_symbol(cnav_v27_part_t *part, u8 ch);
void _cnav_msg_invert(cnav_v27_part_t *part);
bool _cnav_msg_decode(cnav_v27_part_t *part, cnav_msg_t *msg, u32 *delay);

/**
 * Shift MSB bit buffer contents to the left.
 * The method performs in-place shift operation.
 *
 * \param[in] buf   Pointer to buffer head.
 * \param[in] bits  Number of bytes in the buffer.
 * \param[in] shift Number of bits for left shift operation.
 *
 * \return None
 *
 * @\private
 */
void _cnav_bitlshift(void *buf, size_t size, u32 shift)
{
  if (shift > size * CHAR_BIT) {
    /* Quick check: if the shift is larger, than the buffer, zero the data */
    memset(buf, size, 0);
    return;
  }

  unsigned char       *dst = buf;                // Destination byte.
  const unsigned char *src = dst + shift / CHAR_BIT; // First source byte,
                                                 // possibly incomplete.

  size_t   copy_bits  = size * CHAR_BIT - shift; // Number of bits to move
  u32      byte_shift = copy_bits % CHAR_BIT;    // Shift of data
  size_t   full_bytes = copy_bits / CHAR_BIT;    // Number of bytes to move

  if (0 == byte_shift) {
    /* When moving data in character boundaries, use built-in functions: move
     * data, and then zero the tail. */
    memmove(dst, src, full_bytes);
    memset(dst + full_bytes, 0, size - full_bytes);
  } else {
    /* Create an accumulator: it will hold a value of two consecutive bytes */
    u32      acc   = *src++;
    for (size_t i = 0; i < full_bytes; ++i) {
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
 * Computes CRC-24Q from a CNAV message buffer.
 * CRC-24Q is computed for 274 bits. For a purpose of 8-bit alignment, the
 * message is assumed to be prepended with four zero bits.
 *
 * \param[in] part Decoder component with payload
 *
 * \return CRC-24Q value.
 *
 * \private
 */
u32 _cnav_compute_crc(cnav_v27_part_t *part)
{
  u32 crc = crc24q_bits(0, part->decoded, GPS_CNAV_MSG_DATA_LENGTH,
                        part->invert);

  return crc;
}

/**
 * Extracts CRC-24Q from a CNAV message buffer.
 * CRC-24Q value is the last 24 bits from 300 bits message buffer.
 *
 * \param[in] part Decoded component with payload.
 *
 * \return CRC24-Q value.
 *
 * \private
 */
u32 _cnav_extract_crc(const cnav_v27_part_t *part)
{
  u32 crc = getbitu(part->decoded, 300-24, 24);
  if (part->invert) {
    crc ^= 0xFFFFfF;
  }
  return crc;
}

/**
 * Helper to rescan for preamble in the received buffer.
 * Occasionally there could be a false lock on message contents if it the
 * preamble sequence is encountered in the message body. For this case, the
 * function performs for a scan for a preamble with a different offset:
 * - When found, the preamble octet is moved into the head of the buffer.
 * - When not found, only 7 bits are left in the buffer.
 *
 * \param[in,out] part Decoded component.
 *
 * \return None
 *
 * \private
 */
void _cnav_rescan_preamble(cnav_v27_part_t *part)
{
  part->preamble_seen = false;

  if (part->n_decoded > GPS_CNAV_PREAMBLE_LENGTH + 1) {
    for (size_t i = 1, j = part->n_decoded - GPS_CNAV_PREAMBLE_LENGTH;
         i < j; ++i) {
      u32 c = getbitu(part->decoded, i, GPS_CNAV_PREAMBLE_LENGTH);
      if (GPS_CNAV_PREAMBLE1 == c || GPS_CNAV_PREAMBLE2 == c) {
        part->preamble_seen = true;
        part->invert = (GPS_CNAV_PREAMBLE2 == c);
        /* We shift the accumulated bits to the beginning of the buffer */
        _cnav_bitlshift(part->decoded,
                        sizeof(part->decoded),
                        i);
        part->n_decoded -= i;
        break;
      }
    }
  }
  if (!part->preamble_seen && part->n_decoded >= GPS_CNAV_PREAMBLE_LENGTH) {
    _cnav_bitlshift(part->decoded,
                    sizeof(part->decoded),
                    part->n_decoded - GPS_CNAV_PREAMBLE_LENGTH + 1);
    part->n_decoded = GPS_CNAV_PREAMBLE_LENGTH - 1;
  }
}

/**
 * Feed a symbol into Viterbi decoder instance.
 *
 * The method uses block Viterbi decoder. It first accumulates initial number of
 * symbols, and after that runs decoding every time the buffer is full. Only
 * some of the decoded symbols are used.
 *
 * \param[in,out] part Decoder object
 * \param[in]     s    Symbol (0x00 - Hard 0, 0xFF - Hard 1)
 *
 * \return None
 *
 * \private
 */
void _cnav_add_symbol(cnav_v27_part_t *part, u8 ch)
{
  part->symbols[part->n_symbols++] = ch;

  if (part->init) {
    /* Initial step - load more symbols without decoding. */
    if (part->n_symbols < (GPS_L2C_V27_INIT_BITS + GPS_L2C_V27_DECODE_BITS) * 2) {
      return;
    }
    part->init = false;
  }
  else if (part->n_symbols < GPS_L2C_V27_DECODE_BITS * 2) {
    /* Wait until decoding block is accumulated */
    return;
  }

  /* Feed accumulated symbols into the buffer, reset the number of accumulated
   * symbols. */
  v27_update(&part->dec, part->symbols, part->n_symbols / 2);
  part->n_symbols = 0;

  /* Decode N+M bits, where:
   * - N - Number of bits to put into decoded buffer
   * - M - Number of bits in the tail to ignore.
   */
  unsigned char tmp_bits[ (GPS_L2C_V27_DECODE_BITS + GPS_L2C_V27_DELAY_BITS +
                           CHAR_BIT - 1) / CHAR_BIT];

  v27_chainback_likely(&part->dec, tmp_bits,
                       GPS_L2C_V27_DECODE_BITS + GPS_L2C_V27_DELAY_BITS);

  /* Read decoded bits and add them to the decoded buffer */
  u32 to_add = getbitu(tmp_bits, 0, GPS_L2C_V27_DECODE_BITS);
  setbitu(part->decoded, part->n_decoded, GPS_L2C_V27_DECODE_BITS, to_add);
  part->n_decoded += GPS_L2C_V27_DECODE_BITS;

  /* Depending on the decoder state, one of the following actions are
   * possible:
   * - If no message lock
   *   - If no preamble seen - look for preamble
   *   - If preamble seen - collect 300 bits
   *     - If 300 bits are collected - verify CRC
   *       - If CRC is OK - message lock is acquired
   *       - If CRC fails - rescan for preamble
   *         - If found - continue collecting 300 bits
   *         - If not found - continue preamble wait
   * - If message lock
   *   - If 300 bits collected, compute CRC
   *     - If CRC is OK, message can be decoded
   *     - If CRC is not OK, discard data
   */

  bool retry = true;
  while (retry) {
    retry = false;

    if (!part->preamble_seen) {
      /* Rescan for preamble if possible. The first bit is ignored. */
      _cnav_rescan_preamble(part);
    }
    if (part->preamble_seen && GPS_CNAV_MSG_LENGTH <= part->n_decoded) {

      /* We have collected 300 bits starting from message preamble. Now try
       * to compute CRC-24Q */
      u32 crc  = _cnav_compute_crc(part);
      u32 crc2 = _cnav_extract_crc(part);

      if (part->message_lock) {
        /* We have message lock */
        part->crc_ok = (crc == crc2);
        if (part->crc_ok) {
          /* Reset message lock counter */
          part->n_crc_fail = 0;
        } else {
          /* Increment message lock counter */
          part->n_crc_fail++;
          if (part->n_crc_fail > GPS_CNAV_LOCK_MAX_CRC_FAILS) {
            /* CRC has failed too many times - drop the lock. */
            part->n_crc_fail    = 0;
            part->message_lock  = false;
            part->preamble_seen = false;
            /* Try to find a new preamble, reuse data from buffer. */
            retry = true;
          }
        }
      } else if (crc == crc2) {
        /* CRC match - message can be decoded */
        part->message_lock = true;
        part->crc_ok       = true;
        part->n_crc_fail   = 0;
      } else {
        /* There is no message lock and the CRC check fails. Assume there is
         * false positive lock - rescan for preamble. */
        part->crc_ok = false;
        part->preamble_seen = false;

        /* CRC mismatch - try to re-scan for preamble */
        retry = true;
      }
    }
    else
    {
      /* No preamble or preamble and less than 300 bits decoded */
    }
  }
}

/**
 * Invert message bits in the buffer.
 *
 * The method inverts 300 bits of the message data. The data inversion flag is
 * also updated.
 *
 * \param[in,out] part Decoder component with a message buffer.
 *
 * \return None
 */
void _cnav_msg_invert(cnav_v27_part_t *part)
{
  for (size_t i = 0; i < sizeof(part->decoded); i++) {
    part->decoded[i] ^= 0xFFu;
  }
}

/**
 * Performs CNAV message decoding.
 *
 * This function decoded CNAV message, if the following conditions are met:
 * - Buffer contains 300 bits.
 * - First 8 bits are matching direct or inverse preamble.
 * - Message data CRC matches one in the buffer.
 *
 * In case the message starts with inverted preamble, the data is inverted
 * before parsing.
 *
 * \param[in,out] part Decoder component.
 * \param[out]    msg  Container for a decoded message.
 * \param[out]    delay Delay of the message in symbols.
 *
 * \retval true The message has been decoded, and \a msg container is populated.
 * \retval false Not enough data or CRC is not correct.
 *
 * \private
 */
bool _cnav_msg_decode(cnav_v27_part_t *part, cnav_msg_t *msg, u32 *delay)
{
  bool res = false;
  if (GPS_CNAV_MSG_LENGTH <= part->n_decoded) {
    if (part->crc_ok) {
      /* CRC is OK */

      if (part->invert) {
        _cnav_msg_invert(part);
      }

      msg->prn    = getbitu(part->decoded, 8, 6);
      msg->msg_id = getbitu(part->decoded, 14, 6);
      msg->tow    = getbitu(part->decoded, 20, 17);
      msg->alert  = getbitu(part->decoded, 37, 1) ? true : false;

      *delay = (part->n_decoded - GPS_CNAV_MSG_LENGTH + GPS_L2C_V27_DELAY_BITS)
               * 2 + part->n_symbols;

      if (part->invert) {
        _cnav_msg_invert(part);
      }
      res = true;
    } else {
      /* CRC mismatch - no decoding */
    }
    _cnav_bitlshift(part->decoded, sizeof(part->decoded), GPS_CNAV_MSG_LENGTH);
    part->n_decoded -= GPS_CNAV_MSG_LENGTH;
  }

  return res;
}

/**
 * Initialize CNAV decoder.
 *
 * CNAV decoder contains two Viterbi decoders that are used to estimate bit and
 * message boundary.
 *
 * \param[out] dec Decoder structure.
 *
 * \return None
 */
void cnav_msg_decoder_init(cnav_msg_decoder_t *dec)
{
  static const signed char polynomial[2] = {
    GPS_L2C_V27_POLY_A, GPS_L2C_V27_POLY_B };

  memset(dec, 0, sizeof(*dec));
  v27_poly_init(&dec->poly, polynomial);
  v27_init(&dec->part1.dec,
           dec->part1.decisions,
           GPS_L2_V27_HISTORY_LENGTH_BITS,
           &dec->poly,
           0);
  v27_init(&dec->part2.dec,
           dec->part2.decisions,
           GPS_L2_V27_HISTORY_LENGTH_BITS,
           &dec->poly,
           0);
  dec->part1.init = true;
  dec->part2.init = true;
  _cnav_add_symbol(&dec->part2, 0x80);
}

/**
 * Adds a received symbol to decoder.
 *
 * \param[in,out] dec    Decoder object.
 * \param[in]     symbol Symbol value probability, where 0x00 - 100% of 0,
 *                       0xFF - 100% of 1.
 * \param[out]    msg    Buffer for decoded message. The message is available
 *                       only when message lock is acquired and CRC is correct.
 * \param[out]    pdelay Delay of message generation in symbols.
 *
 * \retval true  The message has been decoded. ToW parameter is available.
 * \retval false More data is required.
 */
bool cnav_msg_decoder_add_symbol(cnav_msg_decoder_t *dec,
                                 u8                  symbol,
                                 cnav_msg_t         *msg,
                                 u32                *pdelay)
{
  _cnav_add_symbol(&dec->part1, symbol);
  _cnav_add_symbol(&dec->part2, symbol);

  if (dec->part1.message_lock) {
    /* Flush data in decoder. */
    dec->part2.n_decoded = 0;
    dec->part2.n_symbols = 0;
    return _cnav_msg_decode(&dec->part1, msg, pdelay);
  }
  if (dec->part2.message_lock) {
    /* Flush data in decoder. */
    dec->part1.n_decoded = 0;
    dec->part1.n_symbols = 0;
    return _cnav_msg_decode(&dec->part2, msg, pdelay);
  }

  return false;
}

/** \} */
/** \} */
