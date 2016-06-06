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

#ifndef LIBSWIFTNAV_CNAV_MSG_H
#define LIBSWIFTNAV_CNAV_MSG_H

#include <libswiftnav/common.h>
#include <libfec/fec.h>

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <limits.h>

/** \addtogroup GPS_L2
 * \{ */
/** \addtogroup gps_cnav_decoder
 * \{ */

/** Size of the Viterbi decoder history. */
#define GPS_L2_V27_HISTORY_LENGTH_BITS 64
/** Bits to accumulate before decoding starts. */
#define GPS_L2C_V27_INIT_BITS    (32)
/** Bits to decode at a time. */
#define GPS_L2C_V27_DECODE_BITS  (32)
/** Bits in decoder tail. We ignore them. */
#define GPS_L2C_V27_DELAY_BITS   (32)
/** L2C convolutional encoder constrain length */
#define GPS_L2C_V27_CONSTRAINT_LENGTH 7
/**
 * GPS CNAV message container.
 *
 * @sa cnav_msg_decoder_add_symbol
 */
typedef struct
{
  u8   prn;    /**< SV PRN. 0..31 */
  u8   msg_id; /**< Message id. 0..31 */
  u32  tow;    /**< GPS ToW in 6-second units. Multiply to 6 to get seconds. */
  bool alert;  /**< CNAV message alert flag */
  s8 bit_polarity; /**< Polarity of data bits */
} cnav_msg_t;

/**
 * GPS CNAV decoder component.
 * This component controls symbol decoding string.
 *
 * @sa cnav_msg_decoder_t
 */
typedef struct {
  v27_t          dec;           /**< Viterbi block decoder object */
  v27_decision_t decisions[GPS_L2_V27_HISTORY_LENGTH_BITS];
                                /**< Decision graph */
  unsigned char  symbols[(GPS_L2C_V27_INIT_BITS + GPS_L2C_V27_DECODE_BITS) * 2];
                                /**< Symbol buffer */
  size_t         n_symbols;     /**< Count of symbols in the symbol buffer */
  unsigned char  decoded[GPS_L2C_V27_DECODE_BITS + GPS_L2C_V27_DELAY_BITS];
                                /**< Decode buffer */
  size_t         n_decoded;     /**< Number of bits in the decode buffer */
  bool           preamble_seen; /**< When true, the decode buffer is aligned on
                                 *   preamble. */
  bool           invert;        /**< When true, indicates the bits are inverted */
  bool           message_lock;  /**< When true, indicates the message boundary
                                 *   is found. */
  bool           crc_ok;        /**< Flag that the last message had good CRC */
  size_t         n_crc_fail;    /**< Counter for CRC failures */
  bool           init;          /**< Initial state flag. When true, initial bits
                                 *   do not produce output. */
} cnav_v27_part_t;

/**
 * GPS CNAV message lock and decoder object.
 *
 * Decoder uses two Viterbi decoder objects to ensure the lock is acquired when
 * the input symbol phase is not known.
 */
typedef struct
{
  cnav_v27_part_t part1; /**< Decoder for odd symbol pairs */
  cnav_v27_part_t part2; /**< Decoder for even symbol pairs */
} cnav_msg_decoder_t;

const v27_poly_t *cnav_msg_decoder_get_poly(void);
void cnav_msg_decoder_init(cnav_msg_decoder_t *dec);
bool cnav_msg_decoder_add_symbol(cnav_msg_decoder_t *dec,
                                 unsigned char symbol,
                                 cnav_msg_t *msg,
                                 u32 *delay);

/** \} */
/** \} */

#endif /* LIBSWIFTNAV_CNAV_MSG_H */
