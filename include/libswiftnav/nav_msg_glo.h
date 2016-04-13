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

#ifndef LIBSWIFTNAV_NAV_MSG_GLO_H
#define LIBSWIFTNAV_NAV_MSG_GLO_H

#include <libswiftnav/common.h>
#include <libswiftnav/ephemeris.h>

#define NAV_MSG_GLO_STRING_BITS_LEN 3 /* Buffer 96 nav bits. */

/** Time mark in GLO nav string, GLO ICD, pg. 16 */
#define GLO_TM (0x3E375096)
/** Length of GLO time mark */
#define GLO_TM_LEN 30

/* States of receiver GLO bitstream */
typedef enum {
  INVALID = -1,
  SYNC_TM, /**< Time mark search */
  GET_DATA_BIT /**< Data bit receive */
} glo_receive_machine;

/* The structure is used for GLO receive and decode bitstream */
typedef struct {
  u32 string_bits[NAV_MSG_GLO_STRING_BITS_LEN]; /**< buffer for one GLO string */
  u16 current_head_bit_index; /**< how many bits written into GLO string buffer*/
  u16 nt; /**< tmp container for TOE calculation */
  u8 next_string_id; /**< what is the next string we need for parsing */
  u8 n4; /**< tmp container for TOE calculation */
  u8 hrs; /**< tmp container for TOE calculation */
  u8 min; /**< tmp container for TOE calculation */
  u8 sec; /**< tmp container for TOE calculation */
  glo_receive_machine state; /**< current state of receiver */
  u8 meander_bits_cnt:2; /**< counter for line code bits, MAX is 2 */
  u8 manchester:2; /**< 2 bits line code received */
} nav_msg_glo_t;

void nav_msg_init_glo(nav_msg_glo_t *n);
s8 process_string_glo(nav_msg_glo_t *n, ephemeris_t *e);
u32 extract_word_glo(const nav_msg_glo_t *n, u16 bit_index, u8 n_bits);
s8 nav_msg_update_glo(nav_msg_glo_t *n, bool bit_val);

#endif /* LIBSWIFTNAV_NAV_MSG_GLO_H */
