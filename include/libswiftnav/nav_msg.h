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

#ifndef LIBSWIFTNAV_NAV_MSG_H
#define LIBSWIFTNAV_NAV_MSG_H

#include "common.h"
#include "ephemeris.h"

#define NAV_MSG_SUBFRAME_BITS_LEN 14 /* Buffer 448 nav bits. */

#define TOW_INVALID -1
#define BITSYNC_UNSYNCED -1

#define BIT_POLARITY_NORMAL 0
#define BIT_POLARITY_INVERTED 1
#define BIT_POLARITY_UNKNOWN -1

typedef struct {
  u32 subframe_bits[NAV_MSG_SUBFRAME_BITS_LEN];
  u16 subframe_bit_index;
  bool overrun;
  /** subframe_start_index:
   * - 0 = no preamble found
   * - +x = preamble begins at bit index (x-1)
   * - -x = inverse preamble begins at (1-x)
   */
  s16 subframe_start_index;
  u8 bit_phase;
  s8 bit_phase_ref;  /**< -1 = not synced.*/
  s32 bit_integrate;

  u32 frame_words[3][8];
  u8 next_subframe_id;
  s8 bit_polarity;

  u8 bitsync_count;
  s32 bitsync_prev_corr[20];
  u32 bitsync_histogram[20];
} nav_msg_t;

void nav_msg_init(nav_msg_t *n);
s32 nav_msg_update(nav_msg_t *n, s32 corr_prompt_real, u8 ms);
bool subframe_ready(nav_msg_t *n);
s8 process_subframe(nav_msg_t *n, ephemeris_t *e);

#endif /* LIBSWIFTNAV_NAV_MSG_H */

