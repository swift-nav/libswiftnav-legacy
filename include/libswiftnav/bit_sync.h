/*
 * Copyright (c) 2015,2016 Swift Navigation Inc.
 * Contact: Jacob McNamee <jacob@swiftnav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef LIBSWIFTNAV_BIT_SYNC_H
#define LIBSWIFTNAV_BIT_SYNC_H

#include <libswiftnav/common.h>
#include <libswiftnav/signal.h>

/** \addtogroup bit_sync
 * \{ */

#define BITSYNC_UNSYNCED -1

#define BIT_LENGTH_MAX 20

/** Structure containing bit sync state for a signal */
typedef struct {

  u8 bit_length;      /** length of a single bit */

  u8 bit_phase;
  s8 bit_phase_ref;  /**< -1 = not synced.*/
  s32 bit_integrate;

  u8 bitsync_count;
  s32 bitsync_prev_corr[BIT_LENGTH_MAX];
  u32 bitsync_histogram[BIT_LENGTH_MAX];

} bit_sync_t;

/** \} */

void bit_sync_init(bit_sync_t *b, gnss_signal_t sid);
void bit_sync_set(bit_sync_t *b, s8 bit_phase_ref);
bool bit_sync_update(bit_sync_t *b, s32 corr_prompt_real, u32 ms,
                     s32 *bit_integrate);

#endif /* LIBSWIFTNAV_BIT_SYNC_H */
