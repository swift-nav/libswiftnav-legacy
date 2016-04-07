/*
 * Copyright (C) 2016 Swift Navigation Inc.
 * Contact: Adel Mamin <adelm@exafore.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

/**
  @file

  This utility helps to verify the integrity of samples data from Piksi v3 HW.
  The regular Piksi v3 samples data format in one byte is:
  RF4 RF3 RF2 RF1
   00  00  00  00

  RF1 - GPS L1
  RF2 - GLONASS L1
  RF3 - GLONASS L2
  RF4 - GPS L2

  This utility is expecting RF3 and RF2 to have 4 bits counter,
  which is continuosly incremented with modulo SAMPLE_COUNTER_MODULO.
  RF2 is expected to have two least significant bits
  and RF3 - two most significant bits.
  Therefore, the bit indexes are:
  RF3 RF2
  5&4 3&2
*/

#ifndef LIBSWIFTNAV_COUNTER_CHECKER_H
#define LIBSWIFTNAV_COUNTER_CHECKER_H

#include <stdint.h>
#include <stddef.h>
#include <libswiftnav/common.h>

/** How many bytes we read from samples stream at a time. */
#define SAMPLE_CHUNK_SIZE (1024 * 1024) /* [bytes] */

/** The counter embedded into the sample data stream is expected
    to have this modulo. */
#define SAMPLE_COUNTER_MODULO 13

/** The counter bit size. */
#define SAMPLE_COUNTER_BITS   4 /* [bits] */

/** Teh counter bit offset within a byte. */
#define SAMPLE_COUNTER_OFFSET 2 /* [bits] */

/** How many bytes we read from samples file at a time. */
#define COUNTER_CHECKER_CHUNK_SIZE (1024 * 1024) /* [bytes] */

/** Get \e num of bits at \e offset in \e data */
#define GET_BITS(data, offset, num) \
  ( ((data) >> (offset)) & ((1 << (num)) - 1) )

/** Get \e num of bits at \e offset in \e data */
#define SET_BITS(data, offset, num, bits) \
  ( ( ~( ((1 << (num)) - 1) << (offset) ) & (data) ) | \
  ( ( (bits) & ((1 << (num)) - 1) ) << (offset) ) )

/** Get GPS counter
 * \param data Data byte
 * \return The counter value
 */
uint8_t get_gps_counter(uint8_t data);

/** Get Glonass counter
 * \param data Data byte
 * \return The counter value
 */
uint8_t get_glo_counter(uint8_t data);

uint8_t set_gps_counter(uint8_t data, uint8_t counter);
uint8_t set_glo_counter(uint8_t data, uint8_t counter);

/** A data mismatch descriptor. */
struct mismatch {
  size_t offset;                /**! Data offset [bytes]. */
  uint8_t data;                 /**! Data. */
  uint8_t expected_counter;     /**! The expected counter value. */
  uint8_t actual_counter;       /**! The actual counter value. */
};

/** Mismatch array data */
struct mismatch_data {
  /** The sample data counter mismatch incidents are stored here. */
  struct mismatch data[COUNTER_CHECKER_CHUNK_SIZE];
  /** How many valid entries there are in \e mismatch array. */
  size_t counter;
};

struct callbacks {
  size_t (*read)(void *ptr, size_t size, void *f);
  void (*rewind)(void *f);
  uint8_t (*get_counter)(uint8_t data);
  uint8_t (*set_counter)(uint8_t data, uint8_t counter);
};

typedef size_t (*stream_read_t)(void *ptr, size_t size, void *stream);
typedef void   (*stream_rewind_t)(void *stream);

uint8_t get_rf32_counter(uint8_t data);
uint8_t set_rf32_counter(uint8_t data, uint8_t counter);
uint8_t get_rf41_counter(uint8_t data);
uint8_t set_rf41_counter(uint8_t data, uint8_t counter);

void report_mismatch(const struct mismatch_data *mismatch);

void counter_checker_init(void);

const struct mismatch_data *counter_checker_run(struct callbacks *cbs,
                                               void *stream,
                                               size_t data_size);

#endif /* LIBSWIFTNAV_COUNTER_CHECKER_H */
