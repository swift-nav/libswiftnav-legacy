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
  The regular Piksi v3 samples data format is one byte, which is encoded
  lile this:

  RF4 RF3 RF2 RF1
   00  00  00  00

  , where

  RF1 - GPS L1
  RF2 - GLONASS L1
  RF3 - GLONASS L2
  RF4 - GPS L2

  This utility expects RF3 and RF2 to have 4 bits counter,
  which is continuosly incremented with modulo SAMPLE_COUNTER_MODULO.
  RF2 is expected to have two least significant bits of the counter
  and RF3 - two most significant bits.
  Therefore, the bit indexes are:
  RF3 RF2
  5&4 3&2
*/

#include <stdio.h>
#include <stdint.h>
#include <libswiftnav/counter_checker.h>

/** The data is read in chunks, which are stored here. */
static uint32_t chunk[SAMPLE_CHUNK_SIZE / sizeof(uint32_t)];

/** Get array size */
#define ARR_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))

/** Mismatch data */
static struct mismatch_data mismatch;

/** Checks counter mismatch.
 *
 * \param data_offset File offset in bytes.
 * \param bytes The data bytes to check.
 * \param byte_index The byte index within \e bytes parameter [0,1,2,3]
 * \param expected_counter Next expected counter value.
 *
 * \return The next expected counter value
 */
static inline uint8_t check_counter(size_t data_offset,
                                    uint8_t (*get_counter)(uint8_t data),
                                    uint32_t bytes, uint8_t byte_index,
                                    uint32_t expected_counter)
{
  uint8_t actual_counter;
  uint8_t byte_offset = 8 * byte_index;

  actual_counter = get_counter(bytes >> byte_offset);
  if (actual_counter != expected_counter) {
    if (mismatch.counter <= ARR_SIZE(mismatch.data)) {
      mismatch.data[mismatch.counter].offset = data_offset;
      mismatch.data[mismatch.counter].data = (bytes >> 8 * byte_index) & 0xFF;
      mismatch.data[mismatch.counter].expected_counter = expected_counter;
      mismatch.data[mismatch.counter].actual_counter = actual_counter;
      mismatch.counter++;
    }
  }
  expected_counter = (actual_counter + 1) % SAMPLE_COUNTER_MODULO;

  return expected_counter;
}

/** Report counter mismatch issues. */
void report_mismatch(const struct mismatch_data *mismatch)
{
  size_t i;
  size_t offset_prev = mismatch->data[0].offset;

  if (0 != mismatch->counter) {
    printf("data_offset,data,expected_counter,actual_counter,offset_diff\n");
  }
  for (i = 0; i < mismatch->counter; i++) {
    printf("%ld,", mismatch->data[i].offset);
    printf("0x%02X,", mismatch->data[i].data);
    printf("%d,", mismatch->data[i].expected_counter);
    printf("%d,", mismatch->data[i].actual_counter);
    printf("%ld\n", mismatch->data[i].offset - offset_prev);
    offset_prev = mismatch->data[i].offset;
  }
}

/** Get RF3 RF2 counter
 * \param data Data byte
 * \return The counter value
 */
uint8_t get_rf32_counter(uint8_t data)
{
  uint8_t counter = GET_BITS(data, 2, 4);
  return counter;
}

/** Set RF3 RF2 counter
 * \param data Data byte
 * \param counter The counter value
 * \return The data byte with embedded counter
 */
uint8_t set_rf32_counter(uint8_t data, uint8_t counter)
{
  data = SET_BITS(data, 2, 4, counter);
  return data;
}

/** Get RF4 RF1 counter
 * \param data Data byte
 * \return The counter value
 */
uint8_t get_rf41_counter(uint8_t data)
{
  uint8_t counter;

  // take RF4 part first
  counter = GET_BITS(data, 6, 2) << 2;
  // combine with RF1 part to get full counter
  counter |= GET_BITS(data, 0, 2);

  return counter;
}

/** Set RF4 RF1 counter
 * \param data Data byte
 * \param counter The counter value
 * \return The data byte with embedded counter
 */
uint8_t set_rf41_counter(uint8_t data, uint8_t counter)
{
  data = SET_BITS(data, 0, 2, counter);
  counter >>= 2;
  data = SET_BITS(data, 6, 2, counter);

  return data;
}

/** Initialize sample checker */
void counter_checker_init(void)
{
  mismatch.counter = 0;
}

/** Runs sample checker
 *
 * \param csb Callbacks for data reading and counter extraction
 * \param stream Stream itself.
 * \param data_size Total data size [bytes]
 * \return Mismatch data pointer. Non NULL even if no mismatch found.
           The memory pointer is owned by the callee. Do not free it!
 */
const struct mismatch_data *counter_checker_run(struct callbacks *cbs,
                                               void *stream,
                                               size_t data_size)
{
  size_t i;
  size_t j;
  size_t chunk_size;
  size_t chunk_num;
  uint8_t tmp;
  uint8_t counter;
  size_t data_offset;

  /* initialize the counter from the first byte of the samples data */
  chunk_size = cbs->read(&tmp, 1, stream);
  counter = cbs->get_counter(tmp);
  cbs->rewind(stream);

  data_offset = 0;

  /* now let's start the integrity check of the samples data
     in chunks of size SAMPLE_CHUNK_SIZE */
  chunk_num = data_size / sizeof(chunk);

  /* printf("chunk_num = %ld\n", chunk_num); */

  for (i = 0; i < chunk_num; i++) {
    chunk_size = cbs->read(chunk, sizeof(chunk), stream);

    /* process the chunk in 4 bytes slices */
    for (j = 0; j < chunk_size / sizeof(chunk[0]); j++) {
      uint8_t k;

      for (k = 0; k < sizeof(uint32_t); k++) {
        counter = check_counter(data_offset, cbs->get_counter,
                                chunk[j], k, counter);
        data_offset++;
      }
    }
  }

  chunk_size = data_size % sizeof(chunk);
  chunk_size = cbs->read(chunk, chunk_size, stream);

  /* printf("chunk_size = %ld\n", chunk_size); */

  /* process the last chunk byte by byte */
  for (j = 0; j < chunk_size; j++) {
    counter = check_counter(data_offset, cbs->get_counter,
                            ((uint8_t*)chunk)[j], 0, counter);
    data_offset++;
  }

  return &mismatch;
}

#ifdef STANDALONE_SETUP

/** Standard C library 'fread' function wrapper  */
static size_t file_read(void *ptr, size_t size, void *f)
{
  size_t ret;
  ret = fread(ptr, size, 1, (FILE*)f);
  return ret * size;
}

/** Standard C library 'rewind' function wrapper */
static void file_rewind(void *f)
{
  rewind((FILE*)f);
}

static struct callbacks cbs = {

  file_read,
  file_rewind,

#if defined RF41
  get_rf41_counter,
  set_rf41_counter,
#elif defined RF32
  get_rf32_counter,
  set_rf32_counter
#else
# error "Please, specify, where counter bits are located."
#endif

};

int main(int argc, char* argv[])
{
  FILE *f;
  const char* file_name;
  long file_size;
  int ret = 0;
  const struct mismatch_data *mismatch;


  if (argc < 2) {
    printf("Usage:");
    printf("%s <samples data file name>\n", argv[0]);
    return -1;
  }

  file_name = argv[1];

  f = fopen(file_name, "rb");
  if (NULL == f) {
    printf("Failed to open file %s\n", file_name);
    return -1;
  }

  fseek(f, 0, SEEK_END);
  file_size = ftell(f);
  rewind(f);

  if (0 == file_size) {
    printf("File %s is empty\n", file_name);
    ret = -1;
    goto end;
  }

  counter_checker_init();
  mismatch = counter_checker_run(&cbs, f, file_size);

  if (0 == mismatch->counter) {
    printf("The file is OK!\n");
  } else {
    report_mismatch(mismatch);
  }


end:

  fclose(f);

  return ret;
}

#endif  /* #ifdef STANDALONE_SETUP */
