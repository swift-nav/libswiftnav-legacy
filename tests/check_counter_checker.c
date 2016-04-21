#include <check.h>
#include <libswiftnav/counter_checker.h>
#include <stdio.h>
#include <stdlib.h>

struct recording {
  u8* samples;
  size_t size;
};

static size_t arr_offset = 0;

static size_t arr_read(void *ptr, size_t size, void *f)
{
  memcpy(ptr, (u8*) f + arr_offset, size);
  arr_offset += size;
  return size;
}

static void arr_rewind(void *f)
{
  (void)f;
  arr_offset = 0;
}

static struct recording generate_recording(u8 counter_start,
                         uint8_t set_counter(uint8_t data, uint8_t counter),
                         u8 counter_modulo,
                         size_t size)
{
  size_t i;
  u8* samples;
  u8 counter;
  struct recording recording = {NULL, 0};

  samples = malloc(size);
  if (NULL == samples) {
    return recording;
  }

  counter = counter_start;
  for (i = 0; i < size; i++) {
    samples[i] = set_counter(samples[i], counter);
    counter++;
    counter %= counter_modulo;
  }
  recording.samples = samples;
  recording.size = size;

  return recording;
}

START_TEST(test_counters_rf41)
{
  u8 data = 0xFF;
  u8 i;
  u8 counter;

  for (i = 0; i < 16; i++) {
    data = set_rf41_counter(data, i);
    counter = get_rf41_counter(data);
    fail_if(counter != i, "Unexpected counter value");
  }
  data = set_rf41_counter(data, 0);
  fail_if(0x3C != data, "Unexpected counter value");
}
END_TEST

START_TEST(test_counters_rf32)
{
  u8 data = 0xFF;
  u8 i;
  u8 counter;

  for (i = 0; i < 16; i++) {
    data = set_rf32_counter(data, i);
    counter = get_rf32_counter(data);
    fail_if(counter != i, "Unexpected counter value");
  }
  data = set_rf32_counter(data, 0);
  fail_if(0xC3 != data, "Unexpected counter value");
}
END_TEST

START_TEST(test_counter_checker_rf41)
{
  size_t i;
  struct recording recording;
  const struct mismatch_data *mismatch;
  size_t size;
  size_t disc_counter;

  static struct callbacks cbs = {
    arr_read,
    arr_rewind,
    get_rf41_counter,
    set_rf41_counter,
  };

  // check sample data with different counter values
  size = SAMPLE_CHUNK_SIZE + SAMPLE_CHUNK_SIZE / 2;
  for (i = 0; i < SAMPLE_COUNTER_MODULO; i++) {
    recording = generate_recording(i, set_rf41_counter,
                                   SAMPLE_COUNTER_MODULO,
                                   SAMPLE_CHUNK_SIZE + SAMPLE_CHUNK_SIZE / 2);
    fail_if(NULL == recording.samples, "Could not allocate recording data");
    fail_if(size != recording.size, "Recording data size mismatch");

    arr_offset = 0;
    counter_checker_init();
    mismatch = counter_checker_run(&cbs, recording.samples, recording.size);
    fail_if(NULL == mismatch, "Unexpected return value");
    fail_if(mismatch->counter != 0, "Unexpected counter mismatch found");
    free(recording.samples);
  }

  // introduce counter discontinuity
  size = 2 * SAMPLE_CHUNK_SIZE + SAMPLE_CHUNK_SIZE / 2;
  recording = generate_recording(0, set_rf41_counter, SAMPLE_COUNTER_MODULO,
                                 size);
  fail_if(NULL == recording.samples, "Could not allocate recording data");
  fail_if(size != recording.size, "Recording data size mismatch");

  disc_counter = 0;
  for (i = 1; i < size; i += size / 1024) {
    u8 counter = get_rf41_counter(recording.samples[i - 1]);

    counter += 2;
    counter %= SAMPLE_COUNTER_MODULO;
    recording.samples[i] = set_rf41_counter(recording.samples[i], counter);

    disc_counter++;
  }

  arr_offset = 0;
  counter_checker_init();
  mismatch = counter_checker_run(&cbs, recording.samples, recording.size);

  fail_if(mismatch->counter != 2 * disc_counter, "Counter mismatch failure");
}
END_TEST

START_TEST(test_counter_checker_rf32)
{
  size_t i;
  struct recording recording;
  const struct mismatch_data *mismatch;
  size_t size;
  size_t disc_counter;

  static struct callbacks cbs = {
    arr_read,
    arr_rewind,
    get_rf32_counter,
    set_rf32_counter,
  };

  // check sample data with different counter values
  size = SAMPLE_CHUNK_SIZE + SAMPLE_CHUNK_SIZE / 2;
  for (i = 0; i < SAMPLE_COUNTER_MODULO; i++) {
    recording = generate_recording(i, set_rf32_counter,
                                   SAMPLE_COUNTER_MODULO,
                                   SAMPLE_CHUNK_SIZE + SAMPLE_CHUNK_SIZE / 2);
    fail_if(NULL == recording.samples, "Could not allocate recording data");
    fail_if(size != recording.size, "Recording data size mismatch");

    arr_offset = 0;
    counter_checker_init();
    mismatch = counter_checker_run(&cbs, recording.samples, recording.size);
    fail_if(NULL == mismatch, "Unexpected return value");
    fail_if(mismatch->counter != 0, "Unexpected counter mismatch found");
    free(recording.samples);
  }

  // introduce counter discontinuity
  size = 2 * SAMPLE_CHUNK_SIZE + SAMPLE_CHUNK_SIZE / 2;
  recording = generate_recording(0, set_rf32_counter, SAMPLE_COUNTER_MODULO,
                                 size);
  fail_if(NULL == recording.samples, "Could not allocate recording data");
  fail_if(size != recording.size, "Recording data size mismatch");

  disc_counter = 0;
  for (i = 1; i < size; i += size / 1024) {
    u8 counter = get_rf32_counter(recording.samples[i - 1]);

    counter += 2;
    counter %= SAMPLE_COUNTER_MODULO;
    recording.samples[i] = set_rf32_counter(recording.samples[i], counter);

    disc_counter++;
  }

  arr_offset = 0;
  counter_checker_init();
  mismatch = counter_checker_run(&cbs, recording.samples, recording.size);

  fail_if(mismatch->counter != 2 * disc_counter, "Counter mismatch failure");
}
END_TEST

Suite* counter_checker_suite(void)
{
  Suite *s = suite_create("Counter checker");
  TCase *tc_core = tcase_create("Core");

  tcase_add_test(tc_core, test_counters_rf32);
  tcase_add_test(tc_core, test_counters_rf41);
  tcase_add_test(tc_core, test_counter_checker_rf41);
  tcase_add_test(tc_core, test_counter_checker_rf32);
  suite_add_tcase(s, tc_core);

  return s;
}
