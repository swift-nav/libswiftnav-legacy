#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include "check_utils.h"

#include <libswiftnav/signal.h>

START_TEST(test_signal_from_index)
{
  for (u32 i=0; i<NUM_SATS; i++) {
    gnss_signal_t sid = sid_from_index(i);
    fail_unless(sid_valid(sid), "signal from index %d not valid", i);
    fail_unless(sid_to_index(sid) == i,
                "signal from index %d to index failed", i);
  }
}
END_TEST

START_TEST(test_signal_to_index)
{
  const struct
  {
    gnss_signal_t sid;
    bool valid;
  } test_cases[] = {
      {
        .sid = {
          .constellation = CONSTELLATION_GPS, .band = BAND_L1, .sat = 0
        },
        .valid = false
      },
      {
        .sid = {
          .constellation = CONSTELLATION_GPS, .band = BAND_L1, .sat = 1
        },
        .valid = true
      },
      {
        .sid = {
          .constellation = CONSTELLATION_GPS, .band = BAND_L1, .sat = 32
        },
        .valid = true
      },
      {
        .sid = {
          .constellation = CONSTELLATION_GPS, .band = BAND_L1, .sat = 33
        },
        .valid = false
      },
      {
        .sid = {
          .constellation = CONSTELLATION_SBAS, .band = BAND_L1, .sat = 0
        },
        .valid = false
      },
      {
        .sid = {
          .constellation = CONSTELLATION_SBAS, .band = BAND_L1, .sat = 120
        },
        .valid = true
      },
      {
        .sid = {
          .constellation = CONSTELLATION_SBAS, .band = BAND_L1, .sat = 138
        },
        .valid = true
      },
      {
        .sid = {
          .constellation = CONSTELLATION_SBAS, .band = BAND_L1, .sat = 139
        },
        .valid = false
      },
  };

  for (u32 i=0; i<sizeof(test_cases) / sizeof(test_cases[0]); i++) {
    bool valid = sid_valid(test_cases[i].sid);
    fail_unless(test_cases[i].valid == valid,
                "test signal %d validity incorrect", i);
    if (valid) {
      gnss_signal_t sid = sid_from_index(sid_to_index(test_cases[i].sid));
      fail_unless(sid_is_equal(test_cases[i].sid, sid),
                  "test signal %d index conversion failed");
    }
  }
}
END_TEST

START_TEST(test_signal_compare)
{
  gnss_signal_t sids[NUM_SATS];
  for (u32 i=0; i<NUM_SATS; i++) {
    sids[i] = sid_from_index(i);
  }

  qsort(sids, NUM_SATS, sizeof(gnss_signal_t), cmp_sid_sid);

  for (u32 i=0; i<NUM_SATS-1; i++) {
    fail_unless(!sid_is_equal(sids[i+1], sids[i]),
                "signal index %d not unique", i);
    fail_unless(sid_compare(sids[i+1], sids[i]) > 0,
                "signal index %d not in order", i);
  }
}
END_TEST

Suite* signal_test_suite(void)
{
  Suite *s = suite_create("Signal");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_signal_from_index);
  tcase_add_test(tc_core, test_signal_to_index);
  tcase_add_test(tc_core, test_signal_compare);
  suite_add_tcase(s, tc_core);

  return s;
}

