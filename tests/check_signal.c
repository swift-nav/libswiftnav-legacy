#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "check_utils.h"

#include <libswiftnav/signal.h>

START_TEST(test_signal_from_index)
{
  for (u32 i=0; i<NUM_SIGNALS; i++) {
    gnss_signal_t sid = sid_from_index(i);
    fail_unless(sid_valid(sid), "signal from index %d not valid", i);
    fail_unless(sid_to_index(sid, INDEX_GLOBAL) == i,
                "signal from index %d to index failed", i);
  }
  for (u32 i=0; i<NUM_SIGNALS_GPS; i++) {
    gnss_signal_t sid = sid_from_index(i);
    fail_unless(sid_valid(sid), "signal from index %d not valid", i);
    fail_unless(sid_to_index(sid, INDEX_CONSTELLATION) == i,
                "signal from index %d to index failed", i);
  }
  for (u32 i=0; i<NUM_SIGNALS_SBAS; i++) {
    gnss_signal_t sid = sid_from_index(NUM_SIGNALS_GPS + i);
    fail_unless(sid_valid(sid), "signal from index %d not valid", i);
    fail_unless(sid_to_index(sid, INDEX_CONSTELLATION) == i,
                "signal from index %d to index failed", i);
  }
  for (u32 i=0; i<NUM_SIGNALS_GPS; i++) {
    gnss_signal_t sid = sid_from_index(i);
    fail_unless(sid_valid(sid), "signal from index %d not valid", i);
    fail_unless(sid_to_index(sid, INDEX_SAT_IN_CONS) == i / 2,
                "signal from index %d to index failed", i);
  }
  for (u32 i=0; i<NUM_SIGNALS_SBAS; i++) {
    gnss_signal_t sid = sid_from_index(NUM_SIGNALS_GPS + i);
    fail_unless(sid_valid(sid), "signal from index %d not valid", i);
    fail_unless(sid_to_index(sid, INDEX_SAT_IN_CONS) == i,
                "signal from index %d to index failed: %d", i);
  }
}
END_TEST

START_TEST(test_signal_to_index)
{
  const struct
  {
    gnss_signal_t sid;
    bool valid;
    char str[SID_STR_LEN_MAX];
  } test_cases[] = {
      {
      .sid = {
          .code = CODE_INVALID, .sat = 1
        },
        .valid = false
      },
      {
        .sid = {
          .code = CODE_GPS_L1CA, .sat = 0
        },
        .valid = false,
      },
      {
        .sid = {
          .code = CODE_GPS_L1CA, .sat = 1
        },
        .valid = true,
        .str = "GPS L1CA 1"
      },
      {
        .sid = {
          .code = CODE_GPS_L2CM, .sat = 1
        },
        .valid = true,
        .str = "GPS L2CM 1"
      },
      {
        .sid = {
          .code = CODE_SBAS_L1CA, .sat = 1
        },
        .valid = false
      },
      {
        .sid = {
          .code = CODE_GPS_L1CA, .sat = 32
        },
        .valid = true,
        .str = "GPS L1CA 32"
      },
      {
        .sid = {
          .code = CODE_GPS_L1CA, .sat = 33
        },
        .valid = false
      },
      {
        .sid = {
          .code = CODE_SBAS_L1CA, .sat = 0
        },
        .valid = false
      },
      {
        .sid = {
          .code = CODE_SBAS_L1CA, .sat = 119
        },
        .valid = false
      },
      {
        .sid = {
          .code = CODE_SBAS_L1CA, .sat = 120
        },
        .valid = true,
        .str = "SBAS L1CA 120"
      },
      {
        .sid = {
          .code = CODE_GPS_L1CA, .sat = 120
        },
        .valid = false
      },
      {
        .sid = {
          .code = CODE_SBAS_L1CA, .sat = 138
        },
        .valid = true,
        .str = "SBAS L1CA 138"
      },
      {
        .sid = {
          .code = CODE_SBAS_L1CA, .sat = 139
        },
        .valid = false
      },
  };

  fail_unless(sid_to_index(construct_sid(0, 0), -1) == (u16)INVALID_SID_INDEX,
                    "sid_to_index INVALID_SID_INDEX test failed");

  for (u32 i=0; i<sizeof(test_cases) / sizeof(test_cases[0]); i++) {
    bool valid = sid_valid(test_cases[i].sid);
    fail_unless(test_cases[i].valid == valid,
                "test signal %d validity incorrect", i);
    if (valid) {
      gnss_signal_t sid =
          sid_from_index(sid_to_index(test_cases[i].sid, INDEX_GLOBAL));
      char str[SID_STR_LEN_MAX] = {0};
      u32 ret = sid_to_string(str, SID_STR_LEN_MAX, sid);
      fail_unless(!strcmp(str,
          test_cases[i].str) && ret == strlen(test_cases[i].str),
          "signal to string \"%s\" failed", test_cases[i].str);
      fail_unless(sid_is_equal(test_cases[i].sid, sid),
                  "test signal %d index conversion failed");
    }
  }
}
END_TEST

START_TEST(test_signal_compare)
{
  gnss_signal_t sids[NUM_SIGNALS];
  for (u32 i=0; i<NUM_SIGNALS; i++) {
    sids[i] = sid_from_index(i);
  }

  qsort(sids, NUM_SIGNALS, sizeof(gnss_signal_t), cmp_sid_sid);

  for (u32 i=0; i<NUM_SIGNALS-1; i++) {
    fail_unless(!sid_is_equal(sids[i+1], sids[i]),
                "signal index %d not unique", i);
    fail_unless(sid_compare(sids[i+1], sids[i]) > 0,
                "signal index %d not in order", i);
  }
}
END_TEST

START_TEST(test_signal_construction)
{
  for (u32 i=0; i<NUM_SIGNALS; i++) {
    gnss_signal_t sid = sid_from_index(i);
    fail_unless(sid_valid(sid), "signal from index %d not valid", i);
    gnss_signal_t csid = construct_sid(sid.code, sid.sat);
    fail_unless(sid_valid(csid), "constructed signal from index %d not valid",
        i);
    fail_unless(!sid_compare(sid, csid),
                "signal from index %d to index failed", i);
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
  tcase_add_test(tc_core, test_signal_construction);
  suite_add_tcase(s, tc_core);

  return s;
}

