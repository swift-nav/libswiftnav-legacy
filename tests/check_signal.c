#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "check_utils.h"

#include <libswiftnav/signal.h>
#include <libswiftnav/constants.h>

/* Include singal.c here to have a chance to turn off asserts.
   Otherwise the code lines with asserts cannot be covered by
   tests and it will lower test coverage statistics. */
#define NDEBUG
#include <signal.c>

#define ARRAY_COUNT(arr) ((sizeof(arr) / sizeof(arr[0])))

static const struct code_data_element {
  enum code code;
  u16 sat_count;
} code_data[] = {
  {CODE_GPS_L1CA, NUM_SIGNALS_GPS_L1CA},
  {CODE_GPS_L2CM, NUM_SIGNALS_GPS_L2CM},
  {CODE_SBAS_L1CA, NUM_SIGNALS_SBAS_L1CA},
  {CODE_GLO_L1CA, NUM_SIGNALS_GLO_L1CA},
  {CODE_GLO_L2CA, NUM_SIGNALS_GLO_L2CA},
};

static const struct constellation_data_element {
  enum constellation constellation;
  u16 sat_count;
  u8 code_count;
} constellation_data[] = {
  {CONSTELLATION_GPS, NUM_SIGNALS_GPS, NUM_CODES_GPS},
  {CONSTELLATION_SBAS, NUM_SIGNALS_SBAS, NUM_CODES_SBAS},
  {CONSTELLATION_GLO, NUM_SIGNALS_GLO, NUM_CODES_GLO},
};

START_TEST(test_signal_aggregates)
{
  fail_unless(ARRAY_COUNT(code_data) == CODE_COUNT,
              "missing code entry in code_data[]");

  fail_unless(ARRAY_COUNT(constellation_data) == CONSTELLATION_COUNT,
              "missing constellation entry in constellation_data[]");

  u8 constellation_code_counts[CONSTELLATION_COUNT];
  memset(constellation_code_counts, 0, sizeof(constellation_code_counts));
  for (u32 i=0; i<ARRAY_COUNT(code_data); i++) {
    const struct code_data_element *e = &code_data[i];
    enum constellation constellation = code_to_constellation(e->code);
    constellation_code_counts[constellation]++;
  }
  for (u32 i=0; i<ARRAY_COUNT(constellation_data); i++) {
      const struct constellation_data_element *e = &constellation_data[i];
      fail_unless(e->code_count == constellation_code_counts[e->constellation],
                  "invalid code count definition for code %d",
                  e->constellation);
  }
}
END_TEST

START_TEST(test_signal_from_index)
{
  for (u32 i=0; i<ARRAY_COUNT(code_data); i++) {
    const struct code_data_element *e = &code_data[i];
    enum code code = e->code;
    u16 sat_count = e->sat_count;
    for (u16 code_index = 0; code_index < sat_count; code_index++) {
      gnss_signal_t sid = sid_from_code_index(code, code_index);
      fail_unless(sid_valid(sid), "signal from code index not valid: "
                                  "code %d code index %d",
                                  code, code_index);
      fail_unless(sid_to_code_index(sid) == code_index,
                  "signal from code index to code index failed: "
                  "code %d code index %d",
                  code, code_index);
    }
  }
}
END_TEST

START_TEST(test_signal_properties)
{
  const struct test_case
  {
    gnss_signal_t sid;
    bool valid;
    char str[SID_STR_LEN_MAX];
  } test_cases[] = {
      {
        .sid = {
          .code = CODE_INVALID, .sat = 0
        },
        .valid = false
      },
      {
        .sid = {
          .code = CODE_COUNT, .sat = 0
        },
        .valid = false
      },
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
      {
        .sid = {
          .code = CODE_GLO_L1CA, .sat = 0
        },
        .valid = false,
      },
      {
        .sid = {
          .code = CODE_GLO_L1CA, .sat = 1
        },
        .valid = true,
        .str = "GLO L1CA 1"
      },
      {
        .sid = {
          .code = CODE_GLO_L1CA, .sat = 24
        },
        .valid = true,
        .str = "GLO L1CA 24"
      },
      {
        .sid = {
          .code = CODE_GLO_L1CA, .sat = 25
        },
        .valid = false
      },
      {
        .sid = {
          .code = CODE_GLO_L2CA, .sat = 0
        },
        .valid = false,
      },
      {
        .sid = {
          .code = CODE_GLO_L2CA, .sat = 1
        },
        .valid = true,
        .str = "GLO L2CA 1"
      },
      {
        .sid = {
          .code = CODE_GLO_L2CA, .sat = 24
        },
        .valid = true,
        .str = "GLO L2CA 24"
      },
      {
        .sid = {
          .code = CODE_GLO_L2CA, .sat = 25
        },
        .valid = false
      },
  };

  for (u32 i=0; i<sizeof(test_cases) / sizeof(test_cases[0]); i++) {
    const struct test_case *t = &test_cases[i];
    bool valid = sid_valid(t->sid);
    fail_unless(t->valid == valid, "test signal %d validity incorrect", i);
    if (valid) {
      gnss_signal_t sid =
          sid_from_code_index(t->sid.code, sid_to_code_index(t->sid));
      fail_unless(sid_is_equal(t->sid, sid),
                  "test signal %d code index conversion failed");
      char str[SID_STR_LEN_MAX] = {0};
      u32 ret = sid_to_string(str, sizeof(str), sid);
      fail_unless((strcmp(str, t->str) == 0) && (ret == strlen(t->str)),
                  "signal to string \"%s\" failed", t->str);
      fail_unless(constellation_valid(code_to_constellation(sid.code)),
                  "constellation_valid failed for code %d", sid.code);
    }
    fail_unless(!constellation_valid(CONSTELLATION_COUNT),
                  "constellation_valid failed for constellation %d",
                   CONSTELLATION_COUNT);
  }
}
END_TEST

START_TEST(test_signal_compare)
{
  gnss_signal_t sids[NUM_SIGNALS];
  u32 signal_index = 0;

  for (u32 i=0; i<ARRAY_COUNT(code_data); i++) {
    const struct code_data_element *e = &code_data[i];
    enum code code = e->code;
    u16 sat_count = e->sat_count;
    for (u16 code_index = 0; code_index < sat_count; code_index++) {
      sids[signal_index++] = sid_from_code_index(code, code_index);
    }
  }

  qsort(sids, NUM_SIGNALS, sizeof(gnss_signal_t), cmp_sid_sid);

  for (u32 i=1; i<NUM_SIGNALS; i++) {
    fail_unless(!sid_is_equal(sids[i], sids[i-1]),
                "signal index %d not unique", i);
    fail_unless(sid_compare(sids[i], sids[i-1]) > 0,
                "signal index %d not in order", i);
  }
}
END_TEST

START_TEST(test_signal_construction)
{
  for (u32 i=0; i<ARRAY_COUNT(code_data); i++) {
    const struct code_data_element *e = &code_data[i];
    enum code code = e->code;
    u16 sat_count = e->sat_count;
    for (u16 code_index = 0; code_index < sat_count; code_index++) {
      gnss_signal_t sid = sid_from_code_index(code, code_index);
      gnss_signal_t csid = construct_sid(sid.code, sid.sat);
      fail_unless(sid_valid(csid),
                  "constructed signal not valid: "
                  "code %d code index %d",
                  code, code_index);
      fail_unless(sid_is_equal(sid, csid),
                  "constructed signal mismatch: "
                  "code %d code index %d",
                   code, code_index);
    }
  }
}
END_TEST

START_TEST(test_signal_code_to_carr_freq)
{
  double carr_freq;

  carr_freq = code_to_carr_freq(CODE_GPS_L2CM);
  fail_unless(GPS_L2_HZ == carr_freq);

  /* check unsupported branch for code coverage stats */
  carr_freq = code_to_carr_freq(CODE_GLO_L2CA);
}
END_TEST

START_TEST(test_signal_code_to_lambda)
{
  double lambda;

  lambda = code_to_lambda(CODE_GPS_L2CM);
  fail_unless((GPS_C / GPS_L2_HZ) == lambda);

  lambda = code_to_lambda(CODE_GPS_L1CA);
  fail_unless((GPS_C / GPS_L1_HZ) == lambda);

  /* check unsupported branch for code coverage stats */
  lambda = code_to_lambda(CODE_GLO_L2CA);
}
END_TEST

START_TEST(test_signal_code_to_chip_count)
{
  u16 chip_count;

  chip_count = code_to_chip_count(CODE_GPS_L1CA);
  fail_unless(GPS_L1CA_CHIPS_NUM == chip_count);

  chip_count = code_to_chip_count(CODE_SBAS_L1CA);
  fail_unless(GPS_L1CA_CHIPS_NUM == chip_count);

  chip_count = code_to_chip_count(CODE_GPS_L2CM);
  fail_unless(GPS_L2CM_CHIPS_NUM * 2 == chip_count);

  /* check unsupported branch for code coverage stats */
  chip_count = code_to_chip_count(CODE_GLO_L2CA);
}
END_TEST

START_TEST(test_signal_code_to_chip_rate)
{
  double chip_rate;

  chip_rate = code_to_chip_rate(CODE_GPS_L1CA);
  fail_unless(GPS_CA_CHIPPING_RATE == chip_rate);

  chip_rate = code_to_chip_rate(CODE_SBAS_L1CA);
  fail_unless(GPS_CA_CHIPPING_RATE == chip_rate);

  chip_rate = code_to_chip_rate(CODE_GPS_L2CM);
  fail_unless(GPS_CA_CHIPPING_RATE == chip_rate);

  /* check unsupported branch for code coverage stats */
  chip_rate = code_to_chip_rate(CODE_GLO_L2CA);
}
END_TEST

START_TEST(test_signal_code_requires_direct_acq)
{
  bool requires;

  requires = code_to_chip_rate(CODE_GPS_L1CA);
  fail_unless(true == requires);

  requires = code_requires_direct_acq(CODE_GPS_L2CM);
  fail_unless(false == requires);
}
END_TEST

Suite* signal_test_suite(void)
{
  Suite *s = suite_create("Signal");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_signal_aggregates);
  tcase_add_test(tc_core, test_signal_from_index);
  tcase_add_test(tc_core, test_signal_properties);
  tcase_add_test(tc_core, test_signal_compare);
  tcase_add_test(tc_core, test_signal_construction);
  tcase_add_test(tc_core, test_signal_code_to_carr_freq);
  tcase_add_test(tc_core, test_signal_code_to_chip_count);
  tcase_add_test(tc_core, test_signal_code_to_chip_rate);
  tcase_add_test(tc_core, test_signal_code_requires_direct_acq);
  tcase_add_test(tc_core, test_signal_code_to_lambda);
  suite_add_tcase(s, tc_core);

  return s;
}

