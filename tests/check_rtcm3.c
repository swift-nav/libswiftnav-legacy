
#include <string.h>
#include <check.h>

#include <rtcm3.h>

START_TEST(test_rtcm3_check_frame)
{
  /* Test data taken from RTCM 10403.1 Document Example 4.2 */
  u8 test_data[] = {
    /* Frame Header. */
    0xD3, 0x00, 0x13,
    /* Data Message. */
    0x3E, 0xD7, 0xD3, 0x02, 0x02, 0x98, 0x0E, 0xDE, 0xEF, 0x34, 0xB4, 0xBD,
    0x62, 0xAC, 0x09, 0x41, 0x98, 0x6F, 0x33,
    /* CRC. */
    0x36, 0x0B, 0x98
  };

  /* Check against valid data. */
  s16 ret = rtcm3_check_frame(test_data);
  fail_unless(ret == 19,
      "Should return length 19 for test_data, got %d", ret);

  /* Corrupt preamble byte. */
  test_data[0] = 0x22;
  ret = rtcm3_check_frame(test_data);
  fail_unless(ret == -1,
      "Should return -1 for test_data with invalid preamble, got %d", ret);

  /* Correct preamble byte and corrupt data message. */
  test_data[0] = 0xD3;
  test_data[4] = 0x22;
  ret = rtcm3_check_frame(test_data);
  fail_unless(ret == -2,
      "Should return -2 for test_data with invalid CRC, got %d", ret);
}
END_TEST

START_TEST(test_rtcm3_write_frame)
{
  u8 test_data[] = {
    /* Frame Header. */
    0x22, 0x22, 0x22,
    /* Data Message. */
    0x3E, 0xD7, 0xD3, 0x02, 0x02, 0x98, 0x0E, 0xDE, 0xEF, 0x34, 0xB4, 0xBD,
    0x62, 0xAC, 0x09, 0x41, 0x98, 0x6F, 0x33,
    /* CRC. */
    0x22, 0x22, 0x22
  };

  /* Test data taken from RTCM 10403.1 Document Example 4.2 */
  u8 test_data_expected[] = {
    /* Frame Header. */
    0xD3, 0x00, 0x13,
    /* Data Message. */
    0x3E, 0xD7, 0xD3, 0x02, 0x02, 0x98, 0x0E, 0xDE, 0xEF, 0x34, 0xB4, 0xBD,
    0x62, 0xAC, 0x09, 0x41, 0x98, 0x6F, 0x33,
    /* CRC. */
    0x36, 0x0B, 0x98
  };

  u8 test_empty[] = {
    /* Frame Header. */
    0x22, 0x22, 0x22,
    /* CRC. */
    0x22, 0x22, 0x22
  };

  u8 test_empty_expected[] = {
    /* Frame Header. */
    0xD3, 0x00, 0x00,
    /* CRC. */
    0x47, 0xEA, 0x4B
  };

  s8 ret = rtcm3_write_frame(2222, NULL);

  fail_unless(ret == -1,
      "Should return -1 if length is larger than 1023, got %d", ret);

  fail_unless(rtcm3_write_frame(0, test_empty) == 0,
      "Returned error on test_in_empty");

  fail_unless(memcmp(test_empty, test_empty_expected, 6) == 0,
      "test_empty != test_empty_expected");

  fail_unless(rtcm3_write_frame(sizeof(test_data)-6, test_data) == 0,
      "Returned error on test_in_empty");

  fail_unless(memcmp(test_data, test_data_expected, sizeof(test_data)) == 0,
      "test_data != test_data_expected");
}
END_TEST

Suite* rtcm3_suite(void)
{
  Suite *s = suite_create("RTCMv3");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_rtcm3_check_frame);
  tcase_add_test(tc_core, test_rtcm3_write_frame);
  suite_add_tcase(s, tc_core);

  return s;
}

