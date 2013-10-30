#include <check.h>

#include <edc.h>

const u8 *test_data = (u8*)"123456789";

START_TEST(test_crc16_ccitt)
{
  u16 crc;

  crc = crc16_ccitt(test_data, 0, 0);
  fail_unless(crc == 0,
    "CRC of empty buffer with starting value 0 should be 0, not %d", crc);

  crc = crc16_ccitt(test_data, 0, 22);
  fail_unless(crc == 22,
    "CRC of empty buffer with starting value 22 should be 22, not %d", crc);

  /* Test value taken from python crcmod package tests, see:
   * http://crcmod.sourceforge.net/crcmod.predefined.html */
  crc = crc16_ccitt(test_data, 9, 0);
  fail_unless(crc == 0x31C3,
    "CRC of \"123456789\" should be 0x31C3, not 0x%04X", crc);
}
END_TEST

START_TEST(test_crc24q)
{
  u32 crc;

  crc = crc24q(test_data, 0, 0);
  fail_unless(crc == 0,
    "CRC of empty buffer with starting value 0 should be 0, not %d", crc);

  crc = crc24q(test_data, 0, 22);
  fail_unless(crc == 22,
    "CRC of empty buffer with starting value 22 should be 22, not %d", crc);

  /* Test value taken from python crcmod package tests, see:
   * http://crcmod.sourceforge.net/crcmod.predefined.html */
  crc = crc24q(test_data, 9, 0xB704CE);
  fail_unless(crc == 0x21CF02,
    "CRC of \"123456789\" with init value 0xB704CE should be 0x21CF02, "
    "not 0x%06X", crc);
}
END_TEST

Suite* edc_suite(void)
{
  Suite *s = suite_create("Error Detection and Correction");

  TCase *tc_crc = tcase_create("CRC");
  tcase_add_test(tc_crc, test_crc16_ccitt);
  tcase_add_test(tc_crc, test_crc24q);
  suite_add_tcase(s, tc_crc);

  return s;
}

