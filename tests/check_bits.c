
#include <check.h>

#include <libswiftnav/bits.h>

#include <stdio.h>

START_TEST(test_parity)
{
  fail_unless(parity(0x00000000) == 0);
  fail_unless(parity(0xFFFFFFFF) == 0);
  fail_unless(parity(0x01010101) == 0);
  fail_unless(parity(0x10101010) == 0);
  fail_unless(parity(0x10A010A0) == 0);

  fail_unless(parity(0x10000000) == 1);
  fail_unless(parity(0x00000001) == 1);
  fail_unless(parity(0x70707000) == 1);
  fail_unless(parity(0x0B0B0B00) == 1);
  fail_unless(parity(0x00E00000) == 1);
}
END_TEST

START_TEST(test_getbitu)
{
  u8 test_data[] = {
    0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF
  };

  u32 ret;

  ret = getbitu(test_data, 0, 8);
  fail_unless(ret == 0x01,
      "Test case 1 expected 0x01, got 0x%02X", ret);

  ret = getbitu(test_data, 4, 8);
  fail_unless(ret == 0x12,
      "test case 2 expected 0x12, got 0x%02X", ret);

  ret = getbitu(test_data, 28, 16);
  fail_unless(ret == 0x789A,
      "test case 3 expected 0x789A, got 0x%04X", ret);

  ret = getbitu(test_data, 12, 32);
  fail_unless(ret == 0x3456789A,
      "test case 4 expected 0x3456789A, got 0x%08X", ret);

  ret = getbitu(test_data, 10, 3);
  fail_unless(ret == 0x4,
      "test case 5 expected 0x4, got 0x%01X", ret);

  ret = getbitu(test_data, 10, 13);
  fail_unless(ret == 0x11A2,
      "test case 6 expected 0x11A2, got 0x%04X", ret);
}
END_TEST

START_TEST(test_getbits)
{
  u8 test_data[] = {
    0x00, 0x03, 0x80, 0xFF, 0xFF, 0xFF, 0xFF
  };

  s32 ret;

  ret = getbits(test_data, 0, 8);
  fail_unless(ret == 0,
      "Test case 1 expected 0, got %d", ret);

  ret = getbits(test_data, 13, 3);
  fail_unless(ret == 3,
      "Test case 2 expected 3, got %d", ret);

  ret = getbits(test_data, 14, 3);
  fail_unless(ret == -1,
      "Test case 3 expected -1, got %d", ret);

  ret = getbits(test_data, 14, 4);
  fail_unless(ret == -2,
      "Test case 4 expected -2, got %d", ret);

  ret = getbits(test_data, 24, 32);
  fail_unless(ret == -1,
      "Test case 5 expected -1, got %d", ret);
}
END_TEST

START_TEST(test_setbitu)
{
  u8 test_data[10];

  u32 ret;

  setbitu(test_data, 10, 13, 0x11A2);
  ret = getbitu(test_data, 10, 13);
  fail_unless(ret == 0x11A2,
      "test case 1 expected 0x11A2, got 0x%04X", ret);

  /* TODO: Check that setbitu doesn't write to bits other than those in the bit
   * field. */
}
END_TEST

START_TEST(test_setbits)
{
  u8 test_data[10];

  s32 ret;

  setbits(test_data, 14, 3, -1);
  ret = getbits(test_data, 14, 3);
  fail_unless(ret == -1,
      "Test case 1 expected -1, got %d", ret);

  setbits(test_data, 14, 8, 22);
  ret = getbits(test_data, 14, 8);
  fail_unless(ret == 22,
      "Test case 2 expected 22, got %d", ret);

  setbits(test_data, 24, 32, -1);
  ret = getbits(test_data, 24, 32);
  fail_unless(ret == -1,
      "Test case 3 expected -1, got %d", ret);
}
END_TEST

START_TEST(test_bitshl)
{
  u8 src0[] = { 0xDE, 0xAD, 0xBE, 0xEF };
  u8 res0[] = { 0xBE, 0xEF, 0x00, 0x00 };

  u8 src1[] = { 0xDE, 0xAD, 0xBE, 0xEF };
  u8 res1[] = { 0xEA, 0xDB, 0xEE, 0xF0 };

  u8 src2[] = { 0xDE, 0xAD, 0xBE, 0xEF };
  u8 res2[] = { 0xDB, 0xEE, 0xF0, 0x00 };

  u8 src3[] = { 0xDE, 0xAD, 0xBE, 0xEF };
  u8 res3[] = { 0xB6, 0xFB, 0xBC, 0x00 };

  bitshl(src0, sizeof(src0), 16);
  fail_unless(0 == memcmp(src0, res0, 4), "Byte shift test");

  bitshl(src1, sizeof(src1), 4);
  fail_unless(0 == memcmp(src1, res1, 4), "4-bit shift");

  bitshl(src2, sizeof(src2), 12);
  fail_unless(0 == memcmp(src2, res2, 4), "12-bit shift");

  bitshl(src3, sizeof(src3), 10);
  fail_unless(0 == memcmp(src3, res3, 4), "10-bit shift");
}
END_TEST

START_TEST(test_bitshl2)
{
  u8 src0[] = { 0xDE, 0xAD, 0xBE, 0xEF };
  u8 res0[] = { 0x00, 0x00, 0x00, 0x00 };

  bitshl(src0, sizeof(src0), 64);
  fail_unless(0 == memcmp(src0, res0, 4), "Byte shift test");
}
END_TEST

START_TEST(test_bitcopy)
{
  u8 src0[] = { 0xDE, 0xAD, 0xBE, 0xEF };
  u8 res0[] = { 0xBE, 0xEF, 0xBE, 0xEF };

  u8 src1[] = { 0xDE, 0xAD, 0xBE, 0xEF };
  u8 res1[] = { 0xEA, 0xDB, 0xEE, 0xFF };

  u8 src2[] = { 0xDE, 0xAD, 0xBE, 0xEF, 0xDE, 0xAD, 0xBE, 0xEF, 0xDE, 0xAD};
  // u8 dst2[] = { 0xDE, 0xAD, 0xBE, 0xEF, 0xDE, 0xAD, 0xBE, 0xEF, 0xDE, 0xAD};
  u8 res2[] = { 0xAD, 0xBE, 0xEF, 0xDE, 0xAD, 0xBE, 0xEF, 0xDE, 0xAD, 0xAD};

  bitcopy(src0, 0, src0, 16, 16);
  fail_unless(0 == memcmp(src0, res0, 4), "16-bit copy");

  bitcopy(src1, 0, src1, 4, 28);
  fail_unless(0 == memcmp(src1, res1, 4), "28-bit copy");

  bitcopy(src2, 0, src2, 8, 72);
  fail_unless(0 == memcmp(src2, res2, 4), "72-bit copy");
}
END_TEST

START_TEST(test_count_bits_x)
{
  u8  src8[]  = { 0xDE, 0xAD, 0x12, 0xEF };
  u8  res8[]  = { 6, 5, 2, 7 };

  u16 src16[] = { 0xDE05, 0xADF6, 0xBE32, 0xEF45 };
  u8  res16[]  = { 8, 11, 9, 10 };

  u32 src32[] = { 0xDE051234, 0x00000000, 0x00329300, 0x1F45A6C8 };
  u8  res32[]  = { 13, 0, 7, 15 };

  u64 src64[] = { 0xDE051234432150ED, 0x0000000080000000,
                  0x0032930000392300, 0x10F14350A060C080 };
  u8  res64[]  = { 26, 1, 14, 18 };

  for (unsigned i = 0; i < sizeof(src8)/ sizeof(src8[0]); i++) {
    u8 r = count_bits_u8(src8[i], 1);
    fail_unless(res8[i] == r, "count_bits_u8(0x%x, 1) = %d", src8[i], r);
    r = count_bits_u8(src8[i], 0);
    fail_unless(8-res8[i] == r, "count_bits_u8(0x%x, 0) = %d", src8[i], r);
  }

  for (unsigned i = 0; i < sizeof(src16)/ sizeof(src16[0]); i++) {
    u8 r = count_bits_u16(src16[i], 1);
    fail_unless(res16[i] == r, "count_bits_u16(0x%x, 1) = %d", src16[i], r);
    r = count_bits_u16(src16[i], 0);
    fail_unless(16-res16[i] == r, "count_bits_u16(0x%x, 0) = %d", src16[i], r);
  }

  for (unsigned i = 0; i < sizeof(src32)/ sizeof(src32[0]); i++) {
    u8 r = count_bits_u32(src32[i], 1);
    fail_unless(res32[i] == r, "count_bits_u32(0x%x, 1) = %d", src32[i], r);
    r = count_bits_u32(src32[i], 0);
    fail_unless(32-res32[i] == r, "count_bits_u32(0x%x, 0) = %d", src32[i], r);
  }

  for (unsigned i = 0; i < sizeof(src64)/ sizeof(src64[0]); i++) {
    u8 r = count_bits_u64(src64[i], 1);
    fail_unless(res64[i] == r, "count_bits_u64(0x%lx, 1) = %d", src64[i], r);
    r = count_bits_u64(src64[i], 0);
    fail_unless(64-res64[i] == r, "count_bits_u64(0x%lx, 0) = %d", src64[i], r);
  }
}
END_TEST

Suite* bits_suite(void)
{
  Suite *s = suite_create("Bit Utils");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_parity);
  tcase_add_test(tc_core, test_getbitu);
  tcase_add_test(tc_core, test_getbits);
  tcase_add_test(tc_core, test_setbitu);
  tcase_add_test(tc_core, test_setbits);
  tcase_add_test(tc_core, test_bitshl);
  tcase_add_test(tc_core, test_bitshl2);
  tcase_add_test(tc_core, test_bitcopy);
  tcase_add_test(tc_core, test_count_bits_x);
  suite_add_tcase(s, tc_core);

  return s;
}
