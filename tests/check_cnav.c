#include <libswiftnav/cnav_msg.h>
#include <libswiftnav/edc.h>
#include <libswiftnav/bits.h>

#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <check.h>

#include "check_utils.h"

#define G1 0x79 // 0171
#define G2 0x5B // 0133

// internal CNAV declarations for test access
void _cnav_bitlshift(void *buf, size_t bits, u32 shift);
u32 _cnav_compute_crc(cnav_v27_part_t *part);
u32 _cnav_extract_crc(const cnav_v27_part_t *part);
void _cnav_rescan_preamble(cnav_v27_part_t *part);
void _cnav_add_symbol(cnav_v27_part_t *part, u8 ch);
void _cnav_msg_invert(cnav_v27_part_t *part);
bool _cnav_msg_decode(cnav_v27_part_t *part, cnav_msg_t *msg, u32 *delay);

static inline u32 parity8(u32 v)
{
  return (0x6996 >> ((v ^ (v >> 4)) & 15)) & 1;
}
static inline size_t encode_byte(u32 *pacc, u8 v, u8 *dst)
{
  u32 acc = *pacc;
  for (int n_bit = CHAR_BIT - 1; n_bit >= 0; --n_bit) {
    u32 t = ((v >> n_bit) & 1) << 6;
    acc = (acc >> 1) | t;
    *dst++ = parity8(acc & G1) ? 0xFF : 0x00;
    *dst++ = parity8(acc & G2) ? 0xFF : 0x00;
  }
  *pacc = acc;
  return CHAR_BIT * 2;
}
static inline size_t encode(u32 *pacc, const u8 *src, size_t bytes, u8 *dst)
{
  size_t cnt = 0;
  for (size_t i = 0; i < bytes; ++i) {
    cnt += encode_byte(pacc, src[i], dst + cnt);
  }
  return cnt;
}
static inline size_t encode_bits(u32 *pacc, const u8 *src, size_t bits, u8 *dst)
{
  size_t cnt = 0;
  size_t n_bytes = bits / CHAR_BIT;

  for (size_t i = 0; i < n_bytes; ++i) {
    cnt += encode_byte(pacc, src[i], dst + cnt);
  }
  if (bits - n_bytes * CHAR_BIT != 0) {
    size_t x_bits = bits % CHAR_BIT;
    u8 v = src[n_bytes];
    u32 acc = *pacc;
    for (size_t i = 0; i < x_bits; ++i) {
      u32 t = ((v >> (CHAR_BIT - 1 - i)) & 1) << 6;
      acc = (acc >> 1) | t;
      *dst++ = parity8(acc & G1) ? 0xFF : 0x00;
      *dst++ = parity8(acc & G2) ? 0xFF : 0x00;
      cnt += 2;
    }
    *pacc = acc;
  }

  return cnt;
}

START_TEST(test_cnav_bitlshift)
{
  u8 src0[] = { 0xDE, 0xAD, 0xBE, 0xEF };
  u8 res0[] = { 0xBE, 0xEF, 0x00, 0x00 };

  u8 src1[] = { 0xDE, 0xAD, 0xBE, 0xEF };
  u8 res1[] = { 0xEA, 0xDB, 0xEE, 0xF0 };

  u8 src2[] = { 0xDE, 0xAD, 0xBE, 0xEF };
  u8 res2[] = { 0xDB, 0xEE, 0xF0, 0x00 };

  u8 src3[] = { 0xDE, 0xAD, 0xBE, 0xEF };
  u8 res3[] = { 0xB6, 0xFB, 0xBC, 0x00 };

  _cnav_bitlshift(src0, sizeof(src0), 16);
  fail_unless(0 == memcmp(src0, res0, 4), "Byte shift test");

  _cnav_bitlshift(src1, sizeof(src1), 4);
  fail_unless(0 == memcmp(src1, res1, 4), "4-bit shift");

  _cnav_bitlshift(src2, sizeof(src2), 12);
  fail_unless(0 == memcmp(src2, res2, 4), "12-bit shift");

  _cnav_bitlshift(src3, sizeof(src3), 10);
  fail_unless(0 == memcmp(src3, res3, 4), "10-bit shift");
}
END_TEST

START_TEST(test_cnav_compute_crc)
{
  cnav_v27_part_t part = { .n_decoded = 300 };
  u32 crc = _cnav_compute_crc(&part);
  fail_unless(0 == crc, "CRC: 0x000000 0x%06" PRIX32, crc);

  part.decoded[0] = 1;
  crc = _cnav_compute_crc(&part);
  fail_unless(0x733518 == crc, "CRC: 0x733518 0x%06" PRIX32, crc);
}
END_TEST

START_TEST(test_cnav_extract_crc)
{
  cnav_v27_part_t part = { .n_decoded = 300 };
  u32 crc = _cnav_extract_crc(&part);
  fail_unless(0 == crc, "CRC: 0x000000 0x%06" PRIX32, crc);

  part.decoded[34] = 0xDE;
  part.decoded[35] = 0xAD;
  part.decoded[36] = 0xBE;
  part.decoded[37] = 0xEF;
  crc = _cnav_extract_crc(&part);
  fail_unless(0xEADBEE == crc, "CRC: 0xEADBEE 0x%06" PRIX32, crc);
}
END_TEST

START_TEST(test_cnav_rescan_preamble)
{
  cnav_v27_part_t part = { .n_decoded = 300 };

  _cnav_rescan_preamble(&part);
  fail_unless(7 == part.n_decoded, "Decoded: 7 != %zu", part.n_decoded);
  fail_unless(0 == part.decoded[0], "Decoded[0]: 0x00 != 0x%02" PRIX8, (u8)part.decoded[0]);

  part.n_decoded = 300;
  part.decoded[36] = 0x07;
  part.decoded[37] = 0xF0;
  _cnav_rescan_preamble(&part);
  fail_unless(7 == part.n_decoded, "Decoded: 7 != %zu", part.n_decoded);
  fail_unless(0xFE == part.decoded[0], "Decoded[0]: 0xFE != 0x%02" PRIX8, (u8)part.decoded[0]);

  part.n_decoded = 300;
  part.decoded[10] = 0b10001011u;
  part.decoded[36] = 0x07;
  part.decoded[37] = 0xF0;
  _cnav_rescan_preamble(&part);
  fail_unless(220 == part.n_decoded, "Decoded: 220 != %zu", part.n_decoded);
  fail_unless(0x8B == part.decoded[0], "Decoded[0]: 0x8B != 0x%02" PRIX8, (u8)part.decoded[0]);
  fail_unless(false == part.invert, "Inverted: false != true");

  part.n_decoded = 300;
  part.decoded[10] = 0b01110100u;
  part.decoded[36] = 0x07;
  part.decoded[37] = 0xF0;
  _cnav_rescan_preamble(&part);
  fail_unless(220 == part.n_decoded, "Decoded: 220 != %zu", part.n_decoded);
  fail_unless(0x74 == part.decoded[0], "Decoded[0]: 0x74 != 0x%02" PRIX8, (u8)part.decoded[0]);
  fail_unless(true == part.invert, "Inverted: true != false");
}
END_TEST

START_TEST(test_cnav_decoder_init)
{
  cnav_msg_decoder_t dec;
  cnav_msg_decoder_init(&dec);

  fail_if(dec.part1.invert, "Invert 1 NOK");
  fail_if(dec.part1.message_lock, "Message lock 1 NOK");
  fail_if(dec.part1.preamble_seen, "Preamble 1 NOK");
  fail_if(0 != dec.part1.n_symbols, "Symbols 1 NOK: %zu", dec.part1.n_symbols);
  fail_if(0 != dec.part1.n_decoded, "Decoded 1 NOK: %zu", dec.part1.n_decoded);

  fail_if(dec.part2.invert, "Invert 2 NOK");
  fail_if(dec.part2.message_lock, "Message lock 2 NOK");
  fail_if(dec.part2.preamble_seen, "Preamble 2 NOK");
  fail_if(1 != dec.part2.n_symbols, "Symbols 2 NOK: %zu", dec.part2.n_symbols);
  fail_if(0 != dec.part2.n_decoded, "Decoded 2 NOK: %zu", dec.part2.n_decoded);
}
END_TEST

static size_t encode_buffer_p(u8 *buf, u8 prefix[4])
{
  u8 src_msg[38] = {0};
  memset(src_msg, 0x55, sizeof(src_msg));
  setbitu(src_msg, 0,  8, 0b10001011u); // preamble
  setbitu(src_msg, 8,  6, 22);    // prn
  setbitu(src_msg, 14, 6, 0);     // msg
  setbitu(src_msg, 20, 17, 1000); // tow
  setbitu(src_msg, 37, 1, 0); // x

  u32 crc = crc24q_bits(0, src_msg, 300 - 24, false);
  setbits(src_msg, 276, 24, crc);

  u8 suffix[4+4] = {0};
  memset(suffix, 0x55, sizeof(suffix));

  size_t dst = 0;
  u32 acc = 0;
  dst += encode(&acc, prefix, 4 /* sizeof(prefix)*/, buf + dst);
  dst += encode_bits(&acc, src_msg, 300 /*sizeof(src_msg)*/, buf + dst);
  dst += encode(&acc, suffix, sizeof(suffix), buf + dst);

  return dst;
}

static size_t encode_buffer(u8 *buf)
{
  u8 prefix[4] = {0};
  memset(prefix, 0x55, sizeof(prefix));

  return encode_buffer_p(buf, prefix);
}

START_TEST(test_cnav_decode1)
{
  cnav_msg_decoder_t dec;
  cnav_msg_decoder_init(&dec);
  cnav_msg_t msg;
  size_t i;

  u8 enc[((4 + 8) * CHAR_BIT + 300) * 2];
  size_t dst = encode_buffer(enc);

  fail_if((4 + 8) * CHAR_BIT*2 + 600 != dst,
          "Buffer size mismatch: expected=%zu actual=%zu",
          sizeof(enc), dst);

  bool ret = false;
  for (i = 0; i < dst; ++i) {
    u32 delay = 0;

    ret = cnav_msg_decoder_add_symbol(&dec, enc[i], &msg, &delay);

    if (ret) {

      fail_if (msg.prn    != 22,    "Wrong PRN");
      fail_if (msg.msg_id != 0,     "Wrong MSG");
      fail_if (msg.tow    != 1000,  "Wrong TOW");
      fail_if (msg.alert  != false, "Wrong alert");

      break;
    }
  }

}
END_TEST

START_TEST(test_cnav_decode2)
{
  cnav_msg_decoder_t dec;
  cnav_msg_decoder_init(&dec);
  cnav_msg_t msg;
  size_t i;

  u8 enc[((4 + 8) * CHAR_BIT + 300) * 2 + 1];
  size_t dst = encode_buffer(enc + 1);

  fail_if((4 + 8) * CHAR_BIT*2 + 600 != dst,
          "Buffer size mismatch: expected=%zu actual=%zu",
          sizeof(enc), dst);

  bool ret = false;
  for (i = 0; i < dst; ++i) {
    u32 delay = 0;

    ret = cnav_msg_decoder_add_symbol(&dec, enc[i], &msg, &delay);

    if (ret) {

      fail_if (msg.prn    != 22,    "Wrong PRN");
      fail_if (msg.msg_id != 0,     "Wrong MSG");
      fail_if (msg.tow    != 1000,  "Wrong TOW");
      fail_if (msg.alert  != false, "Wrong alert");

      break;
    }
  }

}
END_TEST

START_TEST(test_cnav_decode_neg)
{
  cnav_msg_decoder_t dec;
  cnav_msg_decoder_init(&dec);
  cnav_msg_t msg;
  size_t i;

  u8 enc[((4 + 8) * CHAR_BIT + 300) * 2];
  size_t dst = encode_buffer(enc);

  fail_if((4 + 8) * CHAR_BIT*2 + 600 != dst,
          "Buffer size mismatch: expected=%zu actual=%zu",
          sizeof(enc), dst);

  for (i = 0; i < sizeof(enc); ++i) {
    enc[i] ^= 0xFFu;
  }

  bool ret = false;

  for (i = 0; i < dst; ++i) {
    u32 delay = 0;
    ret = cnav_msg_decoder_add_symbol(&dec, enc[i], &msg, &delay);

    if (ret) {

      fail_if (msg.prn    != 22,    "Wrong PRN");
      fail_if (msg.msg_id != 0,     "Wrong MSG");
      fail_if (msg.tow    != 1000,  "Wrong TOW");
      fail_if (msg.alert  != false, "Wrong alert");

      break;
    }
  }

}
END_TEST


START_TEST(test_cnav_decode_false)
{
  cnav_msg_decoder_t dec;
  cnav_msg_decoder_init(&dec);
  cnav_msg_t msg;
  size_t i;

  u8 prefix[4] = {0x55, 0x55, 0x8B, 0x8B};
  u8 enc[((4 + 8) * CHAR_BIT + 300) * 2 + 1];
  size_t dst = encode_buffer_p(enc, prefix);

  fail_if((4 + 8) * CHAR_BIT*2 + 600 != dst,
          "Buffer size mismatch: expected=%zu actual=%zu",
          sizeof(enc), dst);

  bool ret = false;
  for (i = 0; i < dst; ++i) {
    u32 delay = 0;

    ret = cnav_msg_decoder_add_symbol(&dec, enc[i], &msg, &delay);

    if (ret) {

      fail_if (msg.prn    != 22,    "Wrong PRN");
      fail_if (msg.msg_id != 0,     "Wrong MSG");
      fail_if (msg.tow    != 1000,  "Wrong TOW");
      fail_if (msg.alert  != false, "Wrong alert");

      break;
    }
  }

}
END_TEST

START_TEST(test_cnav_crc)
{
  u32 crc0, crc1;
  u8 data0[35] = {0x0D, 0xEA, 0xDB, 0xEE, 0xF0};

  crc0 = crc24q(data0, 35, 0);
  _cnav_bitlshift(data0, 35, 4);
  crc1 = crc24q_bits(0, data0, 35 * 8 - 4, false);

  fail_if(crc0 != crc1, "CRC mismatch");
}
END_TEST

Suite* cnav_test_suite(void)
{
  Suite *s = suite_create("LNAV");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_cnav_bitlshift);
  tcase_add_test(tc_core, test_cnav_crc);
  tcase_add_test(tc_core, test_cnav_compute_crc);
  tcase_add_test(tc_core, test_cnav_extract_crc);
  tcase_add_test(tc_core, test_cnav_rescan_preamble);
  tcase_add_test(tc_core, test_cnav_decoder_init);
  tcase_add_test(tc_core, test_cnav_decode1);
  tcase_add_test(tc_core, test_cnav_decode_neg);
  tcase_add_test(tc_core, test_cnav_decode2);
  tcase_add_test(tc_core, test_cnav_decode_false);

  suite_add_tcase(s, tc_core);

  return s;
}
