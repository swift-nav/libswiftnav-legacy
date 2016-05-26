#include <libswiftnav/cnav_msg.h>
#include <libswiftnav/edc.h>
#include <libswiftnav/bits.h>

#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <check.h>

#include "check_utils.h"

#include <cnav_msg.c>

#define G1 0x79 /* convolution encoder coefficient A 0171 */
#define G2 0x5B /* convolution encoder coefficient B 0133 */
#define SUFFIX_SIZE 16 /* Data to feed after the last message to ensure the
                        * decoder completes message. */

static inline size_t encode_byte(u32 *pacc, u8 v, u8 *dst)
{
  u32 acc = *pacc;
  for (int n_bit = CHAR_BIT - 1; n_bit >= 0; --n_bit) {
    u32 t = ((v >> n_bit) & 1) << 6;
    acc = (acc >> 1) | t;
    *dst++ = parity(acc & G1) ? 0xFF : 00;
    *dst++ = parity(acc & G2) ? 0xFF : 00;
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
      dst[cnt++] = parity(acc & G1) ? 0xFF : 00;
      dst[cnt++] = parity(acc & G2) ? 0xFF : 00;
    }
    *pacc = acc;
  }

  return cnt;
}

START_TEST(test_cnav_compute_crc)
{
  cnav_v27_part_t part = { .n_decoded = GPS_CNAV_MSG_LENGTH };
  u32 crc = _cnav_compute_crc(&part);
  fail_unless(0 == crc, "CRC: 0x000000 0x%06" PRIX32, crc);

  part.decoded[0] = 1;
  crc = _cnav_compute_crc(&part);
  fail_unless(0x733518 == crc, "CRC: 0x733518 0x%06" PRIX32, crc);
}
END_TEST

START_TEST(test_cnav_extract_crc)
{
  cnav_v27_part_t part = { .n_decoded = GPS_CNAV_MSG_LENGTH };
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
  cnav_v27_part_t part = { .n_decoded = GPS_CNAV_MSG_LENGTH };

  _cnav_rescan_preamble(&part);
  fail_unless(7 == part.n_decoded, "Decoded: 7 != %zu", part.n_decoded);
  fail_unless(0 == part.decoded[0],
              "Decoded[0]: 0x00 != 0x%02" PRIX8, (u8)part.decoded[0]);

  part.n_decoded = GPS_CNAV_MSG_LENGTH;
  part.decoded[36] = 0x07;
  part.decoded[37] = 0xF0;
  _cnav_rescan_preamble(&part);
  fail_unless(7 == part.n_decoded, "Decoded: 7 != %zu", part.n_decoded);
  fail_unless(0xFE == part.decoded[0],
              "Decoded[0]: 0xFE != 0x%02" PRIX8, (u8)part.decoded[0]);

  part.n_decoded = GPS_CNAV_MSG_LENGTH;
  part.decoded[10] = GPS_CNAV_PREAMBLE1;
  part.decoded[36] = 0x07;
  part.decoded[37] = 0xF0;
  _cnav_rescan_preamble(&part);
  fail_unless(220 == part.n_decoded, "Decoded: 220 != %zu", part.n_decoded);
  fail_unless(GPS_CNAV_PREAMBLE1 == part.decoded[0],
              "Decoded[0]: GPS_CNAV_PREAMBLE1 != 0x%02" PRIX8, (u8)part.decoded[0]);
  fail_unless(false == part.invert, "Inverted: false != true");

  part.n_decoded = GPS_CNAV_MSG_LENGTH;
  part.decoded[10] = GPS_CNAV_PREAMBLE2;
  part.decoded[36] = 0x07;
  part.decoded[37] = 0xF0;
  _cnav_rescan_preamble(&part);
  fail_unless(220 == part.n_decoded, "Decoded: 220 != %zu", part.n_decoded);
  fail_unless(GPS_CNAV_PREAMBLE2 == part.decoded[0],
              "Decoded[0]: GPS_CNAV_PREAMBLE2 != 0x%02" PRIX8, (u8)part.decoded[0]);
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

static size_t add_encoded_message(u32 *pacc, u8 *buf)
{
  u8 src_msg[38] = {0};
  memset(src_msg, 0x55, sizeof(src_msg));
  setbitu(src_msg, 0,  8, GPS_CNAV_PREAMBLE1); /* preamble */
  setbitu(src_msg, 8,  6, 22);          /* prn */
  setbitu(src_msg, 14, 6, 0);           /* msg */
  setbitu(src_msg, 20, 17, 1000);       /* tow */
  setbitu(src_msg, 37, 1, 0);           /* flag */

  u32 crc = crc24q_bits(0, src_msg, GPS_CNAV_MSG_DATA_LENGTH, false);
  setbitu(src_msg, GPS_CNAV_MSG_DATA_LENGTH, GPS_CNAV_MSG_CRC_LENGTH, crc);

  return encode_bits(pacc, src_msg, GPS_CNAV_MSG_LENGTH, buf);
}

static size_t encode_buffer_p(u8 *buf, const u8 *prefix, size_t prefix_size)
{
  u8 suffix[SUFFIX_SIZE] = {0};
  memset(suffix, 0x55, sizeof(suffix));

  size_t dst = 0;
  u32 acc = 0;
  dst += encode(&acc, prefix, prefix_size, buf + dst);
  dst += add_encoded_message(&acc, buf + dst);
  dst += encode(&acc, suffix, sizeof(suffix), buf + dst);

  return dst;
}

static size_t encode_buffer(u8 *buf)
{
  u8 prefix[4] = {0};
  memset(prefix, 0x55, sizeof(prefix));

  return encode_buffer_p(buf, prefix, sizeof(prefix));
}

START_TEST(test_cnav_decode1)
{
  cnav_msg_decoder_t dec;
  cnav_msg_decoder_init(&dec);
  cnav_msg_t msg;
  size_t i;

  u8 enc[((4 + SUFFIX_SIZE) * CHAR_BIT + GPS_CNAV_MSG_LENGTH) * 2];
  size_t dst = encode_buffer(enc);

  fail_if((4 + SUFFIX_SIZE) * CHAR_BIT * 2 + GPS_CNAV_MSG_LENGTH * 2 != dst,
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

  fail_if(!ret, "Message has not been decoded");
}
END_TEST

START_TEST(test_cnav_decode2)
{
  cnav_msg_decoder_t dec;
  cnav_msg_decoder_init(&dec);
  cnav_msg_t msg;
  size_t i;

  u8 enc[((4 + SUFFIX_SIZE) * CHAR_BIT + GPS_CNAV_MSG_LENGTH) * 2 + 1];
  size_t dst = encode_buffer(enc + 1);

  fail_if((4 + SUFFIX_SIZE) * CHAR_BIT * 2 + GPS_CNAV_MSG_LENGTH * 2 != dst,
          "Buffer size mismatch: expected=%zu actual=%zu",
          sizeof(enc), dst);

  bool ret = false;
  for (i = 0; i < dst; ++i) {
    u32 delay = 0;

    ret = cnav_msg_decoder_add_symbol(&dec, enc[i], &msg, &delay);

    if (ret) {

      u32 expected_delay;
      fail_if (msg.prn    != 22,    "Wrong PRN");
      fail_if (msg.msg_id != 0,     "Wrong MSG");
      fail_if (msg.tow    != 1000,  "Wrong TOW");
      fail_if (msg.alert  != false, "Wrong alert");
      expected_delay = i + 1; /* symbols fed into the decoder */
      expected_delay -= 600;  /* symbols in L2C message */
      expected_delay -= 64;   /* size of dummy header outside
                                 of data message [symbols] */
      fail_if (expected_delay == delay);

      break;
    }
  }

  fail_if(!ret, "Message has not been decoded");
}
END_TEST

START_TEST(test_cnav_decode_neg)
{
  cnav_msg_decoder_t dec;
  cnav_msg_decoder_init(&dec);
  cnav_msg_t msg;
  size_t i;

  u8 enc[((4 + SUFFIX_SIZE) * CHAR_BIT + GPS_CNAV_MSG_LENGTH) * 2];
  size_t dst = encode_buffer(enc);

  fail_if((4 + SUFFIX_SIZE) * CHAR_BIT * 2 + GPS_CNAV_MSG_LENGTH * 2 != dst,
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

  fail_if(!ret, "Message has not been decoded");
}
END_TEST


START_TEST(test_cnav_decode_false)
{
  cnav_msg_decoder_t dec;
  cnav_msg_decoder_init(&dec);
  cnav_msg_t msg;
  size_t i;

  u8 prefix[] = {0x55, 0x55, 0x55, 0x55, 0x8B, 0x8B};
  u8 suffix[SUFFIX_SIZE] = {0};
  u8 enc[((sizeof(prefix) + sizeof(suffix)) * CHAR_BIT + GPS_CNAV_MSG_LENGTH) * 2 + 1];
  size_t dst = encode_buffer_p(enc, prefix, sizeof(prefix));

  fail_if(((sizeof(prefix) + sizeof(suffix)) * CHAR_BIT + GPS_CNAV_MSG_LENGTH) * 2 != dst,
          "Buffer size mismatch: expected=%zu actual=%zu",
          sizeof(enc), dst);

  u32 counter = 0;
  for (i = 0; i < dst; ++i) {
    u32 delay = 0;

    bool ret = cnav_msg_decoder_add_symbol(&dec, enc[i], &msg, &delay);

    if (ret) {
      counter ++;

      fail_if (msg.prn    != 22,    "Wrong PRN");
      fail_if (msg.msg_id != 0,     "Wrong MSG");
      fail_if (msg.tow    != 1000,  "Wrong TOW");
      fail_if (msg.alert  != false, "Wrong alert");
    }
  }

  fail_if(counter == 0, "Message has not been decoded");
}
END_TEST

START_TEST(test_cnav_decode_unlock)
{
  cnav_msg_decoder_t dec;
  cnav_msg_decoder_init(&dec);
  cnav_msg_t msg;
  size_t i;


  u8 prefix[4];
  u8 suffix[SUFFIX_SIZE];

  memset(prefix, 0x55, sizeof(prefix));
  memset(suffix, 0x55, sizeof(suffix));

  u8 enc[((sizeof(prefix) + sizeof(suffix)) * CHAR_BIT + /* Prefix + suffix */
         GPS_CNAV_MSG_LENGTH + /* Message CRC_OK */
         GPS_CNAV_MSG_LENGTH * (1 + GPS_CNAV_LOCK_MAX_CRC_FAILS) + /* Messages CRC bad */
         GPS_CNAV_MSG_LENGTH /* Message CRC_OK */
         ) * 2 + 1];
  u8 bad_message[38] = {0};
  memset(bad_message, 0x55, sizeof(bad_message));

  size_t dst = 0;
  u32 acc = 0;

  dst += encode(&acc, prefix, sizeof(prefix), &enc[0] + dst);
  dst += add_encoded_message(&acc, &0[enc] + dst);

  for (i = 0; i <= GPS_CNAV_LOCK_MAX_CRC_FAILS; ++i) {
    dst += encode_bits(&acc, bad_message, GPS_CNAV_MSG_LENGTH, &enc[0] + dst);
  }

  dst += add_encoded_message(&acc, &0[enc] + dst);
  dst += encode(&acc, suffix, sizeof(suffix), &enc[0] + dst);

  fail_if(((sizeof(prefix) + sizeof(suffix)) * CHAR_BIT +
           GPS_CNAV_MSG_LENGTH * (GPS_CNAV_LOCK_MAX_CRC_FAILS + 3)) * 2 != dst,
          "Buffer size mismatch: expected=%zu actual=%zu",
          sizeof(enc), dst);


  u32 decoded = 0;
  for (i = 0; i < dst; ++i) {
    u32 delay = 0;

    bool ret = cnav_msg_decoder_add_symbol(&dec, enc[i], &msg, &delay);

    if (ret) {
      decoded++;

      fail_if (msg.prn    != 22,    "Wrong PRN: %" PRIu8, msg.prn);
      fail_if (msg.msg_id != 0,     "Wrong MSG: %" PRIu8, msg.msg_id);
      fail_if (msg.tow    != 1000,  "Wrong TOW: %" PRIu32, msg.tow);
      fail_if (msg.alert  != false, "Wrong alert: %d", msg.alert);
    }
  }

  fail_if(decoded != 2, "Message has not been decoded twice; dec=%"PRIu32, decoded);
}
END_TEST

START_TEST(test_cnav_decode_ok_nok_ok)
{
  cnav_msg_decoder_t dec;
  cnav_msg_decoder_init(&dec);
  cnav_msg_t msg;
  size_t i;


  u8 prefix[4];
  u8 suffix[SUFFIX_SIZE];

  memset(prefix, 0x55, sizeof(prefix));
  memset(suffix, 0x55, sizeof(suffix));

  u8 enc[((sizeof(prefix) + sizeof(suffix)) * CHAR_BIT + /* Prefix + suffix */
         GPS_CNAV_MSG_LENGTH + /* Message CRC_OK */
         GPS_CNAV_MSG_LENGTH + /* Messages CRC bad */
         GPS_CNAV_MSG_LENGTH /* Message CRC_OK */
         ) * 2 + 1];
  u8 bad_message[38] = {0};
  memset(bad_message, 0x55, sizeof(bad_message));

  size_t dst = 0;
  u32 acc = 0;

  dst += encode(&acc, prefix, sizeof(prefix), &enc[0] + dst);
  dst += add_encoded_message(&acc, &0[enc] + dst);

  dst += encode_bits(&acc, bad_message, GPS_CNAV_MSG_LENGTH, &enc[0] + dst);

  dst += add_encoded_message(&acc, &0[enc] + dst);
  dst += encode(&acc, suffix, sizeof(suffix), &enc[0] + dst);

  fail_if(((sizeof(prefix) + sizeof(suffix)) * CHAR_BIT +
           GPS_CNAV_MSG_LENGTH * 3) * 2 != dst,
          "Buffer size mismatch: expected=%zu actual=%zu",
          sizeof(enc), dst);


  u32 decoded = 0;
  for (i = 0; i < dst; ++i) {
    u32 delay = 0;

    bool ret = cnav_msg_decoder_add_symbol(&dec, enc[i], &msg, &delay);

    if (ret) {
      decoded++;

      fail_if (msg.prn    != 22,    "Wrong PRN: %" PRIu8, msg.prn);
      fail_if (msg.msg_id != 0,     "Wrong MSG: %" PRIu8, msg.msg_id);
      fail_if (msg.tow    != 1000,  "Wrong TOW: %" PRIu32, msg.tow);
      fail_if (msg.alert  != false, "Wrong alert: %d", msg.alert);
    }
  }

  fail_if(decoded != 2, "Message has not been decoded twice; dec=%"PRIu32, decoded);
}
END_TEST


START_TEST(test_cnav_crc)
{
  u32 crc0, crc1;
  u8 data0[35] = {0x0D, 0xEA, 0xDB, 0xEE, 0xF0};

  crc0 = crc24q(data0, 35, 0);
  bitshl(data0, 35, 4);
  crc1 = crc24q_bits(0, data0, 35 * 8 - 4, false);

  fail_if(crc0 != crc1, "CRC mismatch");
}
END_TEST

Suite* cnav_test_suite(void)
{
  Suite *s = suite_create("LNAV");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_cnav_crc);
  tcase_add_test(tc_core, test_cnav_compute_crc);
  tcase_add_test(tc_core, test_cnav_extract_crc);
  tcase_add_test(tc_core, test_cnav_rescan_preamble);
  tcase_add_test(tc_core, test_cnav_decoder_init);
  tcase_add_test(tc_core, test_cnav_decode1);
  tcase_add_test(tc_core, test_cnav_decode_neg);
  tcase_add_test(tc_core, test_cnav_decode2);
  tcase_add_test(tc_core, test_cnav_decode_false);
  tcase_add_test(tc_core, test_cnav_decode_unlock);
  tcase_add_test(tc_core, test_cnav_decode_ok_nok_ok);

  suite_add_tcase(s, tc_core);

  return s;
}
