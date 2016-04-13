/*
 * Copyright (C) 2016 Swift Navigation Inc.
 * Contact: Dmitry Tatarinov <dmitry.tatarinov@exafore.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */
//#define DEBUG 1

#include <stdio.h>
#include <check.h>
#include <math.h>
#include <string.h>
#include <libswiftnav/nav_msg_glo.h>
#include <libswiftnav/logging.h>

nav_msg_glo_t n;
ephemeris_t e;
/*This input strings were taken from collected data,
 * refer to libswiftnav/tests/data/gloframesstream/raw_glo_frame_ascii.log
 * #GLORAWFRAMEA,USB1,0,73.0,SATTIME,1892,300827.000,00000000,8792,13498;
 * 3,55,4,1892,300827.077,52,52,15,
 * 01074396999b05c3a850b5,0,021760a5256204d9c15f66,0,0380269d60899a6d0e3123,0,
 * 04865d1cc0000000344918,0,050d100000000340000895,0,06ab81d85b1019f107b83c,0,
 * 0705d3179ea697fc554079,0,082c00148c06153f8e133e,0,0972222bcd4e97ff14be12,0,
 * 0aad8090a54019cb035d3f,0,0b3e2240201e97fc34fc39,0,0cae940cdc3c1e2786da9b,0,
 * 0d68bf54a4c697f115320b,0,0eaf8449b38c1e228932d8,0,0f815b653eee981314802a,0*f3837a1c
 */
u32 strings_in[6][3] ={
  { 1, 1, 1 }, /* dummy words used in test_nav_msg_update_glo only */
  { 0xc3a850b5, 0x96999b05, 0x010743 }, /* 01074396999b05c3a850b5 */
  { 0xd9c15f66, 0xa5256204, 0x021760 }, /* 021760a5256204d9c15f66 */
  { 0x6d0e3123, 0x9d60899a, 0x038026 }, /* 0380269d60899a6d0e3123 */
  { 0x00344918, 0x1cc00000, 0x04865d }, /* 04865d1cc0000000344918 */
  { 0x40000895, 0x3, 0x050d10 }         /* 050d100000000340000895 */
};

/* RAW strings above correspond to following data from same file:
 * #GLOEPHEMERISA,USB1,0,73.5,SATTIME,1892,300617.000,00000000,8d29,13498;
 * 55,4,1,4,1892,301517000,10783,104,0,0,59,0,
 * -1.4453039062500000e+07,-6.9681713867187500e+06,1.9873773925781250e+07, <-- X, Y, Z
 * -1.4125013351440430e+03,-2.3216266632080078e+03,-1.8360681533813477e+03, <-- Vx, Vy, Vz
 * 0.00000000000000000,0.00000000000000000,-2.79396772384643555e-06, <-- Ax, Ay, Az
 * -9.71024855971336365e-05, <-- tau
 * 5.587935448e-09,
 * 1.81898940354585648e-12, <-- gamma
 * 52200,3,0,0,13*955c64e9 
 */
double X = -1.4453039062500000e+07;
double Y = -6.9681713867187500e+06;
double Z = 1.9873773925781250e+07;
double VX = -1.4125013351440430e+03;
double VY = -2.3216266632080078e+03;
double VZ = -1.8360681533813477e+03;
double AX = 0;
double AY = 0;
double AZ = -2.79396772384643555e-06;
double GAMMA = 1.81898940354585648e-12;
double TAU = -9.71024855971336365e-05;

void e_out(void)
{
  log_debug("GLO Ephemeris:\n");
  log_debug("\tSID: %u (code %u)\n", e.sid.sat, e.sid.code);
  log_debug("\tGPS time: TOE %f, WN %d\n", e.toe.tow, e.toe.wn);
  log_debug("\tURA: %f\n", e.ura);
  log_debug("\tFit interval: %u\n", e.fit_interval);
  log_debug("\tValid: %u\n", e.valid);
  log_debug("\tHealthy: %u\n", e.healthy);
  log_debug("\tgamma: %25.18f\n", e.glo.gamma);
  log_debug("\ttau: %25.18f\n", e.glo.tau);
  log_debug("\tX, Y, Z: %25.18f, %25.18f, %25.18f\n", e.glo.pos[0],
            e.glo.pos[1], e.glo.pos[2]);
  log_debug("\tVX, VY, VZ: %25.18f, %25.18f, %25.18f\n", e.glo.vel[0],
            e.glo.vel[1], e.glo.vel[2]);
  log_debug("\tAX, AY, AZ: %25.18f, %25.18f, %25.18f\n", e.glo.acc[0],
            e.glo.acc[1], e.glo.acc[2]);
  fail_unless(e.glo.pos[0] - X == 0, "dX %25.18f, expected %25.18f",
              e.glo.pos[0], X);
  fail_unless(e.glo.pos[1] - Y == 0, "dY %25.18f, expected %25.18f",
              e.glo.pos[1], Y);
  fail_unless(e.glo.pos[2] - Z == 0, "dZ %25.18f, expected %25.18f",
              e.glo.pos[2], Z);
  fail_unless(e.glo.vel[0] - VX == 0, "dVX %25.18f, expected %25.18f",
              e.glo.vel[0], VX);
  fail_unless(e.glo.vel[1] - VY == 0, "dVY %25.18f, expected %25.18f",
              e.glo.vel[1], VY);
  fail_unless(e.glo.vel[2] - VZ == 0, "dVZ %25.18f, expected %25.18f",
              e.glo.vel[2], VZ);
  fail_unless(e.glo.acc[0] - AX == 0, "dAX %25.18f, expected %25.18f",
              e.glo.acc[0], AX);
  fail_unless(e.glo.acc[1] - AY == 0, "dAY %25.18f, expected %25.18f",
              e.glo.acc[1], AY);
  fail_unless(e.glo.acc[2] - AZ == 0, "dAZ %25.18f, expected %25.18f",
              e.glo.acc[2], AZ);
  fail_unless(e.glo.tau - TAU == 0, "dTAU %25.18f, expected %25.18f", e.glo.tau,
              TAU);
  fail_unless(e.glo.gamma - GAMMA == 0, "dGAMMA %25.18f, expected %25.18f",
              e.glo.gamma, GAMMA);
}

START_TEST(test_extract_glo_word)
{
  u32 ret = 0;
  nav_msg_init_glo(&n);
  n.string_bits[0] = 5;
  n.string_bits[1] = 5;
  n.string_bits[2] = 5;
  ret = extract_word_glo(&n, 1, 32);
  fail_unless(ret == 5, "1. %x, expected %x", ret, 5);
  ret = extract_word_glo(&n, 33, 3);
  fail_unless(ret == 5, "2. %x, expected %x", ret, 5);
  ret = extract_word_glo(&n, 65, 3);
  fail_unless(ret == 5, "3. %x, expected %x", ret, 5);

  n.string_bits[0] = 0x12345678;
  n.string_bits[1] = 0xdeadbeef;
  n.string_bits[2] = 0x87654321;
  ret = extract_word_glo(&n, 1, 32);
  fail_unless(ret == 0x12345678, "4. %x, expected %x", ret, 0x12345678);
  ret = extract_word_glo(&n, 33, 32);
  fail_unless(ret == 0xdeadbeef, "5. %x, expected %x", ret, 0xdeadbeef);
  ret = extract_word_glo(&n, 65, 32);
  fail_unless(ret == 0x87654321, "6. %x, expected %x", ret, 0x87654321);
  ret = extract_word_glo(&n, 49, 4);
  fail_unless(ret == 0xd, "7. %x, expected %x", ret, 0xd);

  n.string_bits[0] = 0xbeef0000;
  n.string_bits[1] = 0x4321dead;
  n.string_bits[2] = 0x00008765;
  ret = extract_word_glo(&n, 17, 32);
  fail_unless(ret == 0xdeadbeef, "8. %x, expected %x", ret, 0xdeadbeef);
  ret = extract_word_glo(&n, 49, 32);
  fail_unless(ret == 0x87654321, "9. %x, expected %x", ret, 0x87654321);
  ret = extract_word_glo(&n, 49, 16);
  fail_unless(ret == 0x4321, "10. %x, expected %x", ret, 0x4321);
}
END_TEST

START_TEST(test_process_string_glo)
{
  nav_msg_init_glo(&n);
  memset(&e, 0, sizeof(e));
  for (u8 i = 1; i < sizeof(strings_in) / sizeof(strings_in[1]); i++) {
    memcpy(n.string_bits, strings_in[i], sizeof(n.string_bits));
    process_string_glo(&n, &e);
  }
  e_out();
}
END_TEST

START_TEST(test_nav_msg_update_glo)
{
  /* the unit test encodes strings_in to generate glo bitstream, calls
   * nav_msg_update_glo to receive and finally decodes received string */
  nav_msg_init_glo(&n);
  memset(&e, 0, sizeof(e));
  /* get string one by one */
  for (u8 i = 0; i < sizeof(strings_in) / sizeof(strings_in[0]); i++) {
    u8 manchester = 0;
    nav_msg_glo_t a;
    s8 ret;
    u8 j;
    nav_msg_init_glo(&a);
    /* write test string to temporary buffer */
    memcpy(a.string_bits, strings_in[i], sizeof(n.string_bits));
    /* transmit data bits, 85 bit */
    for (j = 85; j > 0; j--) {
      bool one_bit = extract_word_glo(&a, j, 1); /* get bit to be transmitted */
      manchester = (one_bit << 1 | one_bit) ^ 2; /* transform to line code */
      /* now pass it to receiver MSB first, receiver must return -1 */
      ret = nav_msg_update_glo(&n, (manchester >> 1) & 1);
      fail_unless(ret == -1, "ret = %d, expected -1", ret);
      /* now LSB, receiver must return -1 */
      ret = nav_msg_update_glo(&n, manchester & 1);
    }
    /* try to decode the string */
    if (ret == 1) {
      fail_unless(
          memcmp(a.string_bits, n.string_bits, sizeof(a.string_bits)) == 0,
          "Received string %x%x%x not equal to trasmitted one %x%x%x",
          n.string_bits[2], n.string_bits[1], n.string_bits[0],
          a.string_bits[2], a.string_bits[1], a.string_bits[0]);

      if (process_string_glo(&n, &e) == 1)
        e_out();
    }
    /* now pass time mark bit by bit to receiver (MSB first),
     * no line code needed */
    for (u8 j = 30; j > 0; j--) {
      ret = nav_msg_update_glo(&n, (GLO_TM >> (j - 1)) & 1);
      fail_unless(ret == -1, "ret = %d, expected -1", ret);
    }
  }
}
END_TEST

Suite* glo_decoder_test_suite(void)
{
  Suite *s = suite_create("GLO decoder");
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_extract_glo_word);
  tcase_add_test(tc_core, test_process_string_glo);
  tcase_add_test(tc_core, test_nav_msg_update_glo);
  suite_add_tcase(s, tc_core);

  return s;
}
