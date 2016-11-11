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
#include <libswiftnav/logging.h>
#include <libswiftnav/track.h>
#include <libswiftnav/constants.h>

/* real ephemeris for SV 15 at 8:00, 30-May-2016 */
const ephemeris_t e_in = {
    .sid = {
      .sat = 15,
      .code = 0
    },
    .toe = {
      .tow = 115200.0,
      .wn  = 1899
    },
    .ura = 1.0,
    .fit_interval = 14400,
    .valid = 1,
    .health_bits = 0,
    .kepler = {
      .tgd       =        -0.00000001071020960808,
      .crc       =       222.34375000000000000000,
      .crs       =        62.87500000000000000000,
      .cuc       =         0.00000356882810592651,
      .cus       =         0.00000734068453311920,
      .cic       =         0.00000002793967723846,
      .cis       =         0.00000012852251529694,
      .dn        =         0.00000000562166273625,
      .m0        =         1.32453305390110598339,
      .ecc       =         0.00812081422191113234,
      .sqrta     =      5153.59741973876953125000,
      .omega0    =        -0.25515026323261330576,
      .omegadot  =        -0.00000000872250618454,
      .w         =         0.48344690815566465636,
      .inc       =         0.93105119838528094256,
      .inc_dot   =         0.00000000008536069847,
      .af0       =        -0.00031310319900512695,
      .af1       =        -0.00000000000170530257,
      .af2       =         0.00000000000000000000,
      .toc = {
        .tow = 115200.0,
        .wn  = 1899
      },
      .iodc      = 40,
      .iode      = 40
    }
};

/* Real measurements from tracking (Piksi v3 board) channel for L1C/A, SV 15. */
const channel_measurement_t l1ca_meas_in= {
    .sid = {
      .sat = 15,
      .code = 0
    },
    .code_phase_chips =      1022.89805441629141569138,
    .code_phase_rate =   1023002.99259400367736816406,
    .carrier_phase = -116700123.04239280521869659424,
    .carrier_freq =      4405.08837890625000000000,
    .time_of_week_ms = 112250800,
    .rec_time_delta =        -0.00008995220125786164,
    .snr =        51.37894058227539062500,
    .lock_counter = 47858
};

/* Real measurements from tracking (Piksi v3 board) channel for L2CM, SV 15. */
const channel_measurement_t l2cm_meas_in= {
    .sid = {
      .sat = 15,
      .code = 1
    },
    .code_phase_chips =     20459.99452955205924808979,
    .code_phase_rate =   1023002.87969136238098144531,
    .carrier_phase = -67041246.17858920246362686157,
    .carrier_freq =      3433.61230468750000000000,
    .time_of_week_ms = 112250800,
    .rec_time_delta =        -0.00039568301886792450,
    .snr =        38.63439178466796875000,
    .lock_counter = 13049
};

START_TEST(test_l2c_meas_data)
{
  const channel_measurement_t *input_meas[1];
  navigation_measurement_t out_l1ca, out_l2cm;
  navigation_measurement_t *output_meas_l1ca[1] = { &out_l1ca };
  navigation_measurement_t *output_meas_l2c[1] = { &out_l2cm };
  const ephemeris_t* e[1] = { &e_in }; /* Use one ephemeris
                                        for both calculations */
  input_meas[0] = &l1ca_meas_in;
  calc_navigation_measurement(1, input_meas, output_meas_l1ca, 0, e);
  log_debug(" ***** L1CA: *****\n");
  log_debug("raw_pseudorange = %30.20f\n", out_l1ca.raw_pseudorange);
  log_debug("pseudorange = %30.20f\n", out_l1ca.pseudorange);
  log_debug("carrier_phase = %30.20f\n", out_l1ca.carrier_phase);
  log_debug("raw_doppler = %30.20f\n", out_l1ca.raw_doppler);
  log_debug("doppler = %30.20f\n", out_l1ca.doppler);
  log_debug("sat_pos = %30.20f, %30.20f, %30.20f\n", out_l1ca.sat_pos[0],out_l1ca.sat_pos[1],out_l1ca.sat_pos[2]);
  log_debug("sat_vel = %30.20f, %30.20f, %30.20f\n", out_l1ca.sat_vel[0],out_l1ca.sat_vel[1],out_l1ca.sat_vel[2]);
  log_debug("snr = %30.20f\n", out_l1ca.snr);
  log_debug("lock_time = %30.20f\n", out_l1ca.lock_time);
  log_debug("tow = %30.20f, wn = %d\n", out_l1ca.tot.tow, out_l1ca.tot.wn);
  log_debug("sat = %u, code = %u\n", (unsigned int)out_l1ca.sid.sat, (unsigned int)out_l1ca.sid.code);
  log_debug("lock_counter = %u\n", (unsigned int)out_l1ca.lock_counter);

  input_meas[0] = &l2cm_meas_in;
  calc_navigation_measurement(1, input_meas, output_meas_l2c, 0, e);
  log_debug(" \n***** L2CM: *****\n");
  log_debug("raw_pseudorange = %30.20f\n", out_l2cm.raw_pseudorange);
  log_debug("pseudorange = %30.20f\n", out_l2cm.pseudorange);
  log_debug("carrier_phase = %30.20f\n", out_l2cm.carrier_phase);
  log_debug("raw_doppler = %30.20f\n", out_l2cm.raw_doppler);
  log_debug("doppler = %30.20f\n", out_l2cm.doppler);
  log_debug("sat_pos = %30.20f, %30.20f, %30.20f\n", out_l2cm.sat_pos[0],out_l2cm.sat_pos[1],out_l2cm.sat_pos[2]);
  log_debug("sat_vel = %30.20f, %30.20f, %30.20f\n", out_l2cm.sat_vel[0],out_l2cm.sat_vel[1],out_l2cm.sat_vel[2]);
  log_debug("snr = %30.20f\n", out_l2cm.snr);
  log_debug("lock_time = %30.20f\n", out_l2cm.lock_time);
  log_debug("tow = %30.20f, wn = %d\n", out_l2cm.tot.tow, out_l2cm.tot.wn);
  log_debug("sat = %u, code = %u\n", (unsigned int)out_l2cm.sid.sat, (unsigned int)out_l2cm.sid.code);
  log_debug("lock_counter = %u\n", (unsigned int)out_l2cm.lock_counter);

  double check_value = fabs(sqrt(pow(out_l1ca.sat_pos[0],2)  +
                                 pow(out_l1ca.sat_pos[1],2)  +
                                 pow(out_l1ca.sat_pos[2],2)) -
                            sqrt(pow(out_l2cm.sat_pos[0],2)  +
                                 pow(out_l2cm.sat_pos[1],2)  +
                                 pow(out_l2cm.sat_pos[2],2)));
  fail_unless(check_value < 0.02,
              "Position difference out of range: %10.8f m, expected < 0.02\n",
              check_value);

  check_value = fabs(out_l1ca.pseudorange - out_l2cm.pseudorange);
  fail_unless(check_value < 1e-6,
              "Pseudorange difference out of range: %15.13f m, expected < 1e-6\n",
              check_value);

  check_value = fabs(out_l1ca.doppler/out_l2cm.doppler - GPS_L1_HZ/GPS_L2_HZ);
  fail_unless(check_value < 0.002,
              "Doppler rate out of range: %15.10f, expected < 0.002\n",
              check_value);
}
END_TEST

Suite* l2c_meas_test_suite(void)
{
  Suite *s = suite_create("L2C meas data");
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_l2c_meas_data);
  suite_add_tcase(s, tc_core);

  return s;
}
