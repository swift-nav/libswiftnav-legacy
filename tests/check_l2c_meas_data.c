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

/* real ephemeris for SV 5 at 6:00, 17-May-2016 */
const ephemeris_t e_in = {
    .sid = {
      .sat = 5,
      .code = 0
    },
    .toe = {
      .tow = 194400.0,
      .wn  = 1897
    },
    .ura = 2.0,
    .fit_interval = 14400,
    .valid = 1,
    .healthy = 1,
    .kepler = {
      .tgd       = -0.00000001071020960808,
      .crc       = 178.28125000000000000000,
      .crs       = -21.21875000000000000000,
      .cuc       = -0.00000085867941379547,
      .cus       = 0.00001010857522487640,
      .cic       = 0.00000003352761268616,
      .cis       = 0.00000012107193470001,
      .dn        = 0.00000000498092176110,
      .m0        = 1.40317297011784125615,
      .ecc       = 0.00460527464747428894,
      .sqrta     = 5153.79424858093261718750,
      .omega0    = -1.00367334110330941321,
      .omegadot  = -0.00000000820712757409,
      .w         =  0.45273139092536845984,
      .inc       =  0.94620053286729866038,
      .inc_dot   = -0.00000000036144362701,
      .af0       = -0.00012023560702800751,
      .af1       =  0.00000000000306954462,
      .af2       =  0.0,
      .toc = {
        .tow = 194400.0,
        .wn  = 1897
      },
      .iodc      = 36,
      .iode      = 36
    }
};

/* Real measurements from tracking (Piksi v3 board) channel for L1C/A, SV 5. */
const channel_measurement_t l1ca_meas_in= {
    .sid = {
      .sat = 5,
      .code = 0
    },
    .code_phase_chips = 1022.99972138693556189537,
    .code_phase_rate  = 1023001.17311632633209228516,
    .carrier_phase    = 149059.76833494100719690323,
    .carrier_freq     = 1725.01623535156250000000,
    .time_of_week_ms  = 196039320,
    .receiver_time    = 16.75269665408805153106,
    .snr              = 48.40019989013671875000,
    .lock_counter     = 12594
};

/* Real measurements from tracking (Piksi v3 board) channel for L2CM, SV 5. */
const channel_measurement_t l2cm_meas_in= {
    .sid = {
      .sat = 5,
      .code = 1
    },
    .code_phase_chips = 299.96691101603209972382,
    .code_phase_rate  = 1023001.08898675441741943359,
    .carrier_phase    = 113927.78789247456006705761,
    .carrier_freq     = 1342.59655761718750000000,
    .time_of_week_ms  = 196039320, /* Note: since TOW is not decoded from L2C
                                     stream use it from l1C/A */
    .receiver_time    = 16.75298984654088130242,
    .snr              = 41.51419830322265625000,
    .lock_counter     = 49238
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
              "Position difference out of range: %10.8f m, expected < 0.02",
              check_value);

  check_value = fabs(out_l1ca.pseudorange - out_l2cm.pseudorange);
  fail_unless(check_value < 1e-6,
              "Pseudorange difference out of range: %15.13f m, expected < 1e-6",
              check_value);

  check_value = fabs(out_l1ca.doppler/out_l2cm.doppler - GPS_L1_HZ/GPS_L2_HZ);
  fail_unless(check_value < 0.002,
              "Doppler rate out of range: %15.10f, expected < 0.002",
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
