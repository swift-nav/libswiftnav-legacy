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
#define DEBUG 1
#include <stdio.h>
#include <check.h>
#include <math.h>
#include <libswiftnav/logging.h>
#include <libswiftnav/track.h>
#include <libswiftnav/constants.h>

/* real ephemeris for SV 5 at 14:00, 18-May-2016 */
const ephemeris_t e_in = {
    .sid = {
      .sat = 17,
      .code = 0
    },
    .toe = {
      .tow = 309600.0,
      .wn  = 1897
    },
    .ura = 2.0,
    .fit_interval = 14400,
    .valid = 1,
    .healthy = 1,
    .kepler = {
      .tgd       =    -0.00000001071020960808,
      .crc       =   192.18750000000000000000,
      .crs       =    18.15625000000000000000,
      .cuc       =     0.00000095553696155548,
      .cus       =     0.00001021660864353180,
      .cic       =     0.00000025704503059387,
      .cis       =     0.00000025145709514618,
      .dn        =     0.00000000400695261995,
      .m0        =    -1.42679837002772047505,
      .ecc       =     0.01096432015765458345,
      .sqrta     =  5153.69151878356933593750,
      .omega0    =  -3.04227561478505581505,
      .omegadot  =  -0.00000000774889420112,
      .w         =  -1.92851380270647410065,
      .inc       =   0.97664361232661101031,
      .inc_dot   =   0.00000000031429880609,
      .af0       =  -0.00020375801250338554,
      .af1       =  -0.00000000000022737368,
      .af2       =   0.00000000000000000000,
      .toc = {
        .tow = 309600.0,
        .wn  = 1897
      },
      .iodc      = 40,
      .iode      = 40
    }
};

/* Real measurements from tracking (Piksi v3 board) channel for L1C/A, SV 17. */
const channel_measurement_t l1ca_meas_in= {
    .sid = {
      .sat = 17,
      .code = 0
    },
    .code_phase_chips = 1022.97803381178528070450,
    .code_phase_rate =   1022999.43775457143783569336,
    .carrier_phase =   -170000.21982894162647426128,
    .carrier_freq =     -1231.57067871093750000000,
    .time_of_week_ms = 305329331,
    .receiver_time =        12.52904523270440328986,
    .snr =        44.42032241821289062500,
    .lock_counter = 28729

};

/* Real measurements from tracking (Piksi v3 board) channel for L2CM, SV 17. */
const channel_measurement_t l2cm_meas_in= {
    .sid = {
      .sat = 17,
      .code = 1
    },
    .code_phase_chips =     20459.99507837137207388878,
    .code_phase_rate =   1022999.53924232721328735352,
    .carrier_phase =   -130393.67425439949147403240,
    .carrier_freq =      -958.27520751953125000000,
    .time_of_week_ms = 305328546,
    .receiver_time =        12.51804520251572405698,
    .snr =        34.12757873535156250000,
    .lock_counter = 33636
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
  log_debug(/*check_value < 0.02,*/
              "Position difference out of range: %10.8f m, expected < 0.02\n",
              check_value);

  check_value = fabs(out_l1ca.pseudorange - out_l2cm.pseudorange);
  log_debug(/*check_value < 1e-6,*/
              "Pseudorange difference out of range: %15.13f m, expected < 1e-6\n",
              check_value);

  check_value = fabs(out_l1ca.doppler/out_l2cm.doppler - GPS_L1_HZ/GPS_L2_HZ);
  log_debug(/*check_value < 0.002,*/
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
