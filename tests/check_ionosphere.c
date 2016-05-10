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

#include <check.h>
#include <stdio.h>

#include  <libswiftnav/constants.h>
#include  <libswiftnav/ionosphere.h>

START_TEST(test_calc_ionosphere)
{
  gps_time_t t = {.wn = 1875, .tow = 479820};
  ionosphere_t i = {.a0 = 0.1583e-7, .a1 = -0.7451e-8,
                    .a2 = -0.5960e-7, .a3 = 0.1192e-6,
                    .b0 = 0.1290e6, .b1 = -0.2130e6,
                    .b2 = 0.6554e5, .b3 = 0.3277e6};
  double lat_u = -35.3 * D2R, lon_u = 149.1 * D2R;
  double a = 0.0 * D2R, e = 15.0 * D2R;
  double d_true = 7.202;

  const double d_tol = 1e-3;

  double d_l1 = calc_ionosphere(&t, lat_u, lon_u, a, e, &i);
  double d_err = fabs(d_l1 - d_true);

  fail_unless(d_err < d_tol,
      "Distance didn't match hardcoded correct value %0.5f. Saw: %.5f\n",
      d_true, d_l1);

  t.wn = 1042;
  t.tow = 593100;
  i.a0 = 0.3820e-7;
  i.a1 = 0.1490e-7;
  i.a2 = -0.1790e-6;
  i.a3 = 0.0;
  i.b0 = 0.1430e6;
  i.b1 = 0.0;
  i.b2 = -0.3280e6;
  i.b3 = 0.1130e6;
  lat_u = 40.0 * D2R;
  lon_u = 260.0 * D2R;
  a = 210.0 * D2R;
  e = 20.0 * D2R;
  d_true = 23.784;

  d_l1 = calc_ionosphere(&t, lat_u, lon_u, a, e, &i);
  d_err = fabs(d_l1 - d_true);

  fail_unless(d_err < d_tol,
      "Distance didn't match hardcoded correct values %0.5f. Saw: %.5f\n",
      d_true, d_l1);
}
END_TEST

START_TEST(test_decode_iono_parameters)
{
  #define tol 1e-12
  struct {
    u32 frame_words[8];
    ionosphere_t result;
  } t_case = {
      .frame_words = {
          /* 4th SF real data at 11-May-2016 */
          0x1e0300c9,0x7fff8c24,0x23fbdc2,0,0,0,0,0
      },
      .result = { /* reference data provided by u-blox receiver */
        .a0 = 0.0000000111758,
        .a1 = 0.0000000223517,
        .a2 = -0.0000000596046,
        .a3 = -0.0000001192092,
        .b0 = 98304.0,
        .b1 = 131072.0,
        .b2 = -131072.0,
        .b3 = -589824.0,
      }
  };
  ionosphere_t i;
  decode_iono_parameters(t_case.frame_words, &i);
  fail_unless(fabs(i.a0 - t_case.result.a0) < tol,
              "alfa 0 == %30.20f, expected %30.20f, tolerance = %30.20f",
              i.a0, t_case.result.a0, tol);
  fail_unless(fabs(i.a1 - t_case.result.a1) < tol,
              "alfa 1 == %30.20f, expected %30.20f, tolerance = %30.20f",
              i.a1, t_case.result.a1, tol);
  fail_unless(fabs(i.a2 - t_case.result.a2) < tol,
              "alfa 2 == %30.20f, expected %30.20f, tolerance = %30.20f",
              i.a2, t_case.result.a2, tol);
  fail_unless(fabs(i.a3 - t_case.result.a3) < tol,
              "alfa 3 == %30.20f, expected %30.20f, tolerance = %30.20f",
              i.a3, t_case.result.a3, tol);
  fail_unless(fabs(i.b0 - t_case.result.b0) < tol,
              "beta 0 == %30.20f, expected %30.20f, tolerance = %30.20f",
              i.b0, t_case.result.b0, tol);
  fail_unless(fabs(i.b1 - t_case.result.b1) < tol,
              "beta 1 == %30.20f, expected %30.20f, tolerance = %30.20f",
              i.b1, t_case.result.b1, tol);
  fail_unless(fabs(i.b2 - t_case.result.b2) < tol,
              "beta 2 == %30.20f, expected %30.20f, tolerance = %30.20f",
              i.b2, t_case.result.b2, tol);
  fail_unless(fabs(i.b3 - t_case.result.b3) < tol,
              "beta 3 == %30.20f, expected %30.20f, tolerance = %30.20f",
              i.b3, t_case.result.b3, tol);

}
END_TEST

Suite* ionosphere_suite(void)
{
  Suite *s = suite_create("Ionosphere");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_calc_ionosphere);
  tcase_add_test(tc_core, test_decode_iono_parameters);
  suite_add_tcase(s, tc_core);

  return s;
}
