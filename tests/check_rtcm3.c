
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <check.h>
#include "check_utils.h"

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

START_TEST(test_rtcm3_read_write_header)
{
  u8 buff[22];

  gps_time_t t = {
    .wn = 22,
    .tow = 22.222
  };

  rtcm3_write_header(buff, 1234, 2269, t, 1, 22, 1, 6);

  u16 type, id;
  u8 sync, n_sat, div_free, smooth;
  double tow;

  rtcm3_read_header(buff, &type, &id, &tow, &sync, &n_sat, &div_free, &smooth);

  fail_unless(type == 1234, "type decode error, decoded %d, expected 1234", type);
  fail_unless(id == 2269, "id decode error, decoded %d, expected 2269", id);
  fail_unless(fabs(tow - t.tow) < 1e-3, "TOW decode error, decoded %f, expected %f", tow, t.tow);
  fail_unless(sync == 1, "id decode error, decoded %d, expected 1", id);
  fail_unless(n_sat == 22, "n_sat decode error, decoded %d, expected 22", n_sat);
  fail_unless(div_free == 1, "div_free decode error, decoded %d, expected 1", div_free);
  fail_unless(smooth == 6, "smooth decode error, decoded %d, expected 6", smooth);
}
END_TEST


START_TEST(test_rtcm3_encode_decode)
{
  navigation_measurement_t nm_orig[22];
  navigation_measurement_t nm[22];

  seed_rng();

  for (u8 i=0; i<22; i++) {
    nm[i].prn = i;
    nm[i].raw_pseudorange = frand(19e6, 21e6);
    nm[i].carrier_phase = frand(-5e5, 5e5);
    nm[i].lock_time = frand(0, 1000);
    nm[i].snr = frand(0, 20);
  }

  memcpy(nm_orig, nm, sizeof(nm));

  gps_time_t t = {
    .wn = 1234,
    .tow = frand(0, 604800)
  };

  u8 buff[355];

  rtcm3_encode_1002(buff, 1234, t, 22, nm, 0);

  navigation_measurement_t nm_out[22];
  double tow_out;
  u8 sync, n_sat = 222;
  u16 id;

  s8 ret = rtcm3_decode_1002(buff, &id, &tow_out, &n_sat, 0, &sync);

  fail_unless(ret >= 0, "rtcm3_decode_1002 returned an error (%d)", ret);
  fail_unless(id == 1234, "decoded station id as %d, expected 1234", id);
  fail_unless(n_sat == 22, "decoded n_sat as %d, expected 22", n_sat);
  fail_unless(fabs(tow_out - t.tow) < 1e-3, "decoded TOW as %f, expected %f, error %f",
      tow_out, t.tow, tow_out - t.tow);

  ret = rtcm3_decode_1002(buff, &id, &tow_out, &n_sat, nm_out, &sync);

  for (u8 i=0; i<22; i++) {
    fail_unless(nm[i].prn == nm_out[i].prn, "[%d] PRNs not equal - "
        "decoded %d, expected %d", i, nm_out[i].prn, nm[i].prn);

    double pr_err = nm[i].raw_pseudorange - nm_out[i].raw_pseudorange;
    fail_unless(fabs(pr_err) < 0.02, "[%d] pseudorange error > 0.04m - "
        "decoded %f, expected %f, error %f", i, nm_out[i].raw_pseudorange, nm[i].raw_pseudorange, pr_err);

    double carr_err = nm[i].carrier_phase - nm_out[i].carrier_phase;
    fail_unless(fabs(carr_err) < 0.003, "carrier phase error (fractional part) > 0.003 cycles - "
        "[%d] decoded %f, expected %f, error %f", i, nm_out[i].carrier_phase, nm[i].carrier_phase, carr_err);

    double snr_err = nm[i].snr - nm_out[i].snr;
    /* Calculate error bound on SNR given logarithmic error bound on CNR. */
    double err_bound = nm[i].snr * (pow(10.0, 1.0 / 40.0) - 1);
    fail_unless(fabs(snr_err) < err_bound, "SNR error > 0.003 - "
        "[%d] decoded %f, expected %f, error %f, bound %f", i, nm_out[i].snr, nm[i].snr, snr_err, err_bound);

    fail_unless((nm_out[i].lock_time == 0) && (nm[i].lock_time == 0),
        "lock time should be zero when adjusting int. amb. - [%d] decoded %f",
        i, nm_out[i].lock_time, nm[i].lock_time);

    double cp_adj = nm[i].carrier_phase - nm_orig[i].carrier_phase;
    fail_unless(fmod(cp_adj, 1.0) == 0,
        "carrier phase adjusted by non integer amount %f -> %f (%f)",
        nm_orig[i].carrier_phase, nm[i].carrier_phase, cp_adj);
  }

  /* Re-encode after adjustment, now there should be no further adjustment and
   * the lock time should be correct. */

  for (u8 i=0; i<22; i++)
    nm[i].lock_time = frand(0, 1000);

  rtcm3_encode_1002(buff, 1234, t, 22, nm, 0);
  rtcm3_decode_1002(buff, &id, &tow_out, &n_sat, nm_out, &sync);

  for (u8 i=0; i<22; i++) {
    double cp_adj = nm_out[i].carrier_phase - nm[i].carrier_phase;
    fail_unless(cp_adj < 0.003, "carrier phase re-adjusted %f -> %f (%f)",
        nm[i].carrier_phase, nm_out[i].carrier_phase, cp_adj);

    fail_unless(nm_out[i].lock_time <= nm[i].lock_time,
        "lock time error, should always be less than input lock time - [%d] decoded %f, expected %f",
        i, nm_out[i].lock_time, nm[i].lock_time);
  }
}
END_TEST


Suite* rtcm3_suite(void)
{
  Suite *s = suite_create("RTCMv3");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_rtcm3_check_frame);
  tcase_add_test(tc_core, test_rtcm3_write_frame);
  tcase_add_test(tc_core, test_rtcm3_read_write_header);
  tcase_add_test(tc_core, test_rtcm3_encode_decode);
  suite_add_tcase(s, tc_core);

  return s;
}

