#include <check.h>
#include <stdio.h>

#include <libswiftnav/observation.h>

navigation_measurement_t nm1 = {
  .sid = {.sat = 1},
  .raw_pseudorange = 11,
  .carrier_phase = 12,
  .raw_doppler = 13,
  .snr = 14,
  .sat_pos = {1, 2, 3},
  .sat_vel = {4, 5, 6},
};
navigation_measurement_t nm1_2 = {
  .sid = {.sat = 1},
  .raw_pseudorange = 111,
  .carrier_phase = 112,
  .raw_doppler = 113,
  .snr = 114,
  .sat_pos = {7, 8, 9},
  .sat_vel = {10, 11, 12},
};
navigation_measurement_t nm2 = {
  .sid = {.sat = 2},
  .raw_pseudorange = 21,
  .carrier_phase = 22,
  .raw_doppler = 23,
  .snr = 224,
  .sat_pos = {1, 2, 3},
  .sat_vel = {4, 5, 6},
};
navigation_measurement_t nm2_2 = {
  .sid = {.sat = 2},
  .raw_pseudorange = 221,
  .carrier_phase = 222,
  .raw_doppler = 223,
  .snr = 24,
  .sat_pos = {13, 14, 15},
  .sat_vel = {16, 17, 18},
};
navigation_measurement_t nm3 = {
  .sid = {.sat = 3},
  .raw_pseudorange = 31,
  .carrier_phase = 32,
  .raw_doppler = 33,
  .snr = 34,
  .sat_pos = {1, 2, 3},
  .sat_vel = {4, 5, 6},
};
navigation_measurement_t nm4 = {
  .sid = {.sat = 4},
  .raw_pseudorange = 41,
  .carrier_phase = 42,
  .raw_doppler = 43,
  .snr = 44,
  .sat_pos = {1, 2, 3},
  .sat_vel = {4, 5, 6},
};

navigation_measurement_t nms1[2];
navigation_measurement_t nms2[2];

START_TEST(test_single_diff_1)
{
    sdiff_t sds_out[6];

    /* Test for when they are interleaved */
    nms1[0] = nm1; nms1[1] = nm3;
    nms2[0] = nm2; nms2[1] = nm4;

    u8 num_match = single_diff(2, nms1, 2, nms2, sds_out);
    fail_unless(num_match == 0);

    /* Test when one set follows the other */
    nms1[0] = nm1; nms1[1] = nm2;
    nms2[0] = nm3; nms2[1] = nm4;

    num_match = single_diff(2, nms1, 2, nms2, sds_out);
    fail_unless(num_match == 0);

    /* Test it the other way */
    nms2[0] = nm1; nms2[1] = nm2;
    nms1[0] = nm3; nms1[1] = nm4;

    num_match = single_diff(2, nms1, 2, nms2, sds_out);
    fail_unless(num_match == 0);
}
END_TEST

START_TEST(test_single_diff_2)
{
    sdiff_t sds_out[3];

    /* Test for when they both have two */
    nms1[0] = nm1; nms1[1] = nm2;
    nms2[0] = nm1; nms2[1] = nm3;

    u8 num_match = single_diff(2, nms1, 2, nms2, sds_out);

    fail_unless(num_match == 1);
    fail_unless(sds_out[0].sid.sat == 1);
    fail_unless(sds_out[0].pseudorange == 0);
    fail_unless(sds_out[0].carrier_phase == 0);
    fail_unless(sds_out[0].doppler == 0);
    fail_unless(sds_out[0].snr == 14);
    fail_unless(memcmp(sds_out[0].sat_pos, nm1.sat_pos, sizeof(nm1.sat_pos)) == 0);
    fail_unless(memcmp(sds_out[0].sat_vel, nm1.sat_vel, sizeof(nm1.sat_vel)) == 0);

    /* Test for both with two the other way */
    nms1[0] = nm1; nms1[1] = nm3;
    nms2[0] = nm1_2; nms2[1] = nm2;

    num_match = single_diff(2, nms1, 2, nms2, sds_out);

    fail_unless(num_match == 1);
    fail_unless(sds_out[0].sid.sat == 1);
    fail_unless(sds_out[0].pseudorange == -100);
    fail_unless(sds_out[0].carrier_phase == -100);
    fail_unless(sds_out[0].doppler == -100);
    fail_unless(sds_out[0].snr == 14);
    fail_unless(memcmp(sds_out[0].sat_pos, nm1_2.sat_pos, sizeof(nm1_2.sat_pos)) == 0);
    fail_unless(memcmp(sds_out[0].sat_vel, nm1_2.sat_vel, sizeof(nm1_2.sat_vel)) == 0);

    /* Test when one has only one */
    nms1[0] = nm1_2;
    nms2[0] = nm1; nms2[1] = nm2;

    num_match = single_diff(1, nms1, 2, nms2, sds_out);

    fail_unless(num_match == 1);
    fail_unless(sds_out[0].sid.sat == 1);
    fail_unless(sds_out[0].pseudorange == 100);
    fail_unless(sds_out[0].carrier_phase == 100);
    fail_unless(sds_out[0].doppler == 100);
    fail_unless(sds_out[0].snr == 14);
    fail_unless(memcmp(sds_out[0].sat_pos, nm1.sat_pos, sizeof(nm1.sat_pos)) == 0);
    fail_unless(memcmp(sds_out[0].sat_vel, nm1.sat_vel, sizeof(nm1.sat_vel)) == 0);

    /* Test when the other has only one */
    nms1[0] = nm1; nms1[1] = nm2;
    nms2[0] = nm1;

    num_match = single_diff(2, nms1, 1, nms2, sds_out);

    fail_unless(num_match == 1);
    fail_unless(sds_out[0].sid.sat == 1);
    fail_unless(sds_out[0].pseudorange == 0);
    fail_unless(sds_out[0].carrier_phase == 0);
    fail_unless(sds_out[0].doppler == 0);
    fail_unless(sds_out[0].snr == 14);
    fail_unless(memcmp(sds_out[0].sat_pos, nm1.sat_pos, sizeof(nm1.sat_pos)) == 0);
    fail_unless(memcmp(sds_out[0].sat_vel, nm1.sat_vel, sizeof(nm1.sat_vel)) == 0);

    /* Test when they both match */
    nms1[0] = nm1_2; nms1[1] = nm2;
    nms2[0] = nm1; nms2[1] = nm2_2;

    num_match = single_diff(2, nms1, 2, nms2, sds_out);

    fail_unless(num_match == 2);

    fail_unless(sds_out[0].sid.sat == 1);
    fail_unless(sds_out[0].pseudorange == 100);
    fail_unless(sds_out[0].carrier_phase == 100);
    fail_unless(sds_out[0].doppler == 100);
    fail_unless(sds_out[0].snr == 14);
    fail_unless(memcmp(sds_out[0].sat_pos, nm1.sat_pos, sizeof(nm1.sat_pos)) == 0);
    fail_unless(memcmp(sds_out[0].sat_vel, nm1.sat_vel, sizeof(nm1.sat_vel)) == 0);

    fail_unless(sds_out[1].sid.sat == 2);
    fail_unless(sds_out[1].pseudorange == -200);
    fail_unless(sds_out[1].carrier_phase == -200);
    fail_unless(sds_out[1].doppler == -200);
    fail_unless(sds_out[1].snr == 24);
    fail_unless(memcmp(sds_out[1].sat_pos, nm2_2.sat_pos, sizeof(nm2_2.sat_pos)) == 0);
    fail_unless(memcmp(sds_out[1].sat_vel, nm2_2.sat_vel, sizeof(nm2_2.sat_vel)) == 0);
}
END_TEST

START_TEST(test_single_diff_3)
{
    sdiff_t sds_out[3];

    /* Test for when they both have two */
    nms1[0] = nm1; nms1[1] = nm3;
    nms2[0] = nm2; nms2[1] = nm3;

    u8 num_match = single_diff(2, nms1, 2, nms2, sds_out);
    fail_unless(num_match == 1);
    fail_unless(sds_out[0].sid.sat == 3);

    /* Test for both with two the other way */
    nms1[0] = nm2; nms1[1] = nm3;
    nms2[0] = nm1; nms2[1] = nm3;

    num_match = single_diff(2, nms1, 2, nms2, sds_out);

    fail_unless(num_match == 1);
    fail_unless(sds_out[0].sid.sat == 3);

    /* Test when one has only one */
    nms1[0] = nm2;
    nms2[0] = nm1; nms2[1] = nm2;

    num_match = single_diff(1, nms1, 2, nms2, sds_out);

    fail_unless(num_match == 1);
    fail_unless(sds_out[0].sid.sat == 2);

    /* Test when the other has only one */
    nms1[0] = nm1; nms1[1] = nm2;
    nms2[0] = nm2;

    num_match = single_diff(2, nms1, 1, nms2, sds_out);

    fail_unless(num_match == 1);
    fail_unless(sds_out[0].sid.sat == 2);
}
END_TEST

Suite* observation_test_suite(void)
{
  Suite *s = suite_create("Observation Handling");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_single_diff_1);
  tcase_add_test(tc_core, test_single_diff_2);
  tcase_add_test(tc_core, test_single_diff_3);
  suite_add_tcase(s, tc_core);

  return s;
}
