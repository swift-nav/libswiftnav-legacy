#include <check.h>
#include <stdio.h>
#include <math.h>
#include "check_utils.h"

#include "pvt.h"

static navigation_measurement_t nm1 = {
  .prn = 9,
  .pseudorange = 23946993.888943646,
  .sat_pos = {-19477278.087422125, -7649508.9457812719, 16674633.163554827}
};

static navigation_measurement_t nm2 = {
  .prn = 1,
  .pseudorange = 22932174.156858064,
  .sat_pos = {-9680013.5408340245, -15286326.354385279, 19429449.383770257},
};

static navigation_measurement_t nm3 = {
  .prn = 2,
  .pseudorange = 24373231.648055989,
  .sat_pos = {-19858593.085281931, -3109845.8288993631, 17180320.439503901},
};

static navigation_measurement_t nm4 = {
  .prn = 3,
  .pseudorange = 24779663.252316438,
  .sat_pos = {6682497.8716542246, -14006962.389166718, 21410456.275678463},
};

static navigation_measurement_t nm5 = {
  .prn = 4,
  .pseudorange = 26948717.022331879,
  .sat_pos = {7415370.9916331079, -24974079.044485383, -3836019.0262199985},
};

static navigation_measurement_t nm6 = {
  .prn = 5,
  .pseudorange = 23327405.435463827,
  .sat_pos = {-2833466.1648670658, -22755197.793894723, 13160322.082875408},
};

static navigation_measurement_t nm7 = {
  .prn = 6,
  .pseudorange = 27371419.016328193,
  .sat_pos = {14881660.383624561, -5825253.4316490609, 21204679.68313824},
};

static navigation_measurement_t nm8 = {
  .prn = 7,
  .pseudorange = 26294221.697782904,
  .sat_pos = {12246530.477279386, -22184711.955107089, 7739084.2855069181},
};

static navigation_measurement_t nm9 = {
  .prn = 8,
  .pseudorange = 25781999.479948733,
  .sat_pos = {-25360766.249484103, -1659033.490658124, 7821492.0398916304},
};


START_TEST(test_pvt_failed_repair)
{
  u8 n_used = 5;
  gnss_solution soln;
  dops_t dops;

  navigation_measurement_t nms[9] = {nm1, nm2, nm3, nm4, nm5, nm6, nm7, nm8};

  calc_PVT(n_used, nms, &soln, &dops);
  /* PVT repair requires at least 6 measurements. */
  fail_unless(soln.valid == 0, "Solution should be invalid!");
}
END_TEST
START_TEST(test_pvt_repair)
{
  u8 n_used = 6;
  gnss_solution soln;
  dops_t dops;

  navigation_measurement_t nms[9] = {nm1, nm2, nm3, nm4, nm5, nm6, nm7, nm8, nm9};

  s8 code = calc_PVT(n_used, nms, &soln, &dops);
  fail_unless(code == 1, "Return code should be 1 (pvt repair). Saw: %d\n", code);
  fail_unless(soln.n_used == n_used - 1, "PVT solver failed to repair solution.");
}
END_TEST

Suite* pvt_test_suite(void)
{
  Suite *s = suite_create("PVT Solver");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_pvt_repair);
  tcase_add_test(tc_core, test_pvt_failed_repair);
  suite_add_tcase(s, tc_core);

  return s;
}

