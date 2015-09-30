#include <check.h>
#include <stdio.h>
#include <math.h>
#include "check_utils.h"

#include "pvt.h"

static navigation_measurement_t nm1 = {
  .sid.prn = 9,
  .pseudorange = 23946993.888943646,
  .sat_pos = {-19477278.087422125, -7649508.9457812719, 16674633.163554827}
};

static navigation_measurement_t nm2 = {
  .sid.prn = 1,
  .pseudorange = 22932174.156858064,
  .sat_pos = {-9680013.5408340245, -15286326.354385279, 19429449.383770257},
};

static navigation_measurement_t nm3 = {
  .sid.prn = 2,
  .pseudorange = 24373231.648055989,
  .sat_pos = {-19858593.085281931, -3109845.8288993631, 17180320.439503901},
};

static navigation_measurement_t nm4 = {
  .sid.prn = 3,
  .pseudorange = 24779663.252316438,
  .sat_pos = {6682497.8716542246, -14006962.389166718, 21410456.275678463},
};

static navigation_measurement_t nm5 = {
  .sid.prn = 4,
  .pseudorange = 26948717.022331879,
  .sat_pos = {7415370.9916331079, -24974079.044485383, -3836019.0262199985},
};

static navigation_measurement_t nm6 = {
  .sid.prn = 5,
  .pseudorange = 23327405.435463827,
  .sat_pos = {-2833466.1648670658, -22755197.793894723, 13160322.082875408},
};

static navigation_measurement_t nm7 = {
  .sid.prn = 6,
  .pseudorange = 27371419.016328193,
  .sat_pos = {14881660.383624561, -5825253.4316490609, 21204679.68313824},
};

static navigation_measurement_t nm8 = {
  .sid.prn = 7,
  .pseudorange = 26294221.697782904,
  .sat_pos = {12246530.477279386, -22184711.955107089, 7739084.2855069181},
};

static navigation_measurement_t nm9 = {
  .sid.prn = 8,
  .pseudorange = 25781999.479948733,
  .sat_pos = {-25360766.249484103, -1659033.490658124, 7821492.0398916304},
};


START_TEST(test_pvt_failed_repair)
{
  u8 n_used = 5;
  gnss_solution soln;
  dops_t dops;

  navigation_measurement_t nms[9] = {nm1, nm2, nm3, nm4, nm5, nm6, nm7, nm8};

  calc_PVT(n_used, nms, false, &soln, &dops);
  /* PVT repair requires at least 6 measurements. */
  fail_unless(soln.valid == 0, "Solution should be invalid!");
}
END_TEST

START_TEST(test_pvt_repair)
{
  u8 n_used = 6;
  gnss_solution soln;
  dops_t dops;

  navigation_measurement_t nms[9] =
    {nm1, nm2, nm3, nm4, nm5, nm6, nm7, nm8, nm9};

  s8 code = calc_PVT(n_used, nms, false, &soln, &dops);
  fail_unless(code == 1,
    "Return code should be 1 (pvt repair). Saw: %d\n", code);
  fail_unless(soln.n_used == n_used - 1,
    "PVT solver failed to repair solution.");
}
END_TEST

START_TEST(test_disable_pvt_raim)
{
  u8 n_used = 6;
  gnss_solution soln;
  dops_t dops;

  navigation_measurement_t nms[9] =
    {nm1, nm2, nm3, nm4, nm5, nm6, nm7, nm8, nm9};

  /* disable raim check */
  s8 code = calc_PVT(n_used, nms, true, &soln, &dops);
  fail_unless(code == 2,
    "Return code should be 2 (raim not used). Saw: %d\n", code);
  fail_unless(soln.valid == 1,
    "Solution should be valid!");
}
END_TEST

START_TEST(test_dops)
{
  u8 n_used = 6;
  gnss_solution soln;
  dops_t dops = {.pdop = 22, .gdop = 22, .tdop = 22, .hdop = 22, .vdop = 22};
  dops_t truedops = {.pdop = 2.69955, .gdop = 3.07696, .tdop = 1.47652,
                     .hdop = 1.75922, .vdop = 2.04761};

  const double dop_tol = 1e-3;
  
  navigation_measurement_t nms[6] =
    {nm1, nm2, nm3, nm4, nm5, nm6};

  /* disable raim check */
  s8 code = calc_PVT(n_used, nms, false, &soln, &dops);
  fail_unless(code >= 0,
    "Return code should be >=0 (success). Saw: %d\n", code);
  fail_unless(soln.valid == 1,
    "Solution should be valid!");
  fail_unless(fabs(dops.pdop * dops.pdop -
                   (dops.vdop * dops.vdop + dops.hdop * dops.hdop))
              < dop_tol,
              "HDOP^2 + VDOP^2 != PDOP^2.  Saw: %.5f, %.5f, %.5f, %.5f, %.5f\n",
              dops.pdop, dops.gdop, dops.tdop, dops.hdop, dops.vdop);
  double dop_err = fabs(dops.pdop - truedops.pdop)
                 + fabs(dops.gdop - truedops.gdop)
                 + fabs(dops.tdop - truedops.tdop)
                 + fabs(dops.hdop - truedops.hdop)
                 + fabs(dops.vdop - truedops.vdop);
  fail_unless(dop_err < dop_tol,
              "DOPs don't match hardcoded correct values.  "
              "Saw: %.5f, %.5f, %.5f, %.5f, %.5f\n",
              dops.pdop, dops.gdop, dops.tdop, dops.hdop, dops.vdop);
}
END_TEST


Suite* pvt_test_suite(void)
{
  Suite *s = suite_create("PVT Solver");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_pvt_repair);
  tcase_add_test(tc_core, test_pvt_failed_repair);
  tcase_add_test(tc_core, test_disable_pvt_raim);
  tcase_add_test(tc_core, test_dops);
  suite_add_tcase(s, tc_core);

  return s;
}

