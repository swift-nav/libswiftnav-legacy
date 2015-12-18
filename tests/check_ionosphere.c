
#include <check.h>

#include <constants.h>
#include <ionosphere.h>

START_TEST(test_calc_ionosphere)
{
  gps_time_t t = {.wn = 1875, .tow = 479820};
  ionosphere_t i = {.a0 = 0.1583e-7, .a1 = -0.7451e-8, .a2 = -0.5960e-7, .a3 = 0.1192e-6,
                    .b0 = 0.1290e6, .b1 = -0.2130e6, .b2 = 0.6554e5, .b3 = 0.3277e6};
  double lat_u = -35.3 * D2R, lon_u = 149.1 * D2R;
  double a = 0.0 * D2R, e = 15.0 * D2R;
  double d_true = 7.202;

  const double d_tol = 1e-3;

  double d_l1 = calc_ionosphere(t, lat_u, lon_u, a, e, &i);
  double d_err = fabs(d_l1 - d_true);

  fail_unless(d_err < d_tol,
      "Distance didn't match hardcoded correct values. Saw: %.5f\n", d_l1);

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

  d_l1 = calc_ionosphere(t, lat_u, lon_u, a, e, &i);
  d_err = fabs(d_l1 - d_true);

  fail_unless(d_err < d_tol,
      "Distance didn't match hardcoded correct values. Saw: %.5f\n", d_l1);
}
END_TEST

Suite* ionosphere_suite(void)
{
  Suite *s = suite_create("Ionosphere");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_calc_ionosphere);
  suite_add_tcase(s, tc_core);

  return s;
}
