
#include <check.h>

#include <libswiftnav/constants.h>
#include <libswiftnav/troposphere.h>
#include <libswiftnav/time.h>

START_TEST(test_calc_troposphere)
{
	const double d_tol = 1e-4;

	/* some tests against "true" values computed with UNB3M.f */
	/* http://www2.unb.ca/gge/Personnel/Santos/UNB_pack.pdf */

	double lat = 40 * D2R;
	double h = 1300.0;
	double doy = 32.5;
	double el = 45 * D2R;
	double d_true = 2.8567;

	/* GPS week 1669 starts on 1.1.2012, so easier to generate given doy */
	gps_time_t t = {.wn = 1669, .tow = doy*DAY_SECS};
	normalize_gps_time(&t);

	double d_tropo = calc_troposphere(&t, lat, h, el);

	fail_unless(fabs(d_tropo - d_true) < d_tol,
		      "Distance didn't match hardcoded correct values %0.5f. Saw: %.5f\n",
		      d_true, d_tropo);

	lat = -10 * D2R;
	h = 0.0;
	doy = 180.5;
	el = 20 * D2R;
	d_true = 7.4942;

	t.wn = 1669;
	t.tow = doy*DAY_SECS;
	normalize_gps_time(&t);

	d_tropo = calc_troposphere(&t, lat, h, el);

	fail_unless(fabs(d_tropo - d_true) < d_tol,
		      "Distance didn't match hardcoded correct values %0.5f. Saw: %.5f\n",
		      d_true, d_tropo);

	lat = 75 * D2R;
	h = 0.0;
	doy = 50.5;
	el = 10 * D2R;
	d_true = 12.9004;

	t.wn = 1669;
	t.tow = doy*DAY_SECS;
	normalize_gps_time(&t);

	d_tropo = calc_troposphere(&t, lat, h, el);

	fail_unless(fabs(d_tropo - d_true) < d_tol,
		      "Distance didn't match hardcoded correct values %0.5f. Saw: %.5f\n",
		      d_true, d_tropo);
}
END_TEST

Suite* troposphere_suite(void)
{
  Suite *s = suite_create("Troposphere");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_calc_troposphere);
  suite_add_tcase(s, tc_core);

  return s;
}
