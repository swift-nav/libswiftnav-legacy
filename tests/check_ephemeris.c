
#include <check.h>

#include <ephemeris.h>

START_TEST(test_ephemeris_equal)
{
  ephemeris_t a;
  ephemeris_t b;

  memset(&a, 0, sizeof(a));
  memset(&b, 0, sizeof(b));

  fail_unless(ephemeris_equal(&a, &b), "Ephemerides should be equal");

  a.tgd = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (tgd)");
  memset(&a, 0, sizeof(a));

  a.tgd = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (tgd)");
  memset(&a, 0, sizeof(a));

  a.crs = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (crs)");
  memset(&a, 0, sizeof(a));

  a.crc = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (crc)");
  memset(&a, 0, sizeof(a));

  a.cuc = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (cuc)");
  memset(&a, 0, sizeof(a));

  a.cus = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (cus)");
  memset(&a, 0, sizeof(a));

  a.cic = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (cic)");
  memset(&a, 0, sizeof(a));

  a.cis = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (cis)");
  memset(&a, 0, sizeof(a));

  a.dn = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (dn)");
  memset(&a, 0, sizeof(a));

  a.m0 = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (m0)");
  memset(&a, 0, sizeof(a));

  a.ecc = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (ecc)");
  memset(&a, 0, sizeof(a));

  a.sqrta = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (sqrta)");
  memset(&a, 0, sizeof(a));

  a.omega0 = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (omega0)");
  memset(&a, 0, sizeof(a));

  a.omegadot = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (omegadot)");
  memset(&a, 0, sizeof(a));

  a.w = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (w)");
  memset(&a, 0, sizeof(a));

  a.inc = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (inc)");
  memset(&a, 0, sizeof(a));

  a.inc_dot = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (inc_dot)");
  memset(&a, 0, sizeof(a));

  a.af0 = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (af0)");
  memset(&a, 0, sizeof(a));

  a.af1 = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (af1)");
  memset(&a, 0, sizeof(a));

  a.af2 = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (af2)");
  memset(&a, 0, sizeof(a));

  a.valid = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (valid)");
  memset(&a, 0, sizeof(a));

  a.healthy = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (healthy)");
  memset(&a, 0, sizeof(a));

  a.prn = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (prn)");
  memset(&a, 0, sizeof(a));

  a.iode = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (iode)");
  memset(&a, 0, sizeof(a));

  a.toe.wn = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (toe.wn)");
  memset(&a, 0, sizeof(a));

  a.toe.tow = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (toe.tow)");
  memset(&a, 0, sizeof(a));

  a.toc.wn = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (toc.wn)");
  memset(&a, 0, sizeof(a));

  a.toc.tow = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (toc.tow)");
  memset(&a, 0, sizeof(a));
}
END_TEST

Suite* ephemeris_suite(void)
{
  Suite *s = suite_create("Ephemeris");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_ephemeris_equal);
  suite_add_tcase(s, tc_core);

  return s;
}

