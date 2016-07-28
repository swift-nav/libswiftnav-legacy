
#include <check.h>

#include  <libswiftnav/almanac.h>

START_TEST(test_almanac_equal)
{
  almanac_t a;
  almanac_t b;

  memset(&a, 0, sizeof(a));
  memset(&b, 0, sizeof(b));

  fail_unless(almanac_equal(&a, &b), "Almanacs should be equal");

  a.valid = 1;
  fail_unless(!almanac_equal(&a, &b),
    "Almanacs should not be equal (valid)");
  memset(&a, 0, sizeof(a));

  a.health_bits = 0x3f;
  fail_unless(!almanac_equal(&a, &b),
    "Almanacs should not be equal (health_bits)");
  memset(&a, 0, sizeof(a));

  a.sid.sat = 1;
  fail_unless(!almanac_equal(&a, &b),
    "Almanacs should not be equal (sid.sat)");
  memset(&a, 0, sizeof(a));

  a.sid.code = 1;
  fail_unless(!almanac_equal(&a, &b),
    "Almanacs should not be equal (sid.band)");
  memset(&a, 0, sizeof(a));

  a.sid.code = 1;
  fail_unless(!almanac_equal(&a, &b),
    "Almanacs should not be equal (sid.constellation)");
  memset(&a, 0, sizeof(a));

  a.toa.wn = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (toa.wn)");
  memset(&a, 0, sizeof(a));

  a.toa.tow = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (toa.tow)");
  memset(&a, 0, sizeof(a));

  a.ura = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (ura)");
  memset(&a, 0, sizeof(a));

  a.fit_interval = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (fit_interval)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.m0 = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (kepler.m0)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.ecc = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (kepler.ecc)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.sqrta = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (kepler.sqrta)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.omega0 = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (kepler.omega0)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.omegadot = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (kepler.omegadot)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.w = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (kepler.w)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.inc = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (kepler.inc)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.af0 = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (kepler.af0)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.af1 = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (kepler.af1)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_SBAS_L1CA;
  a.xyz.pos[0] = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (xyz.pos[0])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_SBAS_L1CA;
  a.xyz.pos[1] = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (xyz.pos[1])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_SBAS_L1CA;
  a.xyz.pos[2] = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (xyz.pos[2])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_SBAS_L1CA;
  a.xyz.vel[0] = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (xyz.vel[0])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_SBAS_L1CA;
  a.xyz.vel[1] = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (xyz.vel[1])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_SBAS_L1CA;
  a.xyz.vel[2] = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (xyz.vel[2])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_SBAS_L1CA;
  a.xyz.acc[0] = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (xyz.acc[0])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_SBAS_L1CA;
  a.xyz.acc[1] = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (xyz.acc[1])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_SBAS_L1CA;
  a.xyz.acc[2] = 1;
  fail_unless(!almanac_equal(&a, &b),
      "Almanacs should not be equal (xyz.acc[2])");
  memset(&a, 0, sizeof(a));
}
END_TEST

Suite* almanac_suite(void)
{
  Suite *s = suite_create("Almanac");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_almanac_equal);
  suite_add_tcase(s, tc_core);

  return s;
}
