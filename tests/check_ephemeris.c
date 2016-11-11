
#include <check.h>

#include  <libswiftnav/ephemeris.h>

START_TEST(test_ephemeris_equal)
{
  ephemeris_t a;
  ephemeris_t b;

  memset(&a, 0, sizeof(a));
  memset(&b, 0, sizeof(b));

  fail_unless(ephemeris_equal(&a, &b), "Ephemerides should be equal");

  a.valid = 1;
  fail_unless(!ephemeris_equal(&a, &b),
    "Ephemerides should not be equal (valid)");
  memset(&a, 0, sizeof(a));

  a.health_bits = 1;
  fail_unless(!ephemeris_equal(&a, &b),
    "Ephemerides should not be equal (healthy)");
  memset(&a, 0, sizeof(a));

  a.sid.sat = 1;
  fail_unless(!ephemeris_equal(&a, &b),
    "Ephemerides should not be equal (sid.sat)");
  memset(&a, 0, sizeof(a));

  a.sid.code = 1;
  fail_unless(!ephemeris_equal(&a, &b),
    "Ephemerides should not be equal (sid.band)");
  memset(&a, 0, sizeof(a));

  a.sid.code = 1;
  fail_unless(!ephemeris_equal(&a, &b),
    "Ephemerides should not be equal (sid.constellation)");
  memset(&a, 0, sizeof(a));

  a.toe.wn = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (toe.wn)");
  memset(&a, 0, sizeof(a));

  a.toe.tow = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (toe.tow)");
  memset(&a, 0, sizeof(a));

  a.ura = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (ura)");
  memset(&a, 0, sizeof(a));

  a.fit_interval = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (fit_interval)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.tgd = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.tgd)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.crs = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.crs)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.crc = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.crc)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.cuc = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.cuc)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.cus = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.cus)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.cic = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.cic)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.cis = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.cis)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.dn = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.dn)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.m0 = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.m0)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.ecc = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.ecc)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.sqrta = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.sqrta)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.omega0 = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.omega0)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.omegadot = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.omegadot)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.w = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.w)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.inc = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.inc)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.inc_dot = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.inc_dot)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.af0 = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.af0)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.af1 = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.af1)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.af2 = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.af2)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.iode = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.iode)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.iodc = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.iodc)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.toc.wn = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.toc.wn)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GPS_L1CA;
  a.kepler.toc.tow = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (kepler.toc.tow)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_SBAS_L1CA;
  a.xyz.pos[0] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (xyz.pos[0])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_SBAS_L1CA;
  a.xyz.pos[1] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (xyz.pos[1])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_SBAS_L1CA;
  a.xyz.pos[2] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (xyz.pos[2])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_SBAS_L1CA;
  a.xyz.vel[0] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (xyz.vel[0])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_SBAS_L1CA;
  a.xyz.vel[1] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (xyz.vel[1])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_SBAS_L1CA;
  a.xyz.vel[2] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (xyz.vel[2])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_SBAS_L1CA;
  a.xyz.acc[0] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (xyz.acc[0])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_SBAS_L1CA;
  a.xyz.acc[1] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (xyz.acc[1])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_SBAS_L1CA;
  a.xyz.acc[2] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (xyz.acc[2])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_SBAS_L1CA;
  a.xyz.a_gf0 = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (xyz.a_gf0)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_SBAS_L1CA;
  a.xyz.a_gf1 = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (xyz.a_gf1)");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GLO_L1CA;
  a.glo.pos[0] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (glo.pos[0])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GLO_L1CA;
  a.glo.pos[1] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (glo.pos[1])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GLO_L1CA;
  a.glo.pos[2] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (glo.pos[2])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GLO_L1CA;
  a.glo.vel[0] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (glo.vel[0])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GLO_L1CA;
  a.glo.vel[1] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (glo.vel[1])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GLO_L1CA;
  a.glo.vel[2] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (glo.vel[2])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GLO_L1CA;
  a.glo.acc[0] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (glo.acc[0])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GLO_L1CA;
  a.glo.acc[1] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (glo.acc[1])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GLO_L1CA;
  a.glo.acc[2] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (glo.acc[2])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GLO_L2CA;
  a.glo.pos[0] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (glo.pos[0])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GLO_L2CA;
  a.glo.pos[1] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (glo.pos[1])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GLO_L2CA;
  a.glo.pos[2] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (glo.pos[2])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GLO_L2CA;
  a.glo.vel[0] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (glo.vel[0])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GLO_L2CA;
  a.glo.vel[1] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (glo.vel[1])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GLO_L2CA;
  a.glo.vel[2] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (glo.vel[2])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GLO_L2CA;
  a.glo.acc[0] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (glo.acc[0])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GLO_L2CA;
  a.glo.acc[1] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (glo.acc[1])");
  memset(&a, 0, sizeof(a));

  a.sid.code = CODE_GLO_L2CA;
  a.glo.acc[2] = 1;
  fail_unless(!ephemeris_equal(&a, &b),
      "Ephemerides should not be equal (glo.acc[2])");
  memset(&a, 0, sizeof(a));
}
END_TEST

START_TEST(test_signal_component_health)
{
  return;
  const struct test_case
  {
    ephemeris_t e;
    bool healthy;
  } test_cases[] = {
    {
      .e = {
        .sid = {.code = CODE_GPS_L1CA, .sat = 1},
        .health_bits = 0,
        .ura = 2.0,
        .valid = true},
      .healthy = true
    },
    {
      .e = {
        .sid = {.code = CODE_GPS_L1CA, .sat = 1},
        .health_bits = 0,
        .ura = 200.0,
        .valid = false},
      .healthy = true
    },
    {
      .e = {
        .sid = {.code = CODE_GPS_L2CM, .sat = 32},
        .health_bits = 0,
        .ura = 2000.0,
        .valid = true},
      .healthy = true
    },
    {
      .e = {
        .sid = {.code = CODE_GPS_L1CA, .sat = 1},
        .health_bits = 0,
        .ura = 33333.0,
        .valid = true},
      .healthy = false
    },
    {
      .e = {
        .sid = {.code = CODE_GPS_L1CA, .sat = 1},
        .health_bits = 0,
        .ura = INVALID_GPS_URA_VALUE,
        .valid = true},
      .healthy = false
    },
    {
      .e = {
        .sid = {.code = CODE_GPS_L1CA, .sat = 1},
        .health_bits = 0,
        .ura = -100.0,
        .valid = true},
      .healthy = false
    },
    {
      .e = {
        .sid = {.code = CODE_GPS_L1CA, .sat = 11},
        .health_bits = 0x3F,
        .ura = 1.0,
        .valid = true},
      .healthy = false
    },
    {
      .e = {
        .sid = {.code = CODE_GPS_L1CA, .sat = 12},
        .health_bits = 0x2A,
        .ura = 1.0,
        .valid = true},
      .healthy = false
    },
    {
      .e = {
        .sid = {.code = CODE_GPS_L2CM, .sat = 13},
        .health_bits = 0x2E,
        .ura = 4000.0,
        .valid = true},
      .healthy = false
    },
    {
      .e = {
        .sid = {.code = CODE_GPS_L2CM, .sat = 1},
        .health_bits = 0x2A,
        .ura = 4000.0,
        .valid = true},
      .healthy = true
    },
    {
      .e = {
        .sid = {.code = CODE_GPS_L1CA, .sat = 22},
        .health_bits = 0x2E,
        .ura = 10.0,
        .valid = true},
      .healthy = true
    },
  };

  for (u32 i=0; i<sizeof(test_cases) / sizeof(test_cases[0]); i++) {
    const struct test_case *t = &test_cases[i];
    fail_unless(t->healthy == signal_healthy(&t->e, t->e.sid.code),
      "test signal %d healthy incorrect", i);
  }
}
END_TEST

Suite* ephemeris_suite(void)
{
  Suite *s = suite_create("Ephemeris");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_ephemeris_equal);
  tcase_add_test(tc_core, test_signal_component_health);
  suite_add_tcase(s, tc_core);

  return s;
}
