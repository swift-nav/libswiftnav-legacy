#include <math.h>
#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "check_utils.h"

#include <libswiftnav/track.h>

START_TEST(test_costas_discriminator)
{
  const struct
  {
    float i, q;
    float res;
    bool valid;
  } test_cases[] = {
      {
        .i = 0,
        .res = 0,
        .valid = true
      },
      {
        .i = 0.1,
        .q = 0,
        .res = 0,
        .valid = true
      },
      {
        .i = 1.0,
        .q = M_PI,
        .res = 0.20095336902385319,
        .valid = true
      },
      {
        .i = 2.0,
        .q = -M_PI,
        .res = -0.20095336902385319,
        .valid = false
      },
      {
        .i = 1.0,
        .q = -M_PI,
        .res = -0.20095336902385319,
        .valid = true
      },
  };

  for (u32 i=0; i<sizeof(test_cases) / sizeof(test_cases[0]); i++) {
    float res  = costas_discriminator(test_cases[i].i, test_cases[i].q);
    /*
    Result is expected to be close to the predicted value only when
    the test is expected to pass.
    */

    fail_unless(((fabs(test_cases[i].res - res) < 1e-6) && test_cases[i].valid) ||
                ((fabs(test_cases[i].res - res) > 1e-6) && !test_cases[i].valid),
                "result is incorrect (%f vs %f)", test_cases[i].res, res);
  }
}
END_TEST

Suite* track_test_suite(void)
{
  Suite *s = suite_create("Track");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_costas_discriminator);
  suite_add_tcase(s, tc_core);

  return s;
}

