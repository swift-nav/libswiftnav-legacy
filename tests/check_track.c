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

START_TEST(test_aided_tl_adjustment)
{
  aided_tl_state_t s;
  aided_tl_init(&s, 50,    /* loop_freq [Hz] */
                0,         /* Doppler on code_freq [Hz] */
                1 /*code_bw*/, 0.7 /*code_zeta*/, 1 /*code_k*/,
                1200 /*carr_to_code*/,
                0 /* Doppler on carr_freq [Hz]*/,
                13 /*carr_bw*/, 0.7 /*carr_zeta*/, 1 /*carr_k*/,
                5 /*carr_freq_b1*/);
  float carr_freq_prev = s.carr_freq;
  float err = 15;
  aided_tl_adjust(&s, err);
  fail_unless((carr_freq_prev + err) == s.carr_freq,
               "incorrect adjusment result");
}
END_TEST

Suite* track_test_suite(void)
{
  Suite *s = suite_create("Track");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_costas_discriminator);
  tcase_add_test(tc_core, test_aided_tl_adjustment);
  suite_add_tcase(s, tc_core);

  return s;
}

