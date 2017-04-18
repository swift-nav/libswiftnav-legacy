#include <check.h>
#include <libswiftnav/track.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define BW   1000
#define CN0_0 40
#define CUTOFF_FREQ 0.1
#define LOOP_FREQ 1000

static s8* generate_I_input(u32 length)
{
  s8* input;
  u32 ii = 0;

  input = malloc(length);
  if (NULL == input) {
    return NULL;
  }

  for (ii = 0; ii < length; ii++) {
    input[ii] = 100;
  }

  return input;
}

static s8* generate_Q_input(u32 length)
{
  s8* input;
  u32 ii = 0;

  input = malloc(length);
  if (NULL == input) {
    return NULL;
  }

  for (ii = 0; ii < length; ii++) {
    input[ii] = 50;
  }

  return input;
}

START_TEST(test_cn0)
{
  cn0_est_state_t s;
  cn0_est_params_t p;
  s8* signal_I;
  s8* signal_Q;
  u32 ii = 0;
  u32 test_length = 1000;
  float cn0 = 0.0;

  signal_I = generate_I_input(test_length);
  fail_if(NULL == signal_I, "Could not allocate I data");
  signal_Q = generate_Q_input(test_length);
  fail_if(NULL == signal_Q, "Could not allocate Q data");

  cn0_est_init(&s, BW, CN0_0, CUTOFF_FREQ, LOOP_FREQ);
  cn0_est_compute_params(&p, BW, CUTOFF_FREQ, LOOP_FREQ);

  for(ii = 0; ii < test_length; ii++) {
	  cn0 = cn0_est(&s, &p, signal_I[ii], signal_Q[ii]);
  }

  fail_if(cn0 < 30.0);

  free(signal_I);
  free(signal_Q);
}
END_TEST

Suite* cn0_suite(void)
{
  Suite *s = suite_create("CN0");
  TCase *tc_core = tcase_create("Core");

  tcase_add_test(tc_core, test_cn0);
  suite_add_tcase(s, tc_core);

  return s;
}
