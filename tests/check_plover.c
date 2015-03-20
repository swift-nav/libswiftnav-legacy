#include <check.h>
#include <stdio.h>
#include "amb_kf.h"
#include "single_diff.h"
#include "check_utils.h"

#include "plover.h"

START_TEST(hello_world) {
  plover_test();
}
END_TEST

Suite* plover_test_suite(void)
{
  Suite *s = suite_create("Generated code suite");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, hello_world);
  suite_add_tcase(s, tc_core);

  return s;
}
