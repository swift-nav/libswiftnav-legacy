#include <stdio.h>

#include <check.h>

#include "linear_algebra.h"
#include "plover/qr.h"


START_TEST(test_ok)
{
  double m1[] = {1, 0, 2,
                 1, 1, -1,
                 0, 0, 1,
                 0, 0, 22 };
  double b[] = {5, 4, 1.1, 22};
  double soln[3];
  double residual;

  s8 code = qr_solve(4, 3, m1, b, soln, &residual);
  
  fail_unless(code == 0, "Solver error code: %d\n", code);
  fail_unless(residual < 0.1, "Residual too large: %f\n", residual);
}
END_TEST

/* Deficient */
START_TEST(test_bad)
{
  double m1[] = {1, 0, 0,
                 0, 1, 0,
                 0, 0, 0 };
  double b[] = {1,1,1};
  double soln[3];
  double residual;

  s8 code = qr_solve(3, 3, m1, b, soln, &residual);
  
  fail_unless(code == -1, "Wrong solver code: %d\n", code);
}
END_TEST

Suite* qr_test_suite(void)
{
  Suite *s = suite_create("Generated QR solver");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_ok);
  tcase_add_test(tc_core, test_bad);
  suite_add_tcase(s, tc_core);
  return s;
}
