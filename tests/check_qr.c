#include <stdio.h>
#include <stdlib.h>

#include <check.h>

#include "linear_algebra.h"
#include "printing_utils.h"
#include "check_utils.h"
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

double randd() {
  int max = 22;
  return (double)(rand()) / RAND_MAX * max - (double)max / 2;
}

void test_qr(int rows, int cols)
{
  double mat[rows*cols];
  double vec[cols];

  for(int r = 0; r < rows; r++) {
    for(int c = 0; c < cols; c++) {
      mat[r*cols + c] = randd();
    }
  }
  for(int c = 0; c < cols; c++) {
    vec[c] = randd();
  }

  double m1[rows*cols];
  memcpy(m1, mat, sizeof(m1));

  //printf("%d %d\n", rows, cols);
  //print_double_mtx(mat, rows, cols);

  double out[rows];
  matrix_multiply(rows, cols, 1, mat, vec, out);
  double soln[cols];
  double residual;
  int code = qr_solve(rows, cols, mat, out, soln, &residual);
  fail_unless(residual < 0.0000001, "Residual too large: %f\n", residual);
  if (!arr_within_epsilon(cols, soln, vec) && code != -1) {
    printf("original matrix: ");
    print_double_mtx(m1, rows, cols);
    printf("reduced R matrix: ");
    print_double_mtx(mat, rows, cols);
    printf("solution: ");
    print_double_mtx(soln, 1, cols);
    printf("vec: ");
    print_double_mtx(vec, 1, cols);
    fail_unless(false, "Solution wrong\n");
  }
}

/* Performs 10000 trials. */
START_TEST(test_rand)
{
  srand(0);
  for (int i = 0; i < 10000; i++) {
    int rows = rand() % 25 + 1;
    int cols = rand() % rows + 1;
    test_qr(rows, cols);
  }
}
END_TEST

Suite* qr_test_suite(void)
{
  Suite *s = suite_create("Generated QR solver");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_ok);
  tcase_add_test(tc_core, test_bad);
  tcase_add_test(tc_core, test_rand);
  suite_add_tcase(s, tc_core);
  return s;
}
