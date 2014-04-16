
#include <check.h>
#include <stdio.h>
#include <float_kf.h>
#include <amb_kf.h>
#include <linear_algebra.h>
#include <cblas.h>

START_TEST(test_udu)
{
  double M1[2 * 2] = {
    1, 0,
    0, 1
  };
  double U[2 * 2];
  double D[2];
  udu(2, M1, U, D);
  MAT_PRINTF(U,2,2);
  VEC_PRINTF(D,2);

  reconstruct_udu(2, U, D, M1);
  MAT_PRINTF(M1, 2, 2);

  double M2[2 * 2] = {
    1, 2,
    2, 1
  };
  udu(2, M2, U, D);
  MAT_PRINTF((double *) U, 2, 2);
  VEC_PRINTF((double *) D, 2);

  reconstruct_udu(2, U, D, M2);
  MAT_PRINTF((double *) M2, 2, 2);

  double x[2] = {1,3};
  double y[2] = {0,0};
  cblas_dgemv(CblasRowMajor, CblasNoTrans,
              2, 2, 1, (double *) M1,
              2, (double *) x, 1,
              0, y, 1);
  VEC_PRINTF((double *) y, 2);
}
END_TEST

Suite* float_kf_suite(void)
{
  Suite *s = suite_create("Float KF");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_udu);
  suite_add_tcase(s, tc_core);

  return s;
}

