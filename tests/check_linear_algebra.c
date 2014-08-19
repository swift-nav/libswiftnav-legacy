#include <math.h>
#include <stdlib.h>
#include <check.h>

#include <stdio.h>
#include "check_utils.h"

#include <linear_algebra.h>

#define LINALG_TOL 1e-10
#define LINALG_NUM 22
#define MATRIX_MIN -1e3
#define MATRIX_MAX 1e3
#define mrand frand(MATRIX_MIN, MATRIX_MAX)
#define MSIZE_MAX 64

/* TODO: matrix_multiply, matrix_add_sc, matrix_copy, all vector functions */

START_TEST(test_matrix_inverse_2x2) {
  u32 i, j, t;
  double A[4];
  double B[4];
  double I[4];

  seed_rng();
  /* 2x2 inverses */
  for (t = 0; t < LINALG_NUM; t++) {
    do {
      for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
          A[2*i + j] = mrand;
    } while (matrix_inverse(2, A, B) < 0);
    matrix_multiply(2, 2, 2, A, B, I);
    fail_unless(fabs(I[0] - 1) < LINALG_TOL,
                "Matrix differs from identity: %lf", I[0]);
    fail_unless(fabs(I[3] - 1) < LINALG_TOL,
                "Matrix differs from identity: %lf", I[3]);
  }
  for (i = 0; i < 2; i++)
    for (j = 0; j < 2; j++)
      if (j == 0)
        A[2*i + j] = 22;
      else
        A[2*i + j] = 1;
  s32 mi = matrix_inverse(2, A, B);
  fail_unless(mi < 0, "Singular matrix not detected.");
}
END_TEST

START_TEST(test_matrix_inverse_3x3) {
  u32 i, j, t;
  double A[9];
  double B[9];
  double I[9];

  seed_rng();
  /* 3x3 inverses */
  for (t = 0; t < LINALG_NUM; t++) {
    do {
      for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
          A[3*i + j] = mrand;
    } while (matrix_inverse(3, A, B) < 0);
    matrix_multiply(3, 3, 3, A, B, I);
    fail_unless(fabs(I[0] - 1) < LINALG_TOL,
                "Matrix differs from identity: %lf", I[0]);
    fail_unless(fabs(I[4] - 1) < LINALG_TOL,
                "Matrix differs from identity: %lf", I[4]);
    fail_unless(fabs(I[8] - 1) < LINALG_TOL,
                "Matrix differs from identity: %lf", I[8]);
  }
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      if (j == 0)
        A[3*i + j] = 33;
      else
        A[3*i + j] = 1;
  s32 mi = matrix_inverse(3, A, B);
  fail_unless(mi < 0, "Singular matrix not detected.");
}
END_TEST

START_TEST(test_matrix_inverse_4x4) {
  u32 i, j, t;
  double A[16];
  double B[16];
  double I[16];

  seed_rng();
  /* 4x4 inverses */
  for (t = 0; t < LINALG_NUM; t++) {
    do {
      for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
          A[4*i + j] = mrand;
    } while (matrix_inverse(4, A, B) < 0);
    matrix_multiply(4, 4, 4, A, B, I);
    fail_unless(fabs(I[0] - 1) < LINALG_TOL,
                "Matrix differs from identity: %lf", I[0]);
    fail_unless(fabs(I[5] - 1) < LINALG_TOL,
                "Matrix differs from identity: %lf", I[5]);
    fail_unless(fabs(I[10] - 1) < LINALG_TOL,
                "Matrix differs from identity: %lf", I[10]);
    fail_unless(fabs(I[15] - 1) < LINALG_TOL,
                "Matrix differs from identity: %lf", I[15]);
  }
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      if (j == 0)
        A[4*i + j] = 44;
      else
        A[4*i + j] = 1;
  s32 mi = matrix_inverse(4, A, B);
  fail_unless(mi < 0, "Singular matrix not detected.");
}
END_TEST

START_TEST(test_matrix_inverse_5x5) {
  u32 i, j, t;
  double A[25];
  double B[25];
  double I[25];

  seed_rng();
  /* 5x5 inverses */
  for (t = 0; t < LINALG_NUM; t++) {
    do {
      for (i = 0; i < 5; i++)
        for (j = 0; j < 5; j++)
          A[5*i + j] = mrand;
    } while (matrix_inverse(5, A, B) < 0);
    matrix_multiply(5, 5, 5, A, B, I);
    fail_unless(fabs(I[0] - 1) < LINALG_TOL,
                "Matrix differs from identity: %lf", I[0]);
    fail_unless(fabs(I[6] - 1) < LINALG_TOL,
                "Matrix differs from identity: %lf", I[6]);
    fail_unless(fabs(I[12] - 1) < LINALG_TOL,
                "Matrix differs from identity: %lf", I[12]);
    fail_unless(fabs(I[18] - 1) < LINALG_TOL,
                "Matrix differs from identity: %lf", I[18]);
    fail_unless(fabs(I[24] - 1) < LINALG_TOL,
                "Matrix differs from identity: %lf", I[24]);
  }
  for (i = 0; i < 5; i++)
    for (j = 0; j < 5; j++)
      if (j == 0)
        A[5*i + j] = 55;
      else
        A[5*i + j] = 1;
  s32 mi = matrix_inverse(5, A, B);
  fail_unless(mi < 0, "Singular matrix not detected.");
}
END_TEST

START_TEST(test_matrix_eye)
{
  double M[10][10];

  matrix_eye(10, (double *)M);

  for (u32 i=0; i<10; i++) {
    for (u32 j=0; j<10; j++) {
      if (i == j) {
        fail_unless(M[i][j] == 1, "Identity diagonal element != 1");
      } else {
        fail_unless(M[i][j] == 0, "Identity off-diagonal element != 0");
      }
    }
  }
}
END_TEST

START_TEST(test_matrix_triu)
{
  double M[4][4] = {
    { 1,  2,  3,  4},
    { 5,  6,  7,  8},
    { 9, 10, 11, 12},
    {13, 14, 15, 16}
  };

  double M_[4][4] = {
    {1, 2,  3,  4},
    {0, 6,  7,  8},
    {0, 0, 11, 12},
    {0, 0,  0, 16}
  };

  matrix_triu(4, (double *)M);

  for (u32 i=0; i<4; i++) {
    for (u32 j=0; j<4; j++) {
      fail_unless(M[i][j] == M_[i][j], "triu result != test matrix");
    }
  }
}
END_TEST

START_TEST(test_matrix_udu_1)
{
  double M[4][4] = {
    {100, 145, 121,  16},
    {145, 221, 183,  24},
    {121, 183, 199,  28},
    { 16,  24,  28,   4}
  };

  double U[4][4] = {{0}};
  double D[4] = {0};

  matrix_udu(4, (double *)M, (double *)U, D);

  double U_[4][4] = {
    {1, 2, 3, 4},
    {0, 1, 5, 6},
    {0, 0, 1, 7},
    {0, 0, 0, 1}
  };

  double D_[4] = {1, 2, 3, 4};

  for (u32 i=0; i<4; i++) {
    for (u32 j=0; j<4; j++) {
      fail_unless(U[i][j] == U_[i][j], "U result != test matrix");
    }
  }
  for (u32 i=0; i<4; i++) {
    fail_unless(D[i] == D_[i], "D result != test D");
  }
}
END_TEST

START_TEST(test_matrix_udu_2)
{
  u32 n = sizerand(MSIZE_MAX);
  double M[n][n];
  double M_orig[n][n];

  for (u32 i=0; i<n; i++) {
    for (u32 j=0; j<=i; j++) {
      M[i][j] = M[j][i] = mrand;
    }
  }
  memcpy(M_orig, M, n * n * sizeof(double));

  double U[n][n];
  memset(U, 0, n * n * sizeof(double));
  double D[n];
  memset(D, 0, n * sizeof(double));

  matrix_udu(n, (double *)M, (double *)U, D);

  /* Check U is unit upper triangular. */
  for (u32 i=0; i<n; i++) {
    for (u32 j=0; j<n; j++) {
      if (i == j) {
        fail_unless(fabs(U[i][j] - 1) < LINALG_TOL,
          "U diagonal element != 1 (was %f)", M[i][j]);
      }
      if (i > j) {
        fail_unless(fabs(U[i][j]) < LINALG_TOL, "U lower triangle element != 0");
      }
    }
  }

  /* Check reconstructed matrix is correct. */
  double M_[n][n];
  memset(M_, 0, n * n * sizeof(double));
  matrix_reconstruct_udu(n, (double *)U, D, (double *)M_);

  for (u32 i=0; i<n; i++) {
    for (u32 j=0; j<n; j++) {
      fail_unless(fabs(M_orig[i][j] - M_[i][j]) < LINALG_TOL * MATRIX_MAX,
        "reconstructed result != original matrix, delta[%d][%d] = %f",
        i, j, fabs(M_orig[i][j] - M_[i][j]));
    }
  }
}
END_TEST

START_TEST(test_matrix_udu_3)
{
  double M[3][3] = {
    {36, 49,  9},
    {49, 77, 15},
    { 9, 15,  3},
  };

  double U[3][3] = {{0}};
  double D[3] = {0};

  matrix_udu(3, (double *)M, (double *)U, D);

  /* Check using formula for 3x3 matrix on Gibbs p. 393 */
  fail_unless(D[2] == M[2][2], "D[2] incorrect");
  fail_unless(U[2][1] == M[2][1] / D[2], "U[2][1] incorrect");
  fail_unless(U[2][0] == M[2][0] / D[2], "U[2][0] incorrect");
  fail_unless(D[1] == (M[1][1] - U[2][1]*U[2][1]*D[2]), "D[1] incorrect");
  fail_unless(U[1][0] == ((M[1][0] - U[2][0]*U[2][1]*D[2])/D[1]),
    "U[1][0] incorrect");
  fail_unless(D[0] == (M[0][0] - U[1][0]*U[1][0]*D[1] - U[2][0]*U[2][0]*D[2]),
    "D[0] incorrect");
}
END_TEST

START_TEST(test_matrix_reconstruct_udu)
{
  double U[4][4] = {
    {1, 2, 3, 4},
    {0, 1, 5, 6},
    {0, 0, 1, 7},
    {0, 0, 0, 1}
  };

  double D[4] = {1, 2, 3, 4};

  double M[4][4] = {{0}};

  double M_[4][4] = {
    {100, 145, 121,  16},
    {145, 221, 183,  24},
    {121, 183, 199,  28},
    { 16,  24,  28,   4}
  };

  matrix_reconstruct_udu(4, (double *)U, D, (double *)M);

  for (u32 i=0; i<4; i++) {
    for (u32 j=0; j<4; j++) {
      fail_unless(M[i][j] == M_[i][j], "reconstructed result != test matrix");
    }
  }
}
END_TEST

START_TEST(test_matrix_add_sc) {
  u32 i, j, t;

  seed_rng();
  for (t = 0; t < LINALG_NUM; t++) {
    u32 n = sizerand(MSIZE_MAX);
    u32 m = sizerand(MSIZE_MAX);
    double A[n*m];
    double B[m*n];
    for (i = 0; i < n; i++)
      for (j = 0; j < m; j++) {
        A[m*i + j] = mrand;
      }
    matrix_add_sc(n, m, A, A, -1, B);
    for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
        fail_unless(fabs(B[m*i + j]) < LINALG_TOL,
                    "Matrix differs from zero: %lf", B[m*i + j]);
  }
}
END_TEST

START_TEST(test_matrix_copy) {
  u32 i, j, t;
  double tmp;

  seed_rng();
  for (t = 0; t < LINALG_NUM; t++) {
    u32 n = sizerand(MSIZE_MAX);
    u32 m = sizerand(MSIZE_MAX);
    double A[n*m];
    double B[m*n];
    for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
        A[m*i + j] = mrand;
    matrix_copy(n, m, A, B);
    for (i = 0; i < n; i++)
      for (j = 0; j < m; j++) {
        tmp = fabs(B[m*i + j] - A[m*i + j]);
        fail_unless(tmp < LINALG_TOL,
                    "Matrix differs from zero: %lf", tmp);
      }
  }
}
END_TEST

START_TEST(test_matrix_transpose) {
  u32 i, j, t;

  seed_rng();
  for (t = 0; t < LINALG_NUM; t++) {
    u32 n = sizerand(MSIZE_MAX);
    u32 m = sizerand(MSIZE_MAX);
    double A[n*m];
    double B[m*n];
    double C[n*m];
    for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
        A[m*i + j] = mrand;
    matrix_transpose(n, m, A, B);
    matrix_transpose(m, n, B, C);
    for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
        fail_unless(fabs(A[m*i + j] - C[m*i + j]) < LINALG_TOL,
                    "Matrix element differs from original: %lf, %lf",
                    A[m*i + j], C[m*i + j]);
  }
}
END_TEST

START_TEST(test_vector_dot) {
  u32 i, t;

  seed_rng();
  for (t = 0; t < LINALG_NUM; t++) {
    u32 n = sizerand(MSIZE_MAX);
    u32 mid;
    if (n % 2 == 0)
      mid = n / 2;
    else
      mid = (n - 1) / 2;

    double A[n], B[n];
    for (i = 0; i < n; i++) {
      A[i] = mrand/1e20;
      if (i < mid)
        B[n-i-1] = -A[i];
      else
        B[n-i-1] = A[i];
    }
    double dot = vector_dot(n, A, B);
    if (n % 2 == 0)
      fail_unless(fabs(dot) < LINALG_TOL,
                  "Dot product differs from zero: %lf", vector_dot(n, A, B));
    else
      fail_unless(fabs(dot - A[mid]*B[mid]) < LINALG_TOL,
                  "Dot product differs from square of middle element "
                  "%lf: %lf (%lf)",
                  A[mid]*B[mid], dot, dot - A[mid]*B[mid]);
  }
}
END_TEST

START_TEST(test_vector_mean) {
  u32 i, t;

  seed_rng();
  for (t = 0; t < LINALG_NUM; t++) {
    u32 n = sizerand(MSIZE_MAX);
    double A[n];
    double test = mrand/1e22;
    for (i = 0; i < n; i++)
      A[i] = test + i;
    double mean = vector_mean(n, A);
    double expect = test + (n - 1.0) / 2.0;
    fail_unless(fabs(mean - expect)
                < LINALG_TOL,
                "Mean differs from expected %lf: %lf (%lf)",
                expect, mean, fabs(mean - expect));
  }
}
END_TEST

START_TEST(test_vector_norm) {
  u32 i, t;

  seed_rng();
  for (t = 0; t < LINALG_NUM; t++) {
    u32 n = sizerand(MSIZE_MAX);
    double test = mrand/1e22;
    double A[n];
    for (i = 0; i < n; i++)
      A[i] = test;
    fail_unless(fabs(vector_norm(n, A)*vector_norm(n, A) - n * test * test)
                < LINALG_TOL * vector_norm(n, A),
                "Norm differs from expected %lf: %lf (%lf)",
                n * test * test, vector_norm(n, A) * vector_norm(n, A),
                fabs(vector_norm(n, A) * vector_norm(n, A) - n * test * test));
  }
}
END_TEST

START_TEST(test_vector_normalize) {
  u32 i, t;

  seed_rng();
  for (t = 0; t < LINALG_NUM; t++) {
    u32 n = sizerand(MSIZE_MAX);
    double A[n];
    for (i = 0; i < n; i++)
      A[i] = mrand;
    vector_normalize(n, A);
    double vnorm = vector_norm(n, A);
    fail_unless(fabs(vnorm - 1) < LINALG_TOL,
                "Norm differs from 1: %lf",
                vector_norm(n, A));
  }
}
END_TEST

START_TEST(test_vector_add_sc) {
  u32 i, t;

  seed_rng();
  for (t = 0; t < LINALG_NUM; t++) {
    u32 n = sizerand(MSIZE_MAX);
    double A[n], B[n];
    for (i = 0; i < n; i++)
      A[i] = mrand;
    vector_add_sc(n, A, A, -1, B);
    for (i = 0; i < n; i++)
      fail_unless(fabs(B[i]) < LINALG_TOL,
                  "Vector element differs from 0: %lf",
                  B[i]);
  }
}
END_TEST

START_TEST(test_vector_add) {
  u32 i, t;

  seed_rng();
  for (t = 0; t < LINALG_NUM; t++) {
    u32 n = sizerand(MSIZE_MAX);
    double A[n], B[n], C[n];
    for (i = 0; i < n; i++) {
      A[i] = mrand;
      B[i] = -A[i];
    }
    vector_add(n, A, B, C);
    for (i = 0; i < n; i++)
      fail_unless(fabs(C[i]) < LINALG_TOL,
                  "Vector element differs from 0: %lf",
                  C[i]);
  }
}
END_TEST

START_TEST(test_vector_subtract) {
  u32 i, t;

  seed_rng();
  for (t = 0; t < LINALG_NUM; t++) {
    u32 n = sizerand(MSIZE_MAX);
    double A[n], B[n], C[n];
    for (i = 0; i < n; i++) {
      A[i] = mrand;
      B[i] = A[i];
    }
    vector_subtract(n, A, B, C);
    for (i = 0; i < n; i++)
      fail_unless(fabs(C[i]) < LINALG_TOL,
                  "Vector element differs from 0: %lf",
                  C[i]);
  }
}
END_TEST

START_TEST(test_vector_cross) {
  u32 i, t;

  seed_rng();
  for (t = 0; t < LINALG_NUM; t++) {
    double A[3], B[3], C[3], D[3];
    for (i = 0; i < 3; i++) {
      A[i] = mrand;
      B[i] = A[i];
    }
    vector_cross(A, B, C);
    for (i = 0; i < 3; i++)
      fail_unless(fabs(C[i]) < LINALG_TOL,
                  "Vector element differs from 0: %lf",
                  C[i]);
    for (i = 0; i < 3; i++) {
      A[i] = mrand;
      B[i] = mrand;
    }
    vector_cross(A, B, C);
    for (i = 0; i < 3; i++) B[i] *= -1;
    vector_cross(B, A, D);
    for (i = 0; i < 3; i++)
      fail_unless(fabs(C[i] - D[i]) < LINALG_TOL,
                  "Vector equality fails: %lf != %lf",
                  C[i], D[i]);
  }
}
END_TEST

START_TEST(test_vector_three) {
  u32 i, t;

  seed_rng();
  for (t = 0; t < LINALG_NUM; t++) {
    double A[3], B[3], C[3], tmp[3];
    double D, E, F, norm;
    for (i = 0; i < 3; i++) {
      A[i] = mrand / 1e20;
      B[i] = mrand / 1e20;
      C[i] = mrand / 1e20;
    }
    /* Check triple product identity */
    vector_cross(A, B, tmp);
    D = vector_dot(3, C, tmp);
    vector_cross(B, C, tmp);
    E = vector_dot(3, A, tmp);
    vector_cross(C, A, tmp);
    F = vector_dot(3, B, tmp);

    norm = (vector_norm(3, A) + vector_norm(3, B) + vector_norm(3, C))/3;
    fail_unless(fabs(E - D) < LINALG_TOL * norm,
                "Triple product failure between %lf and %lf",
                D, E);
    fail_unless(fabs(E - F) < LINALG_TOL * norm,
                "Triple product failure between %lf and %lf",
                E, F);
    fail_unless(fabs(F - D) < LINALG_TOL * norm,
                "Triple product failure between %lf and %lf",
                F, D);
  }
}
END_TEST

/*
START_TEST(test_qrsolve_consistency) {
  u32 i, j, t;
  double norm;

  seed_rng();
  for (t = 0; t < LINALG_NUM; t++) {
    u32 n = sizerand(MSIZE_MAX);
    double x_gauss[n], x_qr[n];
    double A[n*n], Ainv[n*n], b[n];
    do {
      for (i = 0; i < n; i++) {
        b[i] = mrand;
        for (j = 0; j < n; j++)
          A[n*i + j] = mrand;
      }
    } while (matrix_inverse(n, A, Ainv) < 0);
    matrix_multiply(n, n, 1, Ainv, b, x_gauss);
    qrsolve(A, n, n, b, x_qr);

    norm = (vector_norm(n, x_qr) + vector_norm(n, x_gauss)) / 2;
    for (i = 0; i < n; i++)
      fail_unless(fabs(x_qr[i] - x_gauss[i]) < LINALG_TOL * norm,
                  "QR solve failure; difference was %lf for element %u",
                  x_qr[i] - x_gauss[i], i);
  }
}
END_TEST

START_TEST(test_qrsolve_rect) {
  s32 i;
  const double A[8] = {-0.0178505395610981, 1.4638781031761146,
                       -0.8242742209580581, -0.6843477128009663,
                       0.9155272861151404, -0.1651159277864960,
                       -0.9929037180867774, -0.1491537478964264};
  double Q[16], R[8];

  double buf[10] __attribute__((unused)) = {22, 22, 22, 22, 22,
                                            22, 22, 22, 22, 22};

  i = qrdecomp(A, 4, 2, Q, R);

  printf("i returned %d\n", i);

  MAT_PRINTF(A, 4, 2);
  MAT_PRINTF(Q, 4, 4);
  MAT_PRINTF(R, 4, 2);
}
END_TEST
*/

Suite* linear_algebra_suite(void) {
  Suite *s = suite_create("Linear algebra");

  /* Core test case */
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_matrix_eye);
  tcase_add_test(tc_core, test_matrix_triu);
  tcase_add_test(tc_core, test_matrix_reconstruct_udu);
  tcase_add_test(tc_core, test_matrix_udu_1);
  tcase_add_test(tc_core, test_matrix_udu_2);
  tcase_add_test(tc_core, test_matrix_udu_3);
  tcase_add_test(tc_core, test_matrix_add_sc);
  tcase_add_test(tc_core, test_matrix_copy);
  tcase_add_test(tc_core, test_matrix_transpose);
  tcase_add_test(tc_core, test_matrix_inverse_2x2);
  tcase_add_test(tc_core, test_matrix_inverse_3x3);
  tcase_add_test(tc_core, test_matrix_inverse_4x4);
  tcase_add_test(tc_core, test_matrix_inverse_5x5);

  tcase_add_test(tc_core, test_vector_dot);
  tcase_add_test(tc_core, test_vector_mean);
  tcase_add_test(tc_core, test_vector_norm);
  tcase_add_test(tc_core, test_vector_normalize);
  tcase_add_test(tc_core, test_vector_add_sc);
  tcase_add_test(tc_core, test_vector_add);
  tcase_add_test(tc_core, test_vector_subtract);
  tcase_add_test(tc_core, test_vector_cross);
  tcase_add_test(tc_core, test_vector_three);
  /*tcase_add_test(tc_core, test_qrsolve_consistency);*/
  /*tcase_add_test(tc_core, test_qrsolve_rect);*/
  suite_add_tcase(s, tc_core);

  return s;
}

