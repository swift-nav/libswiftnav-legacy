#include "plover/prelude.h"


int ipow(int base, int exp) {
  int result = 1;
  while (exp) {
    if (exp & 1)  result *= base;
    exp >>= 1;
    base *= base;
  }
  return result;
}
double dipow(double base, int exp) {
  double result = 1;
  while (exp) {
    if (exp & 1)  result *= base;
    exp >>= 1;
    base *= base;
  }
  return result;
}
double rand_uniform(void) {
  union { double d; u64 i; } di;
  di.i = 0;
  const int width = 32 - __builtin_clz((u32)RAND_MAX);
  for (int j = 0; j < DBL_MANT_DIG; j += width) {
    di.i = (di.i << width) | rand();
  }
  di.i = (di.i << (64 - DBL_MANT_DIG)) >> (64 - DBL_MANT_DIG);
  di.i |= 0x3ff0000000000000; // zero exponent
  return di.d-1;
}
#define MATRIX_EPSILON (1e-60)
static s32 inv1 (const double * A, double * B);
static double det2 (const double * A);
static s32 inv2 (const double * A, double * B);
static double det3 (const double * A);
static s32 inv3 (const double * A, double * B);
static double det4 (const double * A);
static s32 inv4 (const double * A, double * B);
static void givens (const double a, const double b, double * result);
static double gen_det_qr (const s32 n, const double * U);
static s32 gen_inv_qr (const s32 n, const double * U, double * B);
double rand_normal (void)
{
    double x1;
    double x2;
    double w;
    
    do {
        x1 = 2 * rand_uniform() - 1.0;
        x2 = 2 * rand_uniform() - 1.0;
        w = dipow(x1, 2) + dipow(x2, 2);
    } while(1.0 <= w);
    w = sqrt(-(2 * log(w) / w));
    return x1 * w;
}
double norm (const s32 n, const double * v)
{
    double sum = 0;
    
    for (s32 idx = 0; idx < n; idx++) {
        sum += v[idx] * v[idx];
    }
    return sqrt(sum);
}
void normalize (const s32 n, const double * v, double * result)
{
    double tmp;
    
    tmp = norm(n, v);
    for (s32 idx = 0; idx < n; idx++) {
        result[idx] = v[idx] / tmp;
    }
}
void print_vec (const s32 n, const double * v)
{
    printf("vec(");
    for (s32 i = 0, idx = 0; idx < n; i += 1, idx++) {
        if (0 < i) {
            printf(",");
        }
        printf("% 12lf", v[i]);
    }
    printf(")\n");
}
void print_mat (const s32 n, const s32 m, const double * A)
{
    printf("mat(");
    for (s32 i = 0, idx = 0; idx < n; i += 1, idx++) {
        if (0 < i) {
            printf(";\n    ");
        }
        for (s32 j = 0, idx2 = 0; idx2 < m; j += 1, idx2++) {
            if (0 < j) {
                printf(",");
            }
            printf("% 12lf", A[m * i + j]);
        }
    }
    printf(")\n");
}
s32 matrix_inv (const s32 n, const double * A, double * B)
{
    s32 tmp;
    
    switch (n) {
        
        
      case 0:
        {
            tmp = 0;
            break;
        }
        
        
      case 1:
        {
            tmp = inv1(A, B);
            break;
        }
        
        
      case 2:
        {
            tmp = inv2(A, B);
            break;
        }
        
        
      case 3:
        {
            tmp = inv3(A, B);
            break;
        }
        
        
      case 4:
        {
            tmp = inv4(A, B);
            break;
        }
        
      default:
        {
            tmp = gen_inv_qr(n, A, B);
        }
    }
    return tmp;
}
double det (const s32 n, const double * A)
{
    double tmp;
    
    switch (n) {
        
        
      case 0:
        {
            tmp = 1;
            break;
        }
        
        
      case 1:
        {
            tmp = A[n * 0];
            break;
        }
        
        
      case 2:
        {
            tmp = det2(A);
            break;
        }
        
        
      case 3:
        {
            tmp = det3(A);
            break;
        }
        
        
      case 4:
        {
            tmp = det4(A);
            break;
        }
        
      default:
        {
            tmp = gen_det_qr(n, A);
        }
    }
    return tmp;
}
s32 inv1 (const double * A, double * B)
{
    if (fabs(A[0]) < MATRIX_EPSILON) {
        return -1;
    }
    B[0] = 1 / A[0];
    return 0;
}
double det2 (const double * A)
{
    return A[2 * 0] * A[2 * 1 + 1] - A[2 * 0 + 1] * A[2 * 1];
}
s32 inv2 (const double * A, double * B)
{
    double d;
    
    d = det2(A);
    if (fabs(d) < MATRIX_EPSILON) {
        return -1;
    }
    B[2 * 0] = A[2 * 1 + 1] / d;
    B[2 * 0 + 1] = -(A[2 * 0 + 1] / d);
    B[2 * 1] = -(A[2 * 1] / d);
    B[2 * 1 + 1] = A[2 * 0] / d;
    return 0;
}
double det3 (const double * A)
{
    return -(A[3 * 1] * (A[3 * 0 + 1] * A[3 * 2 + 2] - A[3 * 0 + 2] * A[3 * 2 + 1])) + A[3 * 1 + 1] * (A[3 * 0] * A[3 *
                                                                                                                    2 +
                                                                                                                    2] -
                                                                                                       A[3 * 0 + 2] *
                                                                                                       A[3 * 2]) - A[3 *
                                                                                                                     1 +
                                                                                                                     2] *
        (A[3 * 0] * A[3 * 2 + 1] - A[3 * 0 + 1] * A[3 * 2]);
}
s32 inv3 (const double * A, double * B)
{
    double d;
    
    d = det3(A);
    if (fabs(d) < MATRIX_EPSILON) {
        return -1;
    }
    B[3 * 0] = (A[3 * 1 + 1] * A[3 * 2 + 2] - A[3 * 1 + 2] * A[3 * 2 + 1]) / d;
    B[3 * 1] = -((A[3 * 1] * A[3 * 2 + 2] - A[3 * 1 + 2] * A[3 * 2]) / d);
    B[3 * 2] = (A[3 * 1] * A[3 * 2 + 1] - A[3 * 1 + 1] * A[3 * 2]) / d;
    B[3 * 0 + 1] = -((A[3 * 0 + 1] * A[3 * 2 + 2] - A[3 * 0 + 2] * A[3 * 2 + 1]) / d);
    B[3 * 1 + 1] = (A[3 * 0] * A[3 * 2 + 2] - A[3 * 0 + 2] * A[3 * 2]) / d;
    B[3 * 2 + 1] = -((A[3 * 0] * A[3 * 2 + 1] - A[3 * 0 + 1] * A[3 * 2]) / d);
    B[3 * 0 + 2] = (A[3 * 0 + 1] * A[3 * 1 + 2] - A[3 * 0 + 2] * A[3 * 1 + 1]) / d;
    B[3 * 1 + 2] = -((A[3 * 0] * A[3 * 1 + 2] - A[3 * 0 + 2] * A[3 * 1]) / d);
    B[3 * 2 + 2] = (A[3 * 0] * A[3 * 1 + 1] - A[3 * 0 + 1] * A[3 * 1]) / d;
    return 0;
}
double det4 (const double * A)
{
    return A[4 * 1] * (A[4 * 2 + 1] * (A[4 * 0 + 2] * A[4 * 3 + 3] - A[4 * 0 + 3] * A[4 * 3 + 2]) - A[4 * 2 + 2] *
                       (A[4 * 0 + 1] * A[4 * 3 + 3] - A[4 * 0 + 3] * A[4 * 3 + 1]) + A[4 * 2 + 3] * (A[4 * 0 + 1] *
                                                                                                     A[4 * 3 + 2] -
                                                                                                     A[4 * 0 + 2] *
                                                                                                     A[4 * 3 + 1])) -
        A[4 * 1 + 1] * (A[4 * 2] * (A[4 * 0 + 2] * A[4 * 3 + 3] - A[4 * 0 + 3] * A[4 * 3 + 2]) - A[4 * 2 + 2] * (A[4 *
                                                                                                                   0] *
                                                                                                                 A[4 *
                                                                                                                   3 +
                                                                                                                   3] -
                                                                                                                 A[4 *
                                                                                                                   0 +
                                                                                                                   3] *
                                                                                                                 A[4 *
                                                                                                                   3]) +
                        A[4 * 2 + 3] * (A[4 * 0] * A[4 * 3 + 2] - A[4 * 0 + 2] * A[4 * 3])) + A[4 * 1 + 2] * (A[4 * 2] *
                                                                                                              (A[4 * 0 +
                                                                                                                 1] *
                                                                                                               A[4 * 3 +
                                                                                                                 3] -
                                                                                                               A[4 * 0 +
                                                                                                                 3] *
                                                                                                               A[4 * 3 +
                                                                                                                 1]) -
                                                                                                              A[4 * 2 +
                                                                                                                1] *
                                                                                                              (A[4 *
                                                                                                                 0] *
                                                                                                               A[4 * 3 +
                                                                                                                 3] -
                                                                                                               A[4 * 0 +
                                                                                                                 3] *
                                                                                                               A[4 *
                                                                                                                 3]) +
                                                                                                              A[4 * 2 +
                                                                                                                3] *
                                                                                                              (A[4 *
                                                                                                                 0] *
                                                                                                               A[4 * 3 +
                                                                                                                 1] -
                                                                                                               A[4 * 0 +
                                                                                                                 1] *
                                                                                                               A[4 *
                                                                                                                 3])) -
        A[4 * 1 + 3] * (A[4 * 2] * (A[4 * 0 + 1] * A[4 * 3 + 2] - A[4 * 0 + 2] * A[4 * 3 + 1]) - A[4 * 2 + 1] * (A[4 *
                                                                                                                   0] *
                                                                                                                 A[4 *
                                                                                                                   3 +
                                                                                                                   2] -
                                                                                                                 A[4 *
                                                                                                                   0 +
                                                                                                                   2] *
                                                                                                                 A[4 *
                                                                                                                   3]) +
                        A[4 * 2 + 2] * (A[4 * 0] * A[4 * 3 + 1] - A[4 * 0 + 1] * A[4 * 3]));
}
s32 inv4 (const double * A, double * B)
{
    double d;
    
    d = det4(A);
    if (fabs(d) < MATRIX_EPSILON) {
        return -1;
    }
    B[4 * 0] = (-(A[4 * 2 + 1] * (A[4 * 1 + 2] * A[4 * 3 + 3] - A[4 * 1 + 3] * A[4 * 3 + 2])) + A[4 * 2 + 2] * (A[4 *
                                                                                                                  1 +
                                                                                                                  1] *
                                                                                                                A[4 *
                                                                                                                  3 +
                                                                                                                  3] -
                                                                                                                A[4 *
                                                                                                                  1 +
                                                                                                                  3] *
                                                                                                                A[4 *
                                                                                                                  3 +
                                                                                                                  1]) -
                A[4 * 2 + 3] * (A[4 * 1 + 1] * A[4 * 3 + 2] - A[4 * 1 + 2] * A[4 * 3 + 1])) / d;
    B[4 * 1] = -((-(A[4 * 2] * (A[4 * 1 + 2] * A[4 * 3 + 3] - A[4 * 1 + 3] * A[4 * 3 + 2])) + A[4 * 2 + 2] * (A[4 * 1] *
                                                                                                              A[4 * 3 +
                                                                                                                3] -
                                                                                                              A[4 * 1 +
                                                                                                                3] *
                                                                                                              A[4 *
                                                                                                                3]) -
                  A[4 * 2 + 3] * (A[4 * 1] * A[4 * 3 + 2] - A[4 * 1 + 2] * A[4 * 3])) / d);
    B[4 * 2] = (-(A[4 * 2] * (A[4 * 1 + 1] * A[4 * 3 + 3] - A[4 * 1 + 3] * A[4 * 3 + 1])) + A[4 * 2 + 1] * (A[4 * 1] *
                                                                                                            A[4 * 3 +
                                                                                                              3] - A[4 *
                                                                                                                     1 +
                                                                                                                     3] *
                                                                                                            A[4 * 3]) -
                A[4 * 2 + 3] * (A[4 * 1] * A[4 * 3 + 1] - A[4 * 1 + 1] * A[4 * 3])) / d;
    B[4 * 3] = -((-(A[4 * 2] * (A[4 * 1 + 1] * A[4 * 3 + 2] - A[4 * 1 + 2] * A[4 * 3 + 1])) + A[4 * 2 + 1] * (A[4 * 1] *
                                                                                                              A[4 * 3 +
                                                                                                                2] -
                                                                                                              A[4 * 1 +
                                                                                                                2] *
                                                                                                              A[4 *
                                                                                                                3]) -
                  A[4 * 2 + 2] * (A[4 * 1] * A[4 * 3 + 1] - A[4 * 1 + 1] * A[4 * 3])) / d);
    B[4 * 0 + 1] = -((-(A[4 * 2 + 1] * (A[4 * 0 + 2] * A[4 * 3 + 3] - A[4 * 0 + 3] * A[4 * 3 + 2])) + A[4 * 2 + 2] *
                      (A[4 * 0 + 1] * A[4 * 3 + 3] - A[4 * 0 + 3] * A[4 * 3 + 1]) - A[4 * 2 + 3] * (A[4 * 0 + 1] * A[4 *
                                                                                                                     3 +
                                                                                                                     2] -
                                                                                                    A[4 * 0 + 2] * A[4 *
                                                                                                                     3 +
                                                                                                                     1])) /
                     d);
    B[4 * 1 + 1] = (-(A[4 * 2] * (A[4 * 0 + 2] * A[4 * 3 + 3] - A[4 * 0 + 3] * A[4 * 3 + 2])) + A[4 * 2 + 2] * (A[4 *
                                                                                                                  0] *
                                                                                                                A[4 *
                                                                                                                  3 +
                                                                                                                  3] -
                                                                                                                A[4 *
                                                                                                                  0 +
                                                                                                                  3] *
                                                                                                                A[4 *
                                                                                                                  3]) -
                    A[4 * 2 + 3] * (A[4 * 0] * A[4 * 3 + 2] - A[4 * 0 + 2] * A[4 * 3])) / d;
    B[4 * 2 + 1] = -((-(A[4 * 2] * (A[4 * 0 + 1] * A[4 * 3 + 3] - A[4 * 0 + 3] * A[4 * 3 + 1])) + A[4 * 2 + 1] * (A[4 *
                                                                                                                    0] *
                                                                                                                  A[4 *
                                                                                                                    3 +
                                                                                                                    3] -
                                                                                                                  A[4 *
                                                                                                                    0 +
                                                                                                                    3] *
                                                                                                                  A[4 *
                                                                                                                    3]) -
                      A[4 * 2 + 3] * (A[4 * 0] * A[4 * 3 + 1] - A[4 * 0 + 1] * A[4 * 3])) / d);
    B[4 * 3 + 1] = (-(A[4 * 2] * (A[4 * 0 + 1] * A[4 * 3 + 2] - A[4 * 0 + 2] * A[4 * 3 + 1])) + A[4 * 2 + 1] * (A[4 *
                                                                                                                  0] *
                                                                                                                A[4 *
                                                                                                                  3 +
                                                                                                                  2] -
                                                                                                                A[4 *
                                                                                                                  0 +
                                                                                                                  2] *
                                                                                                                A[4 *
                                                                                                                  3]) -
                    A[4 * 2 + 2] * (A[4 * 0] * A[4 * 3 + 1] - A[4 * 0 + 1] * A[4 * 3])) / d;
    B[4 * 0 + 2] = (-(A[4 * 1 + 1] * (A[4 * 0 + 2] * A[4 * 3 + 3] - A[4 * 0 + 3] * A[4 * 3 + 2])) + A[4 * 1 + 2] *
                    (A[4 * 0 + 1] * A[4 * 3 + 3] - A[4 * 0 + 3] * A[4 * 3 + 1]) - A[4 * 1 + 3] * (A[4 * 0 + 1] * A[4 *
                                                                                                                   3 +
                                                                                                                   2] -
                                                                                                  A[4 * 0 + 2] * A[4 *
                                                                                                                   3 +
                                                                                                                   1])) /
        d;
    B[4 * 1 + 2] = -((-(A[4 * 1] * (A[4 * 0 + 2] * A[4 * 3 + 3] - A[4 * 0 + 3] * A[4 * 3 + 2])) + A[4 * 1 + 2] * (A[4 *
                                                                                                                    0] *
                                                                                                                  A[4 *
                                                                                                                    3 +
                                                                                                                    3] -
                                                                                                                  A[4 *
                                                                                                                    0 +
                                                                                                                    3] *
                                                                                                                  A[4 *
                                                                                                                    3]) -
                      A[4 * 1 + 3] * (A[4 * 0] * A[4 * 3 + 2] - A[4 * 0 + 2] * A[4 * 3])) / d);
    B[4 * 2 + 2] = (-(A[4 * 1] * (A[4 * 0 + 1] * A[4 * 3 + 3] - A[4 * 0 + 3] * A[4 * 3 + 1])) + A[4 * 1 + 1] * (A[4 *
                                                                                                                  0] *
                                                                                                                A[4 *
                                                                                                                  3 +
                                                                                                                  3] -
                                                                                                                A[4 *
                                                                                                                  0 +
                                                                                                                  3] *
                                                                                                                A[4 *
                                                                                                                  3]) -
                    A[4 * 1 + 3] * (A[4 * 0] * A[4 * 3 + 1] - A[4 * 0 + 1] * A[4 * 3])) / d;
    B[4 * 3 + 2] = -((-(A[4 * 1] * (A[4 * 0 + 1] * A[4 * 3 + 2] - A[4 * 0 + 2] * A[4 * 3 + 1])) + A[4 * 1 + 1] * (A[4 *
                                                                                                                    0] *
                                                                                                                  A[4 *
                                                                                                                    3 +
                                                                                                                    2] -
                                                                                                                  A[4 *
                                                                                                                    0 +
                                                                                                                    2] *
                                                                                                                  A[4 *
                                                                                                                    3]) -
                      A[4 * 1 + 2] * (A[4 * 0] * A[4 * 3 + 1] - A[4 * 0 + 1] * A[4 * 3])) / d);
    B[4 * 0 + 3] = -((-(A[4 * 1 + 1] * (A[4 * 0 + 2] * A[4 * 2 + 3] - A[4 * 0 + 3] * A[4 * 2 + 2])) + A[4 * 1 + 2] *
                      (A[4 * 0 + 1] * A[4 * 2 + 3] - A[4 * 0 + 3] * A[4 * 2 + 1]) - A[4 * 1 + 3] * (A[4 * 0 + 1] * A[4 *
                                                                                                                     2 +
                                                                                                                     2] -
                                                                                                    A[4 * 0 + 2] * A[4 *
                                                                                                                     2 +
                                                                                                                     1])) /
                     d);
    B[4 * 1 + 3] = (-(A[4 * 1] * (A[4 * 0 + 2] * A[4 * 2 + 3] - A[4 * 0 + 3] * A[4 * 2 + 2])) + A[4 * 1 + 2] * (A[4 *
                                                                                                                  0] *
                                                                                                                A[4 *
                                                                                                                  2 +
                                                                                                                  3] -
                                                                                                                A[4 *
                                                                                                                  0 +
                                                                                                                  3] *
                                                                                                                A[4 *
                                                                                                                  2]) -
                    A[4 * 1 + 3] * (A[4 * 0] * A[4 * 2 + 2] - A[4 * 0 + 2] * A[4 * 2])) / d;
    B[4 * 2 + 3] = -((-(A[4 * 1] * (A[4 * 0 + 1] * A[4 * 2 + 3] - A[4 * 0 + 3] * A[4 * 2 + 1])) + A[4 * 1 + 1] * (A[4 *
                                                                                                                    0] *
                                                                                                                  A[4 *
                                                                                                                    2 +
                                                                                                                    3] -
                                                                                                                  A[4 *
                                                                                                                    0 +
                                                                                                                    3] *
                                                                                                                  A[4 *
                                                                                                                    2]) -
                      A[4 * 1 + 3] * (A[4 * 0] * A[4 * 2 + 1] - A[4 * 0 + 1] * A[4 * 2])) / d);
    B[4 * 3 + 3] = (-(A[4 * 1] * (A[4 * 0 + 1] * A[4 * 2 + 2] - A[4 * 0 + 2] * A[4 * 2 + 1])) + A[4 * 1 + 1] * (A[4 *
                                                                                                                  0] *
                                                                                                                A[4 *
                                                                                                                  2 +
                                                                                                                  2] -
                                                                                                                A[4 *
                                                                                                                  0 +
                                                                                                                  2] *
                                                                                                                A[4 *
                                                                                                                  2]) -
                    A[4 * 1 + 2] * (A[4 * 0] * A[4 * 2 + 1] - A[4 * 0 + 1] * A[4 * 2])) / d;
    return 0;
}
void givens (const double a, const double b, double * result)
{
    double c;
    
    c = 0.0;
    
    double s;
    
    s = 0.0;
    if (b == 0) {
        c = 1;
        s = 0;
    } else {
        if (fabs(a) < fabs(b)) {
            double tau;
            
            tau = -(a / b);
            s = 1 / sqrt(1 + tau * tau);
            c = s * tau;
        } else {
            double tau;
            
            tau = -(b / a);
            c = 1 / sqrt(1 + tau * tau);
            s = c * tau;
        }
    }
    result[2 * 0] = c;
    result[2 * 0 + 1] = s;
    result[2 * 1] = -s;
    result[2 * 1 + 1] = c;
}
double gen_det_qr (const s32 n, const double * U)
{
    double A [n * n];
    
    for (s32 idx = 0; idx < n; idx++) {
        for (s32 idx2 = 0; idx2 < n; idx2++) {
            A[n * idx + idx2] = U[n * idx + idx2];
        }
    }
    for (s32 j = 1, idx = 0; idx < n; j += 1, idx++) {
        for (s32 i = n, idx2 = 0; idx2 < -j + n; i += -1, idx2++) {
            double rot [2 * 2];
            
            givens(A[n * (i - 2) + (j - 1)], A[n * (i - 1) + (j - 1)], rot);
            for (s32 k = j, idx3 = 0; idx3 < 1 - j + n; k += 1, idx3++) {
                double v [2];
                
                for (s32 idx4 = 0; idx4 < 2; idx4++) {
                    v[idx4] = A[n * (i - 2 + idx4) + (k - 1)];
                }
                for (s32 idx4 = 0; idx4 < 2; idx4++) {
                    double sum = 0;
                    
                    for (s32 idx5 = 0; idx5 < 2; idx5++) {
                        sum += rot[2 * idx5 + idx4] * v[idx5];
                    }
                    A[n * (i - 2 + idx4) + (k - 1)] = sum;
                }
            }
        }
    }
    
    double d;
    
    d = 1.0;
    for (s32 i = 0, idx = 0; idx < n; i += 1, idx++) {
        d = d * A[n * i + i];
    }
    return d;
}
s32 gen_inv_qr (const s32 n, const double * U, double * B)
{
    double A [n * n];
    
    for (s32 idx = 0; idx < n; idx++) {
        for (s32 idx2 = 0; idx2 < n; idx2++) {
            A[n * idx + idx2] = U[n * idx + idx2];
        }
    }
    for (s32 idx = 0; idx < n; idx++) {
        s32 i = idx;
        
        for (s32 idx2 = 0; idx2 < n; idx2++) {
            s32 j = idx2;
            s32 loc;
            
            if (i == j) {
                loc = 1;
            } else {
                loc = 0;
            }
            B[n * idx + idx2] = loc;
        }
    }
    for (s32 j = 1, idx = 0; idx < n; j += 1, idx++) {
        for (s32 i = n, idx2 = 0; idx2 < -j + n; i += -1, idx2++) {
            double rot [2 * 2];
            
            givens(A[n * (i - 2) + (j - 1)], A[n * (i - 1) + (j - 1)], rot);
            for (s32 k = j, idx3 = 0; idx3 < 1 - j + n; k += 1, idx3++) {
                double v [2];
                
                for (s32 idx4 = 0; idx4 < 2; idx4++) {
                    v[idx4] = A[n * (i - 2 + idx4) + (k - 1)];
                }
                for (s32 idx4 = 0; idx4 < 2; idx4++) {
                    double sum = 0;
                    
                    for (s32 idx5 = 0; idx5 < 2; idx5++) {
                        sum += rot[2 * idx5 + idx4] * v[idx5];
                    }
                    A[n * (i - 2 + idx4) + (k - 1)] = sum;
                }
            }
            for (s32 k = 1, idx3 = 0; idx3 < n; k += 1, idx3++) {
                double w [2];
                
                for (s32 idx4 = 0; idx4 < 2; idx4++) {
                    w[idx4] = B[n * (i - 2 + idx4) + (k - 1)];
                }
                for (s32 idx4 = 0; idx4 < 2; idx4++) {
                    double sum = 0;
                    
                    for (s32 idx5 = 0; idx5 < 2; idx5++) {
                        sum += rot[2 * idx5 + idx4] * w[idx5];
                    }
                    B[n * (i - 2 + idx4) + (k - 1)] = sum;
                }
            }
        }
    }
    for (s32 i = 0, idx = 0; idx < n; i += 1, idx++) {
        if (A[n * i + i] == 0) {
            return -1;
        }
        
        double scale;
        
        scale = A[n * i + i];
        for (s32 idx2 = 0; idx2 < n; idx2++) {
            A[n * i + idx2] = A[n * i + idx2] / scale;
        }
        for (s32 idx2 = 0; idx2 < n; idx2++) {
            B[n * i + idx2] = B[n * i + idx2] / scale;
        }
        for (s32 j = 0, idx2 = 0; idx2 < i; j += 1, idx2++) {
            double c;
            
            c = A[n * j + i];
            for (s32 idx3 = 0; idx3 < n; idx3++) {
                A[n * j + idx3] = A[n * j + idx3] - c * A[n * i + idx3];
            }
            for (s32 idx3 = 0; idx3 < n; idx3++) {
                B[n * j + idx3] = B[n * j + idx3] - c * B[n * i + idx3];
            }
        }
    }
    return 0;
}

