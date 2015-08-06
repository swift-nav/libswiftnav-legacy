#include "plover/qr.h"


#include <math.h>
#include <stdio.h>
#include "linear_algebra.h"
static s8 backsolve (const s32 n, const double * U, double * b);
static void givens (const double a, const double b, double * result);
s8 backsolve (const s32 n, const double * U, double * b)
{
    for (s32 i = 0, idx = 0; idx < n; i += 1, idx++) {
        if (fabs(U[n * i + i]) < MATRIX_EPSILON) {
            return -1;
        }
    }
    b[n - 1] = b[n - 1] / U[n * (n - 1) + (n - 1)];
    for (s32 i = n - 1, idx = 0; idx < -1 + n; i += -1, idx++) {
        double sum = 0;
        
        for (s32 idx2 = 0; idx2 < -i + n; idx2++) {
            sum += U[n * (i - 1) + (i + idx2)] * b[i + idx2];
        }
        b[i - 1] = (b[i - 1] - sum) / U[n * (i - 1) + (i - 1)];
    }
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
void newline (void)
{
    puts("");
}
double norm (const s32 n, const double * v)
{
    double sum = 0;
    
    for (s32 idx = 0; idx < n; idx++) {
        sum += v[idx] * v[idx];
    }
    return sqrt(sum);
}
void print_matrix (const s32 n, const s32 m, const double * A)
{
    for (s32 i = 0, idx = 0; idx < n; i += 1, idx++) {
        for (s32 j = 0, idx2 = 0; idx2 < m; j += 1, idx2++) {
            printf("%f ", A[m * i + j]);
        }
        newline();
    }
}
s8 qr_solve (const s32 m, const s32 n, double * A, double * b,
             double * solution, double * residual)
{
    qr_update(m, n, b, A);
    
    s8 code;
    double arg [n * n];
    double arg2 [n];
    
    for (s32 idx = 0; idx < n; idx++) {
        for (s32 idx2 = 0; idx2 < n; idx2++) {
            arg[n * idx + idx2] = A[n * idx + idx2];
        }
    }
    for (s32 idx = 0; idx < n; idx++) {
        arg2[idx] = b[idx];
    }
    code = backsolve(n, arg, arg2);
    for (s32 idx = 0; idx < n; idx++) {
        b[idx] = arg2[idx];
    }
    for (s32 idx = 0; idx < n; idx++) {
        solution[idx] = b[idx];
    }
    
    double arg3 [m + -n];
    
    for (s32 idx = 0; idx < m + -n; idx++) {
        arg3[idx] = b[n + idx];
    }
    *residual = norm(m + -n, arg3);
    return code;
}
void qr_update (const s32 m, const s32 n, double * b, double * A)
{
    for (s32 j = 1, idx = 0; idx < n; j += 1, idx++) {
        for (s32 i = m, idx2 = 0; idx2 < -j + m; i += -1, idx2++) {
            double rot [2 * 2];
            
            givens(A[n * (i - 2) + (j - 1)], A[n * (i - 1) + (j - 1)], rot);
            for (s32 k = j, idx3 = 0; idx3 < 1 + -j + n; k += 1, idx3++) {
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
            
            double v [2];
            
            for (s32 idx3 = 0; idx3 < 2; idx3++) {
                v[idx3] = b[i - 2 + idx3];
            }
            for (s32 idx3 = 0; idx3 < 2; idx3++) {
                double sum = 0;
                
                for (s32 idx4 = 0; idx4 < 2; idx4++) {
                    sum += rot[2 * idx4 + idx3] * v[idx4];
                }
                b[i - 2 + idx3] = sum;
            }
        }
    }
}

