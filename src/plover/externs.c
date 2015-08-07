#include "plover/externs.h"


#include <math.h>
#include <stdio.h>
#include "linear_algebra.h"
void newline (void)
{
    puts("");
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

