#ifndef PLOVER_GENERATED_qr
#define PLOVER_GENERATED_qr

#include "plover/plvstdlib.h"

#include "common.h"
void newline (void);
double norm (const s32 n, const double * v);
void print_matrix (const s32 n, const s32 m, const double * A);
s8 qr_solve (const s32 m, const s32 n, double * A, double * b,
             double * solution, double * residual);
void qr_update (const s32 m, const s32 n, double * b, double * A);


#endif /* PLOVER_GENERATED_qr */
