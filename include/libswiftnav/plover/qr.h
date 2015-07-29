#ifndef PLOVER_GENERATED_qr
#define PLOVER_GENERATED_qr

#include "plover/externs.h"

#include "common.h"
double norm (const s32 n, const double * v);
s8 qr_solve (const s32 m, const s32 n, double * A, double * b,
             double * solution, double * residual);
void qr_update (const s32 m, const s32 n, double * b, double * A);


#endif /* PLOVER_GENERATED_qr */
