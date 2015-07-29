#ifndef PLOVER_GENERATED_qr
#define PLOVER_GENERATED_qr

#include "common.h"
s32 main (void);
void newline (void);
double norm (const s32 n, const double * v);
void print_matrix (const s32 n, const s32 m, const double * A);
double qr_solve (const s32 m, const s32 n, double * A, double * b,
                 double * solution);
void qr_update (const s32 m, const s32 n, double * b, double * A);


#endif /* PLOVER_GENERATED_qr */
