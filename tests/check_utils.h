#include "common.h"

u8 within_epsilon(double a, double b);
u8 arr_within_epsilon(u32 n, const double *a, const double *b);
void seed_rng(void);
double frand(double fmin, double fmax);
void arr_frand(u32 n, double fmin, double fmax, double *v);
u32 sizerand(u32 sizemax);
