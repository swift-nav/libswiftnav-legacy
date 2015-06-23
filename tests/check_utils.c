#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "check_utils.h"

/*#define epsilon 0.0001*/
#define EPSILON 1e-5

u8 within_epsilon(double a, double b) {
  if (fabs(a - b) < EPSILON) {
    return 1;
  }
  return 0;
}

u8 arr_within_epsilon(u32 n, const double *a, const double *b) {
  for (u32 i=0; i < n; i++) {
    if (!within_epsilon(a[i], b[i])) {
      return false;
    }
  }
  return true;
}

void seed_rng(void) {
  srandom(time(NULL));
}

double frand(double fmin, double fmax) {
  double f = (double)random() / RAND_MAX;
  return fmin + f * (fmax - fmin);
}

void arr_frand(u32 n, double fmin, double fmax, double *v)
{
  for (u32 i=0; i < n; i++) {
    v[i] = frand(fmin, fmax);
  }
}

u32 sizerand(u32 sizemax) {
  double f = (double)random() / RAND_MAX;
  return (u32) ceil(f * sizemax);
}
