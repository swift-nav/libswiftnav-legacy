#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "check_utils.h"

#define epsilon 0.0001

u8 within_epsilon(double a, double b) {
  if (fabs(a - b) < epsilon) {
    return 1;
  }
  return 0;
}

void seed_rng(void) {
  srandom(time(NULL));
}

double frand(double fmin, double fmax) {
  double f = (double)random() / RAND_MAX;
  return fmin + f * (fmax - fmin);
}

u32 sizerand(u32 sizemax) {
  double f = (double)random() / RAND_MAX;
  return (u32) ceil(f * sizemax);
}
