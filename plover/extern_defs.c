#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "linear_algebra.h"
//#include "pvt_test.c"

void inverse(int dim, double *arr, double *output)
{
  matrix_inverse(dim, arr, output);
}

double randf()
{
  return (double)rand()/(double)(RAND_MAX);
}

void printNumber(double x)
{
  printf("%g\n", x);
}
