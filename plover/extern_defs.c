#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "linear_algebra.h"
#include <time.h>

void inverse(int dim, double arr[dim][dim], double output[dim][dim])
{
  matrix_inverse(dim, (const double *)arr, (double *)output);
}

int randInt()
{
  srand((int)time(NULL));
  return rand() % 10;
}

double randFloat()
{
  return (double)rand()/(double)(RAND_MAX);
}

void printDouble(double x)
{
  printf("%g\n", x);
}

void printInt(int x)
{
  printf("%i\n", x);
}

