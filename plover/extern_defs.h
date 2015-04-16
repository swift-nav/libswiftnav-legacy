#ifndef LIBSWIFTNAV_PLOVER_EXTERNS_H
#define LIBSWIFTNAV_PLOVER_EXTERNS_H

#include "constants.h"
#include "track.h"

void inverse(int dim, double arr[dim][dim], double output[dim][dim]);
void printDouble(double x);
void printInt(int x);
int randInt();
double randFloat();

#endif /* LIBSWIFTNAV_PLOVER_EXTERNS_H */
