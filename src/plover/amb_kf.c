#include "stdio.h"
#include "constants.h"
#include "math.h"
#include "plover/amb_kf.h"
double simple_amb_measurement(double carrier, double code)
{
  return (carrier + (code / GPS_L1_LAMBDA_NO_VAC));
}
