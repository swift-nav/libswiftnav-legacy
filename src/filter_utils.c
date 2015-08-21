/*
 * Copyright (C) 2014-2015 Swift Navigation Inc.
 * Contact: Ian Horn <ian@swift-nav.com>
 *          Fergus Noble <fergus@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */


#include <string.h>
#include <math.h>
#include <assert.h>

#include "logging.h"
#include "observation.h"
#include "constants.h"
#include "linear_algebra.h"
#include "filter_utils.h"

/** \defgroup filter_utils Utility functions for use in filters
 * \{ */

/** Measure the integer ambiguity just from the code and carrier measurements.
 * The expectation value of carrier + code / lambda is
 * integer ambiguity + bias. Currently, pseudorange bias can get up to 10s of
 * wavelengths for minutes at a time, so averaging carrier + code isn't
 * sufficient for determining the ambiguity. Regardless of bias, this is an
 * important measurement. It is especially useful as a simple initialization
 * of the float filter.
 *
 * \param carrier A carrier phase measurement in units of wavelengths.
 * \param code    A code measurement in the same units as GPS_L1_LAMBDA_NO_VAC.
 * \return An estimate of the integer ambiguity. Its expectation value is the
 *         integer ambiguity plus carrier and code bias.
 */
double simple_amb_measurement(double carrier, double code)
{
  return carrier + code / GPS_L1_LAMBDA_NO_VAC;
}

/*  Presumes that the first alm entry is the reference sat. */
s8 assign_de_mtx(u8 num_sats, const sdiff_t *sats_with_ref_first,
                 const double ref_ecef[3], double *DE)
{
  assert(sats_with_ref_first != NULL);
  assert(ref_ecef != NULL);
  assert(DE != NULL);

  if (num_sats <= 1) {
    log_debug("assign_de_mtx: not enough sats\n");
    return -1;
  }

  assert(num_sats > 1);

  /* Vector to reference satellite. */
  double e_0[3];
  vector_subtract(3, sats_with_ref_first[0].sat_pos, ref_ecef, e_0);
  vector_normalize(3, e_0);

  for (u8 i=1; i<num_sats; i++) {
    /* Vector to satellite i */
    double e_i[3];
    vector_subtract(3, sats_with_ref_first[i].sat_pos, ref_ecef, e_i);
    vector_normalize(3, e_i);
    /* DE row = e_i - e_0 */
    vector_subtract(3, e_i, e_0, &DE[3*(i-1)]);
  }

  return 0;
}

/** Computes log(1 + exp(x)).
 * log(1 + exp(x)) gets inaccurate/blows up for large positive or negative x.
 * This just computes them so that they match the more naive log1p(exp(x))
 * on its domain, and are reasonably accurate elsewhere. It doesn't make
 * much effort to be super accurate (though it's pretty good).
 * \param x   A double
 * \returns   log(1 + exp(x))
 * */
double log1pexp(double x)
{
  if (x > 34) return x;       //  15 for 32 bit floats,  42 for 128 bit floats
  if (x < -37) return exp(x); // -17 for 32 bit floats, -44 for 128 bit floats,
  return log1p(exp(x));
}

/** \} */

