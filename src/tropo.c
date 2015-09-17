/*
 * Copyright (C) 2012 Swift Navigation Inc.
 * Contact: Fergus Noble <fergus@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#include <math.h>

#include "tropo.h"

/* Simple Black model, inspired by GPSTk SimpleTropModel class. */

static double dry_zenith_delay(void)
{
  return 2.235486646978727;
}

static double dry_mapping_function(double elevation)
{
  double d = cos(elevation);

  d /= 1.001012704615527;
  return 1.0 / sqrt(1.0 - d * d);
}

static double wet_zenith_delay(void)
{
  return 0.122382715318184;
}

static double wet_mapping_function(double elevation)
{
  double d = cos(elevation);

  d /= 1.000282213715744;
  return 1.0 / sqrt(1.0 - d * d);
}

double tropo_correction(double elevation)
{
  if (elevation < 0) {
    return 0;
  }

  return dry_zenith_delay() * dry_mapping_function(elevation)
         + wet_zenith_delay() * wet_mapping_function(elevation);
}

