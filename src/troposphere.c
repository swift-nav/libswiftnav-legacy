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

#include "troposphere.h"

/** \defgroup troposphere Tropospheric models
 * Implemenations of tropospheric delay correction models.
 * \{ */

/** Average barometric pressure lookup table [mbar] */
static const double p_avg_lut[5] =
{
  1013.25,
  1017.25,
  1015.75,
  1011.75,
  1013.00
};

/** Average temperature lookup table [K] */
static const double t_avg_lut[5] =
{
  299.65,
  294.15,
  283.15,
  272.15,
  263.65
};

/** Average water vapour pressure lookup table [mbar] */
static const double e_avg_lut[5] =
{
  26.31,
  21.79,
  11.66,
  6.78,
  4.11
};

/** Average temperature lapse rate lookup table [K/m] */
static const double b_avg_lut[5] =
{
  6.30e-3,
  6.05e-3,
  5.58e-3,
  5.39e-3,
  4.53e-3
};

/** Average water vapour pressure height factor lookup table */
static const double l_avg_lut[5] =
{
  2.77,
  3.15,
  2.57,
  1.81,
  1.55
};

/** Amplitude barometric pressure lookup table [mbar] */
static const double p_amp_lut[5] =
{
   0.00,
  -3.75,
  -2.25,
  -1.75,
  -0.50
};

/** Amplitude temperature lookup table [K] */
static const double t_amp_lut[5] =
{
  0.00,
  7.00,
  11.00,
  15.00,
  14.50
};

/** Amplitude water vapour pressure lookup table [mbar] */
static const double e_amp_lut[5] =
{
  0.00,
  8.85,
  7.24,
  5.36,
  3.39
};

/** Amplitude temperature lapse rate lookup table [K/m] */
static const double b_amp_lut[5] =
{
  0.00,
  0.25e-3,
  0.32e-3,
  0.81e-3,
  0.62e-3
};

/** Amplitude water vapour pressure height factor lookup table */
static const double l_amp_lut[5] =
{
  0.00,
  0.33,
  0.46,
  0.74,
  0.30
};

// TODO doc
static double lookup_param(double lat, double[5] lut) {
  if (lat <= 15.0) {
    return lut[0];
  } else if (lat >= 75.0) {
    return lut[4];
  } else {
    u8 i = (lat - 15.0) / 15.0;
    double lat_i = i * 15.0 + 15.0;
    return lut[i] + (lut[i + 1] - lut[i]) / 15.0 * (lat - lat_i);
  }
}

static double calc_param(double lat, double doy, double[5] avg_lut, double[5] amp_lut) {
  double avg = lookup_param(lat, avg_lut);
  double amp = lookup_param(lat, amp_lut);

// TODO check rad deg
// TODO check frac vs int doy
  return avg - amp * cos((doy - 28.0) * 2.0 * M_PI / 365.25);
}

double calc_troposphere(double elevation)
{
  double p = calc_param(lat, doy, p_avg_lut, p_amp_lut);
  double t = lookup_param(lat, t_avg_lut, t_amp_lut);
  double e = lookup_param(lat, e_avg_lut, e_amp_lut);
  double b = lookup_param(lat, b_avg_lut, b_amp_lut);
  double l = lookup_param(lat, l_avg_lut, l_amp_lut);

  const double r = 287.054;
  const double k_1 = 77.60;
  const double k_2 = 16.6;
  const double k_3 = 377600;
  g = ;
  double h = ;

  double dl = l + 1.0;

  double g_m = 9.784 * (1.0 - 2.66e-3 * cos(2.0 * lat) - 2.8e-7 * h);
  double t_m = t * (1.0 - b * r / g_m * dl);

  double zhd = 10e-6 * k_1 * r / g_m * p * powf(1.0 - b * h / t, g / r * b);
  double zwd = 10e-6 * (t_m * k_2 + k_3) * r / (g_m * dl - b * r) * e / t * powf(1.0 - b * h / t, dl * g / (r * b) - 1.0);

  m_h = ;
  m_w = ;

  double ztd = m_h * zhd + m_w * zwd;
}

/* Simple Black model, inspired by GPSTk SimpleTropModel class. */

static double dry_zenith_delay(void)
{
  return 2.235486646978727;
}

static double dry_mapping_function(double elevation)
{
  double d = cos(elevation);
  d /= 1.001012704615527;
  return (1.0 / sqrt(1.0 - d*d));
}

static double wet_zenith_delay(void)
{
  return 0.122382715318184;
}

static double wet_mapping_function(double elevation)
{
  double d = cos(elevation);
  d /= 1.000282213715744;
  return (1.0 / sqrt(1.0 - d*d));
}
/*
double calc_troposphere(double elevation)
{
  if (elevation < 0)
    return 0;

  return (dry_zenith_delay() * dry_mapping_function(elevation)
        + wet_zenith_delay() * wet_mapping_function(elevation));
}
*/
/** \} */
