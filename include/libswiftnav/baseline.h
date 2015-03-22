/*
 * Copyright (C) 2014 Swift Navigation Inc.
 * Contact: Ian Horn <ian@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef LIBSWIFTNAV_BASELINE_H
#define LIBSWIFTNAV_BASELINE_H

#include "common.h"

void amb_from_baseline(u8 num_sats, double *DE, double *dd_meas,
                       double b[3], s32 *N);
void lesq_solution(u8 num_dds, double *dd_meas, s32 *N, double *DE,
                   double b[3], double *resid);

#endif /* LIBSWIFTNAV_BASELINE_H */

