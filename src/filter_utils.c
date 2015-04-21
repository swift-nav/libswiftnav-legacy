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
#include "single_diff.h"
#include "linear_algebra.h"
#include "filter_utils.h"

/*  Presumes that the first alm entry is the reference sat. */
void assign_de_mtx(u8 num_sats, const sdiff_t *sats_with_ref_first,
                   const double ref_ecef[3], double *DE)
{
  DEBUG_ENTRY();

  if (DEBUG) {
    printf("num_sats = %u\nsdiff prns&positions = {\n", num_sats);
    for (u8 i=0; i < num_sats; i++) {
      printf("i = %u, prn = %u, \tpos = [%f, \t%f, \t%f]\n",
             i,
             sats_with_ref_first[i].prn,
             sats_with_ref_first[i].sat_pos[0],
             sats_with_ref_first[i].sat_pos[1],
             sats_with_ref_first[i].sat_pos[2]);
    }
    printf("}\nref_ecef = {%f, \t%f, \t%f}\n", ref_ecef[0], ref_ecef[1], ref_ecef[2]);
  }

  assert(num_sats > 1);
  u8 de_length = num_sats - 1;

  if (num_sats <= 1) {
    log_debug("not enough sats\n");
    DEBUG_EXIT();
    return;
  }

  memset(DE, 0, de_length * 3 * sizeof(double));
  double e0[3];
  double x0 = sats_with_ref_first[0].sat_pos[0] - ref_ecef[0];
  double y0 = sats_with_ref_first[0].sat_pos[1] - ref_ecef[1];
  double z0 = sats_with_ref_first[0].sat_pos[2] - ref_ecef[2];
  double norm0 = sqrt(x0*x0 + y0*y0 + z0*z0);
  e0[0] = x0 / norm0;
  e0[1] = y0 / norm0;
  e0[2] = z0 / norm0;
  for (u8 i=1; i<num_sats; i++) {
    double x = sats_with_ref_first[i].sat_pos[0] - ref_ecef[0];
    double y = sats_with_ref_first[i].sat_pos[1] - ref_ecef[1];
    double z = sats_with_ref_first[i].sat_pos[2] - ref_ecef[2];
    double norm = sqrt(x*x + y*y + z*z);
    DE[3*(i-1)] = x / norm - e0[0];
    DE[3*(i-1) + 1] = y / norm - e0[1];
    DE[3*(i-1) + 2] = z / norm - e0[2];
  }
  if (DEBUG) {
    MAT_PRINTF(DE, (num_sats-1), 3);
  }
  DEBUG_EXIT();
}

