/*
 * Copyright (C) 2014 Swift Navigation Inc.
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

#include "constants.h"
#include "sbp_utils.h"

/** \addtogroup sbp
 * \{ */

/** \defgroup sbp_utils SBP Utils
 * Convert to and from SBP message types and other useful functions.
 * \{ */

void sbp_make_gps_time(sbp_gps_time_t *t_out, gps_time_t *t_in, u8 flags)
{
  t_out->wn = t_in->wn;
  t_out->tow = round(t_in->tow * 1e3);
  t_out->ns = round((t_in->tow - t_out->tow*1e-3) * 1e9);
  t_out->flags = flags;
}

void sbp_make_pos_llh(sbp_pos_llh_t *pos_llh, gnss_solution *soln, u8 flags)
{
  pos_llh->tow = round(soln->time.tow * 1e3);
  pos_llh->lat = soln->pos_llh[0] * R2D;
  pos_llh->lon = soln->pos_llh[1] * R2D;
  pos_llh->height = soln->pos_llh[2];
  /* TODO: fill in accuracy fields. */
  pos_llh->h_accuracy = 0;
  pos_llh->v_accuracy = 0;
  pos_llh->n_sats = soln->n_used;
  pos_llh->flags = flags;
}

void sbp_make_vel_ned(sbp_vel_ned_t *vel_ned, gnss_solution *soln, u8 flags)
{
  vel_ned->tow = round(soln->time.tow * 1e3);
  vel_ned->n = round(soln->vel_ned[0] * 1e3);
  vel_ned->e = round(soln->vel_ned[1] * 1e3);
  vel_ned->d = round(soln->vel_ned[2] * 1e3);
  /* TODO: fill in accuracy fields. */
  vel_ned->h_accuracy = 0;
  vel_ned->v_accuracy = 0;
  vel_ned->n_sats = soln->n_used;
  vel_ned->flags = flags;
}

void sbp_make_dops(sbp_dops_t *dops_out, dops_t *dops_in)
{
  dops_out->pdop = round(dops_in->pdop * 100);
  dops_out->gdop = round(dops_in->gdop * 100);
  dops_out->tdop = round(dops_in->tdop * 100);
  dops_out->hdop = round(dops_in->hdop * 100);
  dops_out->vdop = round(dops_in->vdop * 100);
}

/** \} */
/** \} */

