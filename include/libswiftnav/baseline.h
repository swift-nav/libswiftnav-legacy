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

#include <libswiftnav/common.h>
#include <libswiftnav/constants.h>
#include <libswiftnav/observation.h>

/** \addtogroup baseline
 * \{ */

typedef struct {
  double amb;
  gnss_signal_t sid;
} ambiguity_t;

typedef struct {
  double ambs[MAX_CHANNELS-1];
  gnss_signal_t sids[MAX_CHANNELS];
  u8 n;
} ambiguities_t;

/** \} */

/** Default threshold value for lesq baseline raim check. */
#define DEFAULT_RAIM_THRESHOLD 5.5

void predict_carrier_obs(u8 num_dds, const double *N, const double *DE,
                         const double b[3], double *dd_obs);

void amb_from_baseline(u8 num_dds, const double *DE, const double *dd_obs,
                       const double b[3], s32 *N);

s8 lesq_solution_float(u8 num_dds, const double *dd_obs, const double *N,
                       const double *DE, double b[3], double *resid)
                         __attribute__((warn_unused_result));

s8 least_squares_solve_b_external_ambs(u8 num_dds, const double *ambs,
         const sdiff_t *sdiffs_with_ref_first, const double *dd_measurements,
         const double ref_ecef[3], double b[3],
         bool disable_raim, double raim_threshold);

void diff_ambs(gnss_signal_t ref_sid, u8 num_ambs, const ambiguity_t *amb_set,
               double *dd_ambs);
s8 baseline_(u8 num_sdiffs, const sdiff_t *sdiffs, const double ref_ecef[3],
             u8 num_ambs, const ambiguity_t *single_ambs,
             u8 *num_used, gnss_signal_t *used_sids, double b[3],
             bool disable_raim, double raim_threshold);
s8 baseline(u8 num_sdiffs, const sdiff_t *sdiffs, const double ref_ecef[3],
            const ambiguities_t *ambs, u8 *num_used, gnss_signal_t *used_sids,
            double b[3], bool disable_raim, double raim_threshold);

void ambiguities_init(ambiguities_t *ambs);
s8 lesq_solve_raim(u8 num_dds_u8, const double *dd_obs,
                   const double *N, const double *DE, double b[3],
                   bool disable_raim, double raim_threshold,
                   u8 *n_used, double *residuals, u8 *removed_obs);

#endif /* LIBSWIFTNAV_BASELINE_H */
