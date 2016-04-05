# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from constants cimport *
from libcpp cimport bool
from observation cimport *
from signal cimport *

cdef extern from "libswiftnav/baseline.h":
  float DEFAULT_RAIM_THRESHOLD

  ctypedef struct ambiguity_t:
    double amb
    gnss_signal_t sid

  ctypedef struct ambiguities_t:
    double ambs[MAX_CHANNELS-1]
    gnss_signal_t sids[MAX_CHANNELS]
    u8 n

  enum:
    BASELINE_SUCCESS
    BASELINE_SUCCESS_RAIM_REPAIR
    BASELINE_SUCCESS_NO_RAIM
    BASELINE_FLOAT
    BASELINE_FLOAT_RAIM_REPAIR
    BASELINE_FLOAT_NO_RAIM
    BASELINE_FIXED
    BASELINE_FIXED_RAIM_REPAIR
    BASELINE_FIXED_NO_RAIM
    BASELINE_NOT_ENOUGH_SATS_ROVER
    BASELINE_NOT_ENOUGH_SATS_COMMON
    BASELINE_NOT_ENOUGH_SATS_RAIM
    BASELINE_RAIM_REPAIR_FAIL
    BASELINE_RAIM_REPAIR_MULTI_SOLNS
    BASELINE_NOT_ENOUGH_SATS_FLOAT
    BASELINE_DGELSY_FAIL

  void predict_carrier_obs(u8 num_dds, const double *N, const double *DE,
                         const double b[3], double *dd_obs)
  void amb_from_baseline(u8 num_dds, const double *DE, const double *dd_obs,
                         const double b[3], s32 *N)

  s8 lesq_solution_float(u8 num_dds, const double *dd_obs, const double *N,
                         const double *DE, double b[3], double *resid)
  s8 least_squares_solve_b_external_ambs(u8 num_dds,
                                         const double *ambs,
                                         const sdiff_t *sdiffs_with_ref_first,
                                         const double *dd_measurements,
                                         const double ref_ecef[3],
                                         double b[3],
                                         bool disable_raim,
                                         double raim_threshold)
  void diff_ambs(gnss_signal_t ref_sid, u8 num_ambs,
                 const ambiguity_t *amb_set, double *dd_ambs)

  s8 baseline_(u8 num_sdiffs, const sdiff_t *sdiffs, const double ref_ecef[3],
               u8 num_ambs, const ambiguity_t *single_ambs,
               u8 *num_used, double b[3],
               bool disable_raim, double raim_threshold)
  s8 baseline(u8 num_sdiffs, const sdiff_t *sdiffs, const double ref_ecef[3],
              const ambiguities_t *ambs, u8 *num_used, double b[3],
              bool disable_raim, double raim_threshold)

  void ambiguities_init(ambiguities_t *ambs)

  s8 lesq_solve_raim(u8 num_dds_u8, const double *dd_obs,
                     const double *N, const double *DE, double b[3],
                     bool disable_raim, double raim_threshold,
                     u8 *n_used, double *residuals, u8 *removed_obs)

cdef class Ambiguity:
  cdef ambiguity_t _thisptr

cdef class Ambiguities:
  cdef ambiguities_t _thisptr
