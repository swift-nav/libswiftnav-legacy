# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from amb_kf cimport *
from ambiguity_test cimport *
from baseline cimport *
from common cimport *
from constants cimport *
from observation cimport *
from sats_management cimport *

# TODO (Buro): Note that in everything in DGNSS's state management is
# global. The implementation file will have to be refactored when the
# corresponding part in libswiftnav is refactored.

cdef extern from "libswiftnav/dgnss_management.h":

  ctypedef struct dgnss_settings_t:
    double phase_var_test
    double code_var_test
    double phase_var_kf
    double code_var_kf
    double pos_trans_var
    double vel_trans_var
    double int_trans_var
    double amb_drift_var
    double pos_init_var
    double vel_init_var
    double amb_init_var
    double new_int_var

  ctypedef struct ambiguity_state_t:
    ambiguities_t fixed_ambs
    ambiguities_t float_ambs

  void dgnss_set_settings(double phase_var_test, double code_var_test,
                          double phase_var_kf, double code_var_kf,
                          double amb_drift_var, double amb_init_var,
                          double new_int_var)
  void make_measurements(u8 num_diffs, const sdiff_t *sdiffs, double *raw_measurements)

  void dgnss_init(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3])
  void dgnss_update(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3],
                    bool disable_raim, double raim_threshold)
  void dgnss_rebase_ref(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3],
                        u8 old_prns[MAX_CHANNELS], sdiff_t *corrected_sdiffs)
  nkf_t * get_dgnss_nkf()
  sats_management_t * get_sats_management()
  ambiguity_test_t* get_ambiguity_test()

  s8 dgnss_iar_resolved()
  u32 dgnss_iar_num_hyps()
  u32 dgnss_iar_num_sats()
  s8 dgnss_iar_get_single_hyp(double *hyp)
  void dgnss_reset_iar()

  void dgnss_init_known_baseline(u8 num_sats, sdiff_t *sdiffs, double receiver_ecef[3], double b[3])
  void dgnss_update_ambiguity_state(ambiguity_state_t *s)
  s8 dgnss_baseline(u8 num_sdiffs, const sdiff_t *sdiffs,
                    const double ref_ecef[3], const ambiguity_state_t *s,
                    u8 *num_used, double b[3],
                    bool disable_raim, double raim_threshold)
  void measure_amb_kf_b(u8 num_sdiffs, sdiff_t *sdiffs,
                        const double receiver_ecef[3], double *b)
  void measure_b_with_external_ambs(u8 state_dim, const double *state_mean,
                                    u8 num_sdiffs, sdiff_t *sdiffs,
                                    const double receiver_ecef[3], double *b)
  void measure_iar_b_with_external_ambs(double *state_mean,
                                        u8 num_sdiffs, sdiff_t *sdiffs,
                                        double receiver_ecef[3],
                                        double *b)
  u8 get_amb_kf_de_and_phase(u8 num_sdiffs, sdiff_t *sdiffs,
                             double ref_ecef[3],
                             double *de, double *phase)
  u8 get_iar_de_and_phase(u8 num_sdiffs, sdiff_t *sdiffs,
                          double ref_ecef[3],
                          double *de, double *phase)
  u8 dgnss_iar_pool_contains(double *ambs)
  double dgnss_iar_pool_ll(u8 num_ambs, double *ambs)
  double dgnss_iar_pool_prob(u8 num_ambs, double *ambs)
  u8 get_amb_kf_mean(double *ambs)
  u8 get_amb_kf_cov(double *cov)
  u8 get_amb_kf_prns(u8 *prns)
  u8 get_amb_test_prns(u8 *prns)
  u8 dgnss_iar_MLE_ambs(s32 *ambs)

cdef class DGNSSSettings:
  cdef dgnss_settings_t _thisptr

cdef class AmbiguityState:
  cdef ambiguity_state_t _thisptr
