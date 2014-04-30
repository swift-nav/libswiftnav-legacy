# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from single_diff_c cimport *
from amb_kf_c cimport *
from sats_management_c cimport *

cdef extern from "libswiftnav/dgnss_management.h":
  void make_measurements(u8 num_diffs, sdiff_t *sdiffs, double *raw_measurements)
  # void dgnss_init(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3], double b_init[3])
  void dgnss_init(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3], double dt)
  # void dgnss_update(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3], double dt, u8 filter_choice, double b[3])
  void dgnss_update(u8 num_sats, sdiff_t *sdiffs, double reciever_ecef[3], double dt)
  s8 dgnss_iar_resolved()
  u32 dgnss_iar_num_hyps()
  u32 dgnss_iar_num_sats()
  s8 dgnss_iar_get_single_hyp(double *hyp)
  nkf_t * get_dgnss_kf()
  s32 * get_stupid_filter_ints()
  sats_management_t * get_sats_management()
  void measure_amb_kf_b(double reciever_ecef[3], 
                  		u8 num_sdiffs, sdiff_t *sdiffs,
                  		double *b)
  void measure_b_with_external_ambs(double reciever_ecef[3],
                                  u8 num_sdiffs, sdiff_t *sdiffs,
                                  double *ambs,
                                  double *b)
  u8 get_amb_kf_de_and_phase(u8 num_sdiffs, sdiff_t *sdiffs,
                             double ref_ecef[3],
                             double *de, double *phase)
  u8 get_iar_de_and_phase(u8 num_sdiffs, sdiff_t *sdiffs,
                          double ref_ecef[3],
                          double *de, double *phase)