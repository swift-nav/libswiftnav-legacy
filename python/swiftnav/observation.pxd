# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from common cimport *
from ephemeris cimport ephemeris_t
from time cimport gps_time_t
from signal cimport gnss_signal_t
from track cimport navigation_measurement_t

cdef extern from "libswiftnav/observation.h":
  ctypedef struct sdiff_t:
    double pseudorange
    double carrier_phase
    double doppler
    double sat_pos[3]
    double sat_vel[3]
    double snr
    gnss_signal_t sid

  int cmp_sdiff(const void *a_, const void *b_)
  int cmp_amb(const void *a_, const void *b_)
  int cmp_amb_sdiff(const void *a_, const void *b_)
  int cmp_sdiff_sid(const void *a_, const void *b_)
  int cmp_amb_sid(const void *a_, const void *b_)

  u8 single_diff(u8 n_a, navigation_measurement_t *m_a,
                 u8 n_b, navigation_measurement_t *m_b,
                 sdiff_t *sds)
  int cmp_sid_sdiff(const void *a, const void *b)

  u8 make_propagated_sdiffs_wip(u8 n_local, navigation_measurement_t *m_local,
                                u8 n_remote, navigation_measurement_t *m_remote,
                                double remote_pos_ecef[3], sdiff_t *sds)
  u8 make_propagated_sdiffs(u8 n_local, navigation_measurement_t *m_local,
                            u8 n_remote, navigation_measurement_t *m_remote,
                            double *remote_dists, double remote_pos_ecef[3],
                            sdiff_t *sds)

  s8 make_dd_measurements_and_sdiffs(gnss_signal_t ref_sid, const gnss_signal_t *non_ref_sids,
                                     u8 num_dds, u8 num_sdiffs,
                                     const sdiff_t *sdiffs_in,
                                     double *dd_meas, sdiff_t *sdiffs_out)
  u8 check_lock_counters(u8 n_sds, const sdiff_t *sds, u16 *lock_counters[],
                         gnss_signal_t *sats_to_drop)

  s8 copy_sdiffs_put_ref_first(const gnss_signal_t ref_sid, const u8 num_sdiffs,
                             const sdiff_t *sdiffs,
                             sdiff_t *sdiffs_with_ref_first)
  u8 filter_sdiffs(u8 num_sdiffs, sdiff_t *sdiffs, u8 num_sats_to_drop,
                   gnss_signal_t *sats_to_drop);

  void debug_sdiff(sdiff_t sd)
  void debug_sdiffs(u8 n, sdiff_t *sds)

cdef class SingleDiff:
  cdef sdiff_t _thisptr

cdef mk_sdiff_array(py_sdiffs, u8 n_c_sdiffs, sdiff_t *c_sdiffs)
cdef read_sdiff_array(u8 n_c_sdiffs, sdiff_t *c_sdiffs)
