# Copyright (C) 2012 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

cdef extern from "libswiftnav/coord_system.h":
  void wgsllh2ecef(double llh[3], double ecef[3])
  void wgsecef2llh(double ecef[3], double llh[3])
  void wgsecef2ned(double ecef[3], double ref_ecef[3], double ned[3])
  void wgsecef2ned_d(double ecef[3], double ref_ecef[3], double ned[3])
  void wgsned2ecef(double ned[3], double ref_ecef[3], double ecef[3])
  void wgsned2ecef_d(double ned[3], double ref_ecef[3], double ecef[3])
  void wgsecef2azel(double ecef[3], double ref_ecef[3],
                    double* azimuth, double* elevation)

