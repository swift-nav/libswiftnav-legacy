# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

import numpy as np
cimport numpy as np
cimport float_kf_c
from libc.string cimport memcpy
from almanac cimport *
from almanac_c cimport *
from gpstime cimport *
from gpstime_c cimport *

def udu(M):
  n = M.shape[0]
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] M_ = \
    np.array(M, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] U = \
    np.empty((n,n), dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] D = \
    np.empty(n, dtype=np.double)
  float_kf_c.udu(n, <double *> &M_[0,0], <double *> &U[0,0], <double *> &D[0])
  return UDU_decomposition(U, D)

def reconstruct_udu(ud):
  n = ud.D.shape[0]

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] U = \
    np.array(ud.U, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] D = \
    np.array(ud.D, dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] M = \
    np.empty((n,n), dtype=np.double)

  float_kf_c.reconstruct_udu(n, <double *> &U[0,0], <double *> &D[0], <double *> &M[0,0])
  return M

def predict_forward(KalmanFilter kf, state_mean, state_cov_UDU):
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] state_mean_ = \
    np.array(state_mean, dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] state_cov_U_ = \
    np.array(state_cov_UDU.U, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] state_cov_D_ = \
    np.array(state_cov_UDU.D, dtype=np.double)

  float_kf_c.predict_forward(&(kf.kf),
                             <double *> &state_mean_[0], 
                             <double *> &state_cov_U_[0,0], 
                             <double *> &state_cov_D_[0])
  return state_mean_, UDU_decomposition(state_cov_U_, state_cov_D_)

def update_scalar_measurement(h, R, U, D):
  state_dim = h.shape[0]

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] h_ = \
    np.array(h, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] U_ = \
    np.array(U, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] D_ = \
    np.array(D, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] k = \
    np.empty(state_dim, dtype=np.double)

  float_kf_c.update_scalar_measurement(state_dim,
                                       <double *> &h_[0],
                                       <double> R,
                                       <double *> &U_[0,0],
                                       <double *> &D_[0],
                                       <double *> &k[0])
  return UDU_decomposition(U_, D_), k

def update_for_obs(KalmanFilter kf,
                   intermediate_mean, intermediate_cov_UDU,
                   decor_obs):
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] intermediate_mean_ = \
    np.array(intermediate_mean, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] intermediate_cov_U_ = \
    np.array(intermediate_cov_UDU.U, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] intermediate_cov_D_ = \
    np.array(intermediate_cov_UDU.D, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] decor_obs_ =  \
    np.array(decor_obs, dtype=np.double)

  float_kf_c.update_for_obs(&(kf.kf),
                            <double *> &intermediate_mean_[0], 
                            <double *> &intermediate_cov_U_[0,0], 
                            <double *> &intermediate_cov_D_[0],
                            <double *> &decor_obs_[0])
  return intermediate_mean_, UDU_decomposition(intermediate_cov_U_, intermediate_cov_D_)

def decorrelate(KalmanFilter kf, obs):
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] obs_ =  \
    np.array(obs, dtype=np.double)

  float_kf_c.decorrelate(&(kf.kf),
                        <double *> &obs_[0])

  return obs_

def filter_update(KalmanFilter kf,
                     state_mean, state_cov_UDU,
                     obs):
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] state_mean_ = \
    np.array(state_mean, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] state_cov_U_ = \
    np.array(state_cov_UDU.U, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] state_cov_D_ = \
    np.array(state_cov_UDU.D, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] obs_ =  \
    np.array(obs, dtype=np.double)

  float_kf_c.filter_update(&(kf.kf),
                           <double *> &state_mean_[0], 
                           <double *> &state_cov_U_[0,0],
                           <double *> &state_cov_D_[0], 
                           <double *> &obs_[0])

  return state_mean_, UDU_decomposition(state_cov_U_, state_cov_D_)

class UDU_decomposition:
  def __init__(self, U, D):
    self.U = U
    self.D = D
  def reconstruct(self):
    return reconstruct_udu(self)

cdef class KalmanFilter:
  cdef float_kf_c.kf_t kf
  def __init__(self, 
               np.ndarray[np.double_t, ndim=2, mode="c"] transition_mtx,
               np.ndarray[np.double_t, ndim=2, mode="c"] transition_cov,
               np.ndarray[np.double_t, ndim=2, mode="c"] decor_mtx, 
               np.ndarray[np.double_t, ndim=2, mode="c"] decor_obs_mtx, 
               np.ndarray[np.double_t, ndim=1, mode="c"] decor_obs_cov):
    self.state_dim = transition_mtx.shape[0]
    self.obs_dim = decor_mtx.shape[0] * 2
    self.transition_mtx = transition_mtx
    self.transition_cov = transition_cov
    self.decor_mtx = decor_mtx
    self.decor_obs_mtx = decor_obs_mtx
    self.decor_obs_cov = decor_obs_cov

  def __repr__(self):
    return "<KalmanFilter with state_dim=" + str(self.state_dim) + \
           ", obs_dim=" + str(self.obs_dim) + ">"

  property state_dim:
    def __get__(self):
      return self.kf.state_dim
    def __set__(self, state_dim):
      self.kf.state_dim = state_dim

  property obs_dim:
    def __get__(self):
      return self.kf.obs_dim
    def __set__(self, obs_dim):
      self.kf.obs_dim = obs_dim

  property transition_mtx:
    def __get__(self):
      cdef np.ndarray[np.double_t, ndim=2, mode="c"] transition_mtx = \
        np.empty((self.state_dim, self.state_dim), dtype=np.double)
      memcpy(&transition_mtx[0,0], self.kf.transition_mtx, self.state_dim * self.state_dim * sizeof(double))
      return transition_mtx
    def __set__(self, np.ndarray[np.double_t, ndim=2, mode="c"] transition_mtx):
      memcpy(self.kf.transition_mtx, &transition_mtx[0,0], self.state_dim * self.state_dim * sizeof(double))

  property transition_cov:
    def __get__(self):
      cdef np.ndarray[np.double_t, ndim=2, mode="c"] transition_cov = \
        np.empty((self.state_dim, self.state_dim), dtype=np.double)
      memcpy(&transition_cov[0,0], self.kf.transition_cov, self.state_dim * self.state_dim * sizeof(double))
      return transition_cov
    def __set__(self, np.ndarray[np.double_t, ndim=2, mode="c"] transition_cov):
      memcpy(self.kf.transition_cov, &transition_cov[0,0], self.state_dim * self.state_dim * sizeof(double))

  property decor_mtx:
    def __get__(self):
      cdef np.ndarray[np.double_t, ndim=2, mode="c"] decor_mtx = \
        np.empty((self.obs_dim/2, self.obs_dim/2), dtype=np.double)
      memcpy(&decor_mtx[0,0], self.kf.decor_mtx, self.obs_dim * self.obs_dim * sizeof(double) / 4)
      return decor_mtx
    def __set__(self, np.ndarray[np.double_t, ndim=2, mode="c"] decor_mtx):
      memcpy(self.kf.decor_mtx, &decor_mtx[0,0], self.obs_dim * self.obs_dim * sizeof(double) / 4)

  property decor_obs_mtx:
    def __get__(self):
      cdef np.ndarray[np.double_t, ndim=2, mode="c"] decor_obs_mtx = \
        np.empty((self.obs_dim, self.obs_dim), dtype=np.double)
      memcpy(&decor_obs_mtx[0,0], self.kf.decor_obs_mtx, self.obs_dim * self.obs_dim * sizeof(double))
      return decor_obs_mtx
    def __set__(self, np.ndarray[np.double_t, ndim=2, mode="c"] decor_obs_mtx):
      memcpy(self.kf.decor_obs_mtx, &decor_obs_mtx[0,0], self.obs_dim * self.obs_dim * sizeof(double))

  property decor_obs_cov:
    def __get__(self):
      cdef np.ndarray[np.double_t, ndim=1, mode="c"] decor_obs_cov = \
        np.empty(self.obs_dim, dtype=np.double)
      memcpy(&decor_obs_cov[0], self.kf.decor_obs_cov, self.obs_dim * sizeof(double))
      return decor_obs_cov
    def __set__(self, np.ndarray[np.double_t, ndim=1, mode="c"] decor_obs_cov):
      memcpy(self.kf.decor_obs_cov, &decor_obs_cov[0], self.obs_dim * sizeof(double))

def get_transition_mtx(num_sats, dt):
  state_dim = num_sats + 5
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] transition_mtx = \
        np.empty((state_dim, state_dim), dtype=np.double)
  float_kf_c.assign_transition_mtx(state_dim, dt, &transition_mtx[0,0])
  return transition_mtx

def get_d_mtx(num_sats):
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] D = \
        np.empty((num_sats - 1, num_sats), dtype=np.double)
  float_kf_c.assign_d_mtx(num_sats, &D[0,0])
  return D

def get_e_mtx_from_alms(alms, GpsTime timestamp, ref_ecef):
  n = len(alms)
  cdef almanac_t al[32]
  cdef almanac_t a_
  for i, a in enumerate(alms):
    a_ = (<Almanac> a).almanac
    memcpy(&al[i], &a_, sizeof(almanac_t))
  
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(ref_ecef, dtype=np.double)

  cdef gps_time_t timestamp_ = timestamp.gps_time

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] e_mtx = \
        np.empty((len(alms), 3), dtype=np.double)

  float_kf_c.assign_e_mtx_from_alms(len(alms), &al[0], timestamp_, &ref_ecef_[0], &e_mtx[0,0])

  return e_mtx

def get_de_mtx_from_alms(alms, GpsTime timestamp, ref_ecef):
  n = len(alms)
  cdef almanac_t al[32]
  cdef almanac_t a_
  for i, a in enumerate(alms):
    a_ = (<Almanac> a).almanac
    memcpy(&al[i], &a_, sizeof(almanac_t))
  
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(ref_ecef, dtype=np.double)

  cdef gps_time_t timestamp_ = timestamp.gps_time

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] de_mtx = \
        np.empty((len(alms) - 1, 3), dtype=np.double)

  float_kf_c.assign_de_mtx_from_alms(len(alms), &al[0], timestamp_, &ref_ecef_[0], &de_mtx[0,0])

  return de_mtx

def get_obs_mtx_from_alms(alms, GpsTime timestamp, ref_ecef):
  n = len(alms)
  state_dim = n + 5
  obs_dim = 2 * (n-1)

  cdef almanac_t al[32]
  cdef almanac_t a_
  for i, a in enumerate(alms):
    a_ = (<Almanac> a).almanac
    memcpy(&al[i], &a_, sizeof(almanac_t))
  
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(ref_ecef, dtype=np.double)

  cdef gps_time_t timestamp_ = timestamp.gps_time

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] obs_mtx = \
        np.empty((obs_dim, state_dim), dtype=np.double)

  float_kf_c.assign_obs_mtx_from_alms(len(alms), &al[0], timestamp_, &ref_ecef_[0], &obs_mtx[0,0])

  return obs_mtx

def get_decor_obs_cov(num_diffs, phase_var, code_var):
  
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] decor_obs_cov = \
    np.empty(2 * num_diffs, dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] decor_mtx = \
        np.empty((num_diffs, num_diffs), dtype=np.double)

  float_kf_c.assign_decor_obs_cov(num_diffs, phase_var, code_var, &decor_mtx[0,0], &decor_obs_cov[0])

  return decor_mtx, decor_obs_cov


def get_decor_obs_mtx_from_alms(alms, GpsTime timestamp, ref_ecef, decor_mtx):
  n = len(alms)
  state_dim = n + 5
  obs_dim = 2 * (n-1)

  cdef almanac_t al[32]
  cdef almanac_t a_
  for i, a in enumerate(alms):
    a_ = (<Almanac> a).almanac
    memcpy(&al[i], &a_, sizeof(almanac_t))
  
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] ref_ecef_ = \
    np.array(ref_ecef, dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] decor_mtx_ = \
    np.array(decor_mtx, dtype=np.double)

  cdef gps_time_t timestamp_ = timestamp.gps_time

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] obs_mtx = \
        np.empty((obs_dim, state_dim), dtype=np.double)

  float_kf_c.assign_decor_obs_mtx_from_alms(len(alms), &al[0], timestamp_, &ref_ecef_[0],
                                            &decor_mtx_[0,0], &obs_mtx[0,0])

  return obs_mtx



