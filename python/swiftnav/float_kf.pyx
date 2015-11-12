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
from libc.string cimport memcpy, memcmp, memset
from almanac cimport *
from almanac_c cimport *
from gpstime cimport *
from gpstime_c cimport *
from single_diff_c cimport *
from dgnss_management_c cimport *

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
  float_kf_c.predict_forward(&(kf.kf))

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] state_mean_ = \
    np.array(state_mean, dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] state_cov_U_ = \
    np.array(state_cov_UDU.U, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] state_cov_D_ = \
    np.array(state_cov_UDU.D, dtype=np.double)

  float_kf_c.predict_forward(&(kf.kf))

  memcpy(&state_mean_[0], kf.kf.state_mean, kf.state_dim * sizeof(double))
  memcpy(&state_cov_U_[0,0], kf.kf.state_cov_U, kf.state_dim * kf.state_dim * sizeof(double))
  memcpy(&state_cov_D_[0], kf.kf.state_cov_D, kf.state_dim * sizeof(double))
  
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

def update_for_obs(KalmanFilter kf, decor_obs):
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] intermediate_mean_ = \
    np.empty(kf.state_dim, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] intermediate_cov_U_ = \
    np.empty((kf.state_dim, kf.state_dim), dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] intermediate_cov_D_ = \
    np.empty(kf.state_dim, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] decor_obs_ =  \
    np.array(decor_obs, dtype=np.double)

  float_kf_c.incorporate_obs(&(kf.kf), <double *> &decor_obs_[0])

  memcpy(&intermediate_mean_[0], kf.kf.state_mean, kf.state_dim * sizeof(double))
  memcpy(&intermediate_cov_U_[0,0], kf.kf.state_cov_U, kf.state_dim * kf.state_dim * sizeof(double))
  memcpy(&intermediate_cov_D_[0], kf.kf.state_cov_D, kf.state_dim * sizeof(double))


  return intermediate_mean_, UDU_decomposition(intermediate_cov_U_, intermediate_cov_D_)

def decorrelate(KalmanFilter kf, obs):
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] obs_ =  \
    np.array(obs, dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] decor_obs_ =  \
    np.empty(obs.shape, dtype=np.double)

  float_kf_c.decorrelate(&(kf.kf),
                        <double *> &obs_[0],
                        <double *> &decor_obs_[0])

  return decor_obs_

def filter_update(KalmanFilter kf, obs):
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] state_mean_ = \
    np.empty(kf.state_dim, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] state_cov_U_ = \
    np.empty((kf.state_dim, kf.state_dim), dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] state_cov_D_ = \
    np.empty(kf.state_dim, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] obs_ =  \
    np.array(obs, dtype=np.double)

  float_kf_c.kalman_filter_update(&(kf.kf), <double *> &obs_[0])

  memcpy(&state_mean_[0], kf.kf.state_mean, kf.state_dim * sizeof(double))
  memcpy(&state_cov_U_[0,0], kf.kf.state_cov_U, kf.state_dim * kf.state_dim * sizeof(double))
  memcpy(&state_cov_D_[0], kf.kf.state_cov_D, kf.state_dim * sizeof(double))

  return state_mean_, UDU_decomposition(state_cov_U_, state_cov_D_)

class UDU_decomposition:
  def __init__(self, U, D):
    self.U = U
    self.D = D
  def reconstruct(self):
    return reconstruct_udu(self)

cdef class KalmanFilter:
  def __init__(self,
               np.ndarray[np.double_t, ndim=2, mode="c"] transition_mtx,
               np.ndarray[np.double_t, ndim=2, mode="c"] transition_cov,
               np.ndarray[np.double_t, ndim=2, mode="c"] decor_mtx, 
               np.ndarray[np.double_t, ndim=2, mode="c"] decor_obs_mtx, 
               np.ndarray[np.double_t, ndim=1, mode="c"] decor_obs_cov,
               np.ndarray[np.double_t, ndim=1, mode="c"] state_mean,
               np.ndarray[np.double_t, ndim=2, mode="c"] state_cov_U,
               np.ndarray[np.double_t, ndim=1, mode="c"] state_cov_D,
               ):
    memset(&self.kf, 0, sizeof(kf_t))
    self.state_dim = transition_mtx.shape[0]
    self.obs_dim = decor_mtx.shape[0] * 2
    # self.num_sats = len(prns_with_ref_first)
    # self.prns_with_ref_first = prns_with_ref_first
    self.transition_mtx = transition_mtx
    self.transition_cov = transition_cov
    self.decor_mtx = decor_mtx
    self.decor_obs_mtx = decor_obs_mtx
    self.decor_obs_cov = decor_obs_cov
    self.state_mean = state_mean
    self.state_cov_U = state_cov_U
    self.state_cov_D = state_cov_D

  def __init__(self):
    memset(&self.kf, 0, sizeof(kf_t))

  def __repr__(self):
    return "<KalmanFilter with state_dim=" + str(self.state_dim) + \
           ", obs_dim=" + str(self.obs_dim) + ">"

  def testeq(KalmanFilter self, KalmanFilter other not None):
    state_dim_eq = self.state_dim == other.state_dim
    obs_dim_eq = self.obs_dim == other.obs_dim
    transition_mtx_eq = bool(np.all(self.transition_mtx == other.transition_mtx))
    transition_cov_eq = bool(np.all(self.transition_cov == other.transition_cov))
    decor_mtx_eq = bool(np.all(self.decor_mtx == other.decor_mtx))
    decor_obs_mtx_eq = bool(np.all(self.decor_obs_mtx == other.decor_obs_mtx))
    decor_obs_cov_eq = bool(np.all(self.decor_obs_cov == other.decor_obs_cov))
    state_mean_eq = bool(np.all(self.state_mean == other.state_mean))
    state_cov_U_eq = bool(np.all(self.state_cov_U == other.state_cov_U))
    state_cov_D_eq = bool(np.all(self.state_cov_D == other.state_cov_D))
    eq_dict = {'state_dim'      : state_dim_eq,
               'obs_dim'        : obs_dim_eq,
               'transition_mtx' : transition_mtx_eq,
               'transition_cov' : transition_cov_eq,
               'decor_mtx'      : decor_mtx_eq,
               'decor_obs_mtx'  : decor_obs_mtx_eq,
               'decor_obs_cov'  : decor_obs_cov_eq,
               'state_mean'     : state_mean_eq,
               'state_cov_U'    : state_cov_U_eq,
               'state_cov_D'    : state_cov_D_eq}
    return all(eq_dict.values()), eq_dict



  def cmp(KalmanFilter self, KalmanFilter other not None):
    return memcmp(&self.kf, &other.kf, sizeof(kf_t))

  def __richcmp__(KalmanFilter self, KalmanFilter other not None, int cmp_type):
    if not cmp_type == 2:
      raise NotImplementedError()
    return  0 == memcmp(&self.kf, &other.kf, sizeof(kf_t))

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

  # property num_sats:
  #   def __get__(self):
  #     return self.kf.num_sats
  #   def __set__(self, num_sats):
  #     self.kf.num_sats = num_sats

  # property prns_with_ref_first:
  #   def __get__(self):
  #     cdef np.ndarray[np.uint8_t, ndim=1, mode="c"] prns =\
  #       np.empty(self.num_sats, dtype=np.uint8)
  #     memcpy(&prns[0], self.kf.prns_with_ref_first, self.num_sats * sizeof(u8))
  #     return prns
  #   def __set__(self, np.ndarray[np.uint8_t, ndim=1, mode="c"] prns_with_ref_first):
  #     self.num_sats = len(prns_with_ref_first)
  #     memcpy(self.kf.prns_with_ref_first, &prns_with_ref_first[0], self.num_sats * sizeof(u8))
      

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
        np.empty((self.obs_dim, self.state_dim), dtype=np.double)
      memcpy(&decor_obs_mtx[0,0], self.kf.decor_obs_mtx, self.obs_dim * self.state_dim * sizeof(double))
      return decor_obs_mtx
    def __set__(self, np.ndarray[np.double_t, ndim=2, mode="c"] decor_obs_mtx):
      memcpy(self.kf.decor_obs_mtx, &decor_obs_mtx[0,0], self.obs_dim * self.state_dim * sizeof(double))

  property decor_obs_cov:
    def __get__(self):
      cdef np.ndarray[np.double_t, ndim=1, mode="c"] decor_obs_cov = \
        np.empty(self.obs_dim, dtype=np.double)
      memcpy(&decor_obs_cov[0], self.kf.decor_obs_cov, self.obs_dim * sizeof(double))
      return decor_obs_cov
    def __set__(self, np.ndarray[np.double_t, ndim=1, mode="c"] decor_obs_cov):
      memcpy(self.kf.decor_obs_cov, &decor_obs_cov[0], self.obs_dim * sizeof(double))

  property state_mean:
    def __get__(self):
      cdef np.ndarray[np.double_t, ndim=1, mode="c"] state_mean = \
        np.empty(self.state_dim, dtype=np.double)
      memcpy(&state_mean[0], self.kf.state_mean, self.state_dim * sizeof(double))
      return state_mean
    def __set__(self, np.ndarray[np.double_t, ndim=1, mode="c"] state_mean):
      memcpy(self.kf.state_mean, &state_mean[0], self.state_dim * sizeof(double))

  property state_cov_U:
    def __get__(self):
      cdef np.ndarray[np.double_t, ndim=2, mode="c"] state_cov_U = \
        np.empty((self.state_dim, self.state_dim), dtype=np.double)
      memcpy(&state_cov_U[0,0], self.kf.state_cov_U, self.state_dim * self.state_dim * sizeof(double))
      return state_cov_U
    def __set__(self, np.ndarray[np.double_t, ndim=2, mode="c"] state_cov_U):
      memcpy(self.kf.state_cov_U, &state_cov_U[0,0], self.state_dim * self.state_dim * sizeof(double))

  property state_cov_D:
    def __get__(self):
      cdef np.ndarray[np.double_t, ndim=1, mode="c"] state_cov_D = \
        np.empty(self.state_dim, dtype=np.double)
      memcpy(&state_cov_D[0], self.kf.state_cov_D, self.state_dim * sizeof(double))
      return state_cov_D
    def __set__(self, np.ndarray[np.double_t, ndim=1, mode="c"] state_cov_D):
      memcpy(self.kf.state_cov_D, &state_cov_D[0], self.state_dim * sizeof(double))

  

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

def get_decor_obs_cov(num_diffs, phase_var, code_var):
  
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] decor_obs_cov = \
    np.empty(2 * num_diffs, dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] decor_mtx = \
        np.empty((num_diffs, num_diffs), dtype=np.double)

  float_kf_c.assign_decor_obs_cov(num_diffs, phase_var, code_var, &decor_mtx[0,0], &decor_obs_cov[0])

  return decor_mtx, decor_obs_cov

def least_squares_solve(KalmanFilter kf, measurements):
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] measurements_ = \
    np.array(measurements, dtype=np.double)

  cdef np.ndarray[np.double_t, ndim=1, mode="c"] lsq_state = \
    np.empty(max(kf.state_dim, kf.obs_dim), dtype=np.double)

  float_kf_c.least_squares_solve(&(kf.kf),
                                 &measurements_[0],
                                 &lsq_state[0])
  return lsq_state[:kf.state_dim]





