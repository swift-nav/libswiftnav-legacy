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
  # cdef np.ndarray[np.double_t, ndim=1, mode="c"] R_ = \
  #   np.array([R], dtype=np.double)

  float_kf_c.update_scalar_measurement(state_dim,
                                       <double *> &h_[0],
                                       <double> R,
                                       <double *> &U_[0,0],
                                       <double *> &D_[0],
                                       <double *> &k[0])
  return UDU_decomposition(U_, D_), k

def update_for_obs(decor_obs_mtx, decor_obs_cov,
                   intermediate_mean, intermediate_cov_U, intermediate_cov_D,
                   decor_obs):
  state_dim = intermediate_cov_D.shape[0]
  obs_dim = decor_obs.shape[0]

  cdef np.ndarray[np.double_t, ndim=2, mode="c"] decor_obs_mtx_ = \
    np.array(decor_obs_mtx, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] decor_obs_cov_ = \
    np.array(decor_obs_cov, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] intermediate_mean_ = \
    np.array(intermediate_mean, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=2, mode="c"] intermediate_cov_U_ = \
    np.array(intermediate_cov_U, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] intermediate_cov_D_ = \
    np.array(intermediate_cov_D, dtype=np.double)
  cdef np.ndarray[np.double_t, ndim=1, mode="c"] decor_obs_ =  \
    np.array(decor_obs, dtype=np.double)

  float_kf_c.update_for_obs(state_dim, obs_dim,
                            <double *> &decor_obs_mtx_[0,0],
                            <double *> &decor_obs_cov_[0], 
                            <double *> &intermediate_mean_[0], 
                            <double *> &intermediate_cov_U_[0,0], 
                            <double *> &intermediate_cov_D_[0],
                            <double *> &decor_obs_[0])
  return intermediate_mean_, UDU_decomposition(intermediate_cov_U_, intermediate_cov_D_)

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
               np.ndarray[np.double_t, ndim=2, mode="c"] obs_cov_root_inv, 
               np.ndarray[np.double_t, ndim=2, mode="c"] decor_obs_mtx, 
               np.ndarray[np.double_t, ndim=1, mode="c"] decor_obs_cov):
    self.state_dim = transition_mtx.shape[0]
    self.obs_dim = obs_cov_root_inv.shape[0]
    self.transition_mtx = transition_mtx
    self.transition_cov = transition_cov
    self.obs_cov_root_inv = obs_cov_root_inv
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

  property obs_cov_root_inv:
    def __get__(self):
      cdef np.ndarray[np.double_t, ndim=2, mode="c"] obs_cov_root_inv = \
        np.empty((self.obs_dim, self.obs_dim), dtype=np.double)
      memcpy(&obs_cov_root_inv[0,0], self.kf.obs_cov_root_inv, self.obs_dim * self.obs_dim * sizeof(double))
      return obs_cov_root_inv
    def __set__(self, np.ndarray[np.double_t, ndim=2, mode="c"] obs_cov_root_inv):
      memcpy(self.kf.obs_cov_root_inv, &obs_cov_root_inv[0,0], self.obs_dim * self.obs_dim * sizeof(double))

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
