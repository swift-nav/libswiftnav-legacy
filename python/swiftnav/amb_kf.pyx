# Copyright (C) 2014 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

cimport numpy as np
from almanac cimport *
from gpstime cimport *
from gpstime cimport *
from libc.string cimport memcpy, memcmp, memset
from linear_algebra import *
from observation cimport *
import numpy as np

# def filter_update(KalmanFilter kf, obs):
#   cdef nkf_t kf_ = kf._thisptr
#   cdef np.ndarray[np.double_t, ndim=1, mode="c"] state_mean_ = np.empty(kf_.state_dim, dtype=np.double)
#   cdef np.ndarray[np.double_t, ndim=2, mode="c"] state_cov_U_ = np.empty((kf_.state_dim, kf_.state_dim), dtype=np.double)
#   cdef np.ndarray[np.double_t, ndim=1, mode="c"] state_cov_D_ = np.empty(kf_.state_dim, dtype=np.double)
#   cdef np.ndarray[np.double_t, ndim=1, mode="c"] obs_ =  np.array(obs, dtype=np.double)
#   nkf_update(&(kf_.kf), <double *> &obs_[0])
#   memcpy(&state_mean_[0], kf_.kf.state_mean, kf_.state_dim * sizeof(double))
#   memcpy(&state_cov_U_[0,0], kf_.kf.state_cov_U, kf_.state_dim * kf.state_dim * sizeof(double))
#   memcpy(&state_cov_D_[0], kf_.kf.state_cov_D, kf.state_dim * sizeof(double))
#   return state_mean_, UDU_decomposition(state_cov_U_, state_cov_D_)

cdef class KalmanFilter:

  def __init__(self,
               amb_drift_var,
               np.ndarray[np.double_t, ndim=2, mode="c"] decor_mtx,
               np.ndarray[np.double_t, ndim=2, mode="c"] decor_obs_mtx,
               np.ndarray[np.double_t, ndim=1, mode="c"] decor_obs_cov,
               np.ndarray[np.double_t, ndim=2, mode="c"] null_basis_Q,
               np.ndarray[np.double_t, ndim=1, mode="c"] state_mean,
               np.ndarray[np.double_t, ndim=2, mode="c"] state_cov_U,
               np.ndarray[np.double_t, ndim=1, mode="c"] state_cov_D,
               ):
    memset(&self._thisptr, 0, sizeof(nkf_t))
    self.state_dim = decor_mtx.shape[0]
    self.obs_dim = decor_mtx.shape[0] + max(0, decor_mtx.shape[0] - 3)
    self.amb_drift_var = amb_drift_var
    self.decor_mtx = decor_mtx
    self.decor_obs_mtx = decor_obs_mtx
    self.decor_obs_cov = decor_obs_cov
    self.null_basis_Q = null_basis_Q
    self.state_mean = state_mean
    self.state_cov_U = state_cov_U
    self.state_cov_D = state_cov_D

  def __init__(self):
    memset(&self._thisptr, 0, sizeof(nkf_t))

  def __repr__(self):
    return "<KalmanFilter with state_dim=" + str(self.state_dim) + \
           ", obs_dim=" + str(self.obs_dim) + ">"

  def testeq(KalmanFilter self, KalmanFilter other not None):
    state_dim_eq = self.state_dim == other.state_dim
    obs_dim_eq = self.obs_dim == other.obs_dim
    amb_drift_vareq = self.amb_drift_var == other.amb_drift_var
    decor_mtx_eq = bool(np.all(self.decor_mtx == other.decor_mtx))
    decor_obs_mtx_eq = bool(np.all(self.decor_obs_mtx == other.decor_obs_mtx))
    decor_obs_cov_eq = bool(np.all(self.decor_obs_cov == other.decor_obs_cov))
    null_basis_Q_eq = bool(np.all(self.null_basis_Q == other.null_basis_Q))
    state_mean_eq = bool(np.all(self.state_mean == other.state_mean))
    state_cov_U_eq = bool(np.all(self.state_cov_U == other.state_cov_U))
    state_cov_D_eq = bool(np.all(self.state_cov_D == other.state_cov_D))
    eq_dict = {'state_dim'      : state_dim_eq,
               'obs_dim'        : obs_dim_eq,
               'amb_drift_var'  : amb_drift_vareq,
               'decor_mtx'      : decor_mtx_eq,
               'decor_obs_mtx'  : decor_obs_mtx_eq,
               'decor_obs_cov'  : decor_obs_cov_eq,
               'null_basis_Q'   : null_basis_Q_eq,
               'state_mean'     : state_mean_eq,
               'state_cov_U'    : state_cov_U_eq,
               'state_cov_D'    : state_cov_D_eq}
    return all(eq_dict.values()), eq_dict

  # TODO/FIX (Buro):
  # def cmp(self, other):
  #   cdef nkf_t other_ = other._thisptr
  #   return memcmp(&self._thisptr, &(other_.thisptr), sizeof(nkf_t))

  # def __richcmp__(KalmanFilter self, KalmanFilter other not None, int cmp_type):
  #   if not cmp_type == 2:
  #     raise NotImplementedError()
  #   cdef nkf_t other_ = other._thisptr
  #   return  0 == memcmp(&self._thisptr, &other_._thisptr, sizeof(nkf_t))

  property state_dim:
    def __get__(self):
      # print -1
      return self._thisptr.state_dim
    def __set__(self, state_dim):
      self._thisptr.state_dim = state_dim

  property obs_dim:
    def __get__(self):
      # print -2
      return self._thisptr.obs_dim
    def __set__(self, obs_dim):
      self._thisptr.obs_dim = obs_dim

  property amb_drift_var:
    def __get__(self):
      return self._thisptr.amb_drift_var
    def __set__(self, amb_drift_var):
      self._thisptr.amb_drift_var = amb_drift_var

  property decor_mtx:
    def __get__(self):
      cdef np.ndarray[np.double_t, ndim=2, mode="c"] decor_mtx = \
        np.empty((self.obs_dim/2, self.obs_dim/2), dtype=np.double)
      memcpy(&decor_mtx[0,0], self._thisptr.decor_mtx, self.obs_dim * self.obs_dim * sizeof(double) / 4)
      return decor_mtx
    def __set__(self, np.ndarray[np.double_t, ndim=2, mode="c"] decor_mtx):
      memcpy(self._thisptr.decor_mtx, &decor_mtx[0,0], self.obs_dim * self.obs_dim * sizeof(double) / 4)
      print 15

  property decor_obs_mtx:
    def __get__(self):
      cdef np.ndarray[np.double_t, ndim=2, mode="c"] decor_obs_mtx = \
        np.empty((self.obs_dim, self.state_dim), dtype=np.double)
      memcpy(&decor_obs_mtx[0,0], self._thisptr.decor_obs_mtx, self.obs_dim * self.state_dim * sizeof(double))
      return decor_obs_mtx
    def __set__(self, np.ndarray[np.double_t, ndim=2, mode="c"] decor_obs_mtx):
      memcpy(self._thisptr.decor_obs_mtx, &decor_obs_mtx[0,0], self.obs_dim * self.state_dim * sizeof(double))

  property decor_obs_cov:
    def __get__(self):
      cdef np.ndarray[np.double_t, ndim=1, mode="c"] decor_obs_cov = \
        np.empty(self.obs_dim, dtype=np.double)
      memcpy(&decor_obs_cov[0], self._thisptr.decor_obs_cov, self.obs_dim * sizeof(double))
      return decor_obs_cov
    def __set__(self, np.ndarray[np.double_t, ndim=1, mode="c"] decor_obs_cov):
      memcpy(self._thisptr.decor_obs_cov, &decor_obs_cov[0], self.obs_dim * sizeof(double))

  property null_basis_Q:
    def __get__(self):
      cdef np.ndarray[np.double_t, ndim=2, mode="c"] null_basis_Q = \
        np.empty((self.obs_dim - self.state_dim, self.state_dim), dtype=np.double)
      memcpy(&null_basis_Q[0,0], self._thisptr.null_basis_Q, (self.obs_dim - self.state_dim) * self.state_dim * sizeof(double))
      return null_basis_Q
    def __set__(self, np.ndarray[np.double_t, ndim=2, mode="c"] null_basis_Q):
      memcpy(self._thisptr.null_basis_Q, &null_basis_Q[0,0], (self.obs_dim - self.state_dim) * self.state_dim * sizeof(double))

  property state_mean:
    def __get__(self):
      cdef np.ndarray[np.double_t, ndim=1, mode="c"] state_mean = \
        np.empty(self.state_dim, dtype=np.double)
      memcpy(&state_mean[0], self._thisptr.state_mean, self.state_dim * sizeof(double))
      return state_mean
    def __set__(self, np.ndarray[np.double_t, ndim=1, mode="c"] state_mean):
      memcpy(self._thisptr.state_mean, &state_mean[0], self.state_dim * sizeof(double))

  property state_cov_U:
    def __get__(self):
      cdef np.ndarray[np.double_t, ndim=2, mode="c"] state_cov_U = \
        np.empty((self.state_dim, self.state_dim), dtype=np.double)
      memcpy(&state_cov_U[0,0], self._thisptr.state_cov_U, self.state_dim * self.state_dim * sizeof(double))
      return state_cov_U
    def __set__(self, np.ndarray[np.double_t, ndim=2, mode="c"] state_cov_U):
      memcpy(self._thisptr.state_cov_U, &state_cov_U[0,0], self.state_dim * self.state_dim * sizeof(double))

  property state_cov_D:
    def __get__(self):
      cdef np.ndarray[np.double_t, ndim=1, mode="c"] state_cov_D = \
        np.empty(self.state_dim, dtype=np.double)
      memcpy(&state_cov_D[0], self._thisptr.state_cov_D, self.state_dim * sizeof(double))
      return state_cov_D
    def __set__(self, np.ndarray[np.double_t, ndim=1, mode="c"] state_cov_D):
      memcpy(self._thisptr.state_cov_D, &state_cov_D[0], self.state_dim * sizeof(double))
