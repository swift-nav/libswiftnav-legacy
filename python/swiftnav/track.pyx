# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

# cython: embedsignature=True

from common cimport *
from ephemeris cimport ephemeris_t
from ephemeris cimport Ephemeris
from ephemeris import Ephemeris
from gpstime cimport *
from gpstime import GpsTime
from libc.stdlib cimport malloc, free

# Discriminators

def calc_loop_gains_(float bw, float zeta, float k, float loop_freq):
  """
  Wraps function :libswiftnav:`calc_loop_gains`.

  Parameters
  ----------
  bw : float
    The loop noise bandwidth
  zeta : float
    The damping ratio
  k : float
    The loop gain
  loop_freq : float
    The sampling frequency

  Returns
  -------
  out : (float, float)
    The tuple `(pgain, igain)`.

  """
  cdef float pgain, igain
  calc_loop_gains(bw, zeta, k, loop_freq, &pgain, &igain)
  return (pgain, igain)

def costas_discriminator_(complex p):
  """
  Wraps the function :libswiftnav:`costas_discriminator`.

  Parameters
  ----------
  p : complex, :math:`I_P + Q_P j`
    The prompt correlation. The real component contains the in-phase
    correlation and the imaginary component contains the quadrature
    correlation.

  Returns
  -------
  out : float
    The discriminator value.

  """
  return costas_discriminator_(p.real, p.imag)

def frequency_discriminator_(complex p, complex prev_p):
  """
  Wraps the function :libswiftnav:`frequency_discriminator`.

  Parameters
  ----------
  p : complex, :math:`I_P + Q_P j`
    The prompt correlation. The real component contains the in-phase
    correlation and the imaginary component contains the quadrature
    correlation.
  prev_p : complex, :math:`I_P + Q_P j`
    The prompt correlation. The real component contains the in-phase
    correlation and the imaginary component contains the quadrature
    correlation.

  Returns
  -------
  out : float
    The discriminator value.

  """
  return frequency_discriminator(p.real, p.imag, prev_p.real, prev_p.imag)

cdef class Correlation:

  def __init__(self, **kwargs):
    self._thisptr = kwargs

def dll_discriminator_(cs):
  """
  Wraps the function :libswiftnav:`dll_discriminator`.

  Parameters
  ----------
  cs : [complex], :math:`I_E + Q_E j`
    The early correlation. The real component contains the in-phase
    correlation and the imaginary component contains the quadrature
    correlation. The prompt correlation. The late correlation.

  Returns
  -------
  out : float
    The discriminator value.

  """
  # TODO (Buro): Make this array initialization less janky.
  cdef correlation_t cs_[3]
  cs_[0] = cs[0]._thisptr
  cs_[1] = cs[1]._thisptr
  cs_[2] = cs[2]._thisptr
  return dll_discriminator(cs_)

# Tracking loop: Aided

cdef class AidedLoopFilter:

  def __cinit__(self, y0, pgain, igain, aiding_igain):
    aided_lf_init(&self._thisptr, y0, pgain, igain, aiding_igain)

  def update(self, p_i_error, aiding_error):
    return aided_lf_update(&self._thisptr, p_i_error, aiding_error)

# Tracking loop: Simple

cdef class SimpleLoopFilter:
  """
  Wraps the `libswiftnav` simple first-order loop filter implementation.

  The loop filter state, :libswiftnav:`simple_lf_state_t` is maintained by the
  class instance.

  Parameters
  ----------
  y0 : float
    The initial output variable value.
  pgain : float
    The proportional gain.
  igain : float
    The integral gain.

  """

  def __cinit__(self, y0, pgain, igain):
    simple_lf_init(&self._thisptr, y0, pgain, igain)

  def update(self, error):
    """
    Wraps the function :libswiftnav:`simple_lf_update`.

    Parameters
    ----------
    error : float
      The error input.

    Returns
    -------
    out : float
      The updated value of the output variable.

    """
    return simple_lf_update(&self._thisptr, error)

cdef class SimpleTrackingLoop:
  """
  Wraps the `libswiftnav` simple second-order tracking loop implementation.

  The tracking loop state, :libswiftnav:`simple_tl_state_t` is maintained by
  the class instance.

  For a full description of the loop filter parameters, see
  :libswiftnav:`calc_loop_gains`.

  """

  def __cinit__(self, loop_freq,
                code_freq, code_bw, code_zeta, code_k,
                carr_freq, carr_bw, carr_zeta, carr_k):
    simple_tl_init(&self._thisptr, loop_freq,
                   code_freq, code_bw, code_zeta, code_k,
                   carr_freq, carr_bw, carr_zeta, carr_k)


  def update(self, cs):
    """
    Wraps the function :libswiftnav:`simple_tl_update`.

    Parameters
    ----------
    cs : [complex], :math:`I_E + Q_E j`
      The early correlation. The real component contains the in-phase
      correlation and the imaginary component contains the quadrature
      correlation. The prompt correlation. The late correlation.

    Returns
    -------
    out : (float, float)
      The tuple (code_freq, carrier_freq).

    """
    cdef correlation_t cs_[3]
    cs_[0] = cs[0]._thisptr
    cs_[1] = cs[1]._thisptr
    cs_[2] = cs[2]._thisptr
    simple_tl_update(&self._thisptr, cs_)
    return (self.code_freq, self.carr_freq)

cdef class AidedTrackingLoop:
  """
  Wraps the `libswiftnav` simple second-order tracking loop implementation
  and the aided second-order tracking loop implementation.

  The tracking loop state, :libswiftnav:`aided_tl_state_t` is maintained by
  the class instance.

  TODO, add carrier aiding to the code loop.

  For a full description of the loop filter parameters, see
  :libswiftnav:`calc_loop_gains`.

  """

  def __cinit__(self, **kwargs):
    self._thisptr = kwargs
    aided_tl_init(&self._thisptr,
                  self._thisptr.loop_freq,
                  self._thisptr.code_freq,
                  self._thisptr.code_bw,
                  self._thisptr.code_zeta,
                  self._thisptr.code_k,
                  self._thisptr.carr_to_code,
                  self._thisptr.carr_freq,
                  self._thisptr.carr_bw,
                  self._thisptr.carr_zeta,
                  self._thisptr.carr_k,
                  self._thisptr.carr_freq_b1)


  def retune(self, code_params, carr_params, loop_freq, carr_freq_igain, carr_to_code):
    """
    Retune the tracking loop.

    Parameters
    ----------
    code_params : (float, float, float)
      Code tracking loop parameter tuple, `(bw, zeta, k)`.
    carr_params : (float, float, float)
      Carrier tracking loop parameter tuple, `(bw, zeta, k)`.
    loop_freq : float
      The frequency with which loop updates are performed.
    carr_freq_igain : float
      FLL aiding gain

    """
    self.loop_freq = loop_freq
    self.code_bw, self.code_zeta, self.code_k = code_params
    self.carr_bw, self.carr_zeta, self.carr_k = carr_params
    self.carr_freq_igain = carr_freq_igain
    aided_tl_retune(&self._thisptr, self.loop_freq,
                            self.code_bw, self.code_zeta, self.code_k,
                            self.carr_to_code,
                            self.carr_bw, self.carr_zeta, self.carr_k,
                            self.carr_freq_igain)

  def update(self, cs):
    """
    Wraps the function :libswiftnav:`aided_tl_update`.

    Parameters
    ----------
    cs : [complex], :math:`I_E + Q_E j`
      The early correlation. The real component contains the in-phase
      correlation and the imaginary component contains the quadrature
      correlation. The prompt correlation. The late correlation.

    Returns
    -------
    out : (float, float)
      The tuple (code_freq, carrier_freq).

    """
    cdef correlation_t cs_[3]
    cs_[0] = cs[0]._thisptr
    cs_[1] = cs[1]._thisptr
    cs_[2] = cs[2]._thisptr
    aided_tl_update(&self._thisptr, cs_)
    return (self.code_freq, self.carr_freq)

cdef class CompTrackingLoop:
  """
  Wraps the `libswiftnav` code/carrier phase complimentary filter tracking loop
  implementation.

  The tracking loop state, :libswiftnav:`comp_tl_state_t` is maintained by
  the class instance.

  For a full description of the loop filter parameters, see
  :libswiftnav:`calc_loop_gains`.

  """

  def __cinit__(self, loop_freq,
                code_freq, code_bw, code_zeta, code_k,
                carr_freq, carr_bw, carr_zeta, carr_k,
                tau, cpc, sched):
    comp_tl_init(&self._thisptr, loop_freq,
                    code_freq, code_bw,
                    code_zeta, code_k,
                    carr_freq, carr_bw,
                    carr_zeta, carr_k,
                    tau, cpc, sched)

  def update(self, cs):
    """
    Wraps the function :libswiftnav:`comp_tl_update`.

    Parameters
    ----------
    cs : [complex], :math:`I_E + Q_E j`
      The early correlation. The real component contains the in-phase
      correlation and the imaginary component contains the quadrature
      correlation. The prompt correlation. The late correlation.

    Returns
    -------
    out : (float, float)
      The tuple (code_freq, carrier_freq).

    """
    cdef correlation_t cs_[3]
    cs_[0] = cs[0]._thisptr
    cs_[1] = cs[1]._thisptr
    cs_[2] = cs[2]._thisptr
    comp_tl_update(&self._thisptr, cs_)
    return (self._thisptr.code_freq, self._thisptr.carr_freq)

cdef class LockDetector:
  """
  Wraps the `libswiftnav` PLL lock detector implementation.

  The detector state, :libswiftnav:`lock_detect_t` is maintained by
  the class instance.

  """

  def __cinit__(self, **kwargs):
    self._thisptr = kwargs
    lock_detect_init(&self._thisptr,
                     self._thisptr.k1,
                     self._thisptr.k2,
                     self._thisptr.thislp,
                     self._thisptr.lo)

  def reinit(self, k1, k2, lp, lo):
    lock_detect_reinit(&self._thisptr, k1, k2, lp, lo)

  def update(self, I, Q, DT):
    lock_detect_update(&self._thisptr, I, Q, DT)
    return (self._thisptr.outo, self._thisptr.outp)


cdef class AliasDetector:

  def __cinit__(self, acc_len, time_diff):
    alias_detect_init(&self._thisptr, acc_len, time_diff)

  def first(self, I, Q):
    alias_detect_first(&self._thisptr, I, Q)

  def second(self, I, Q):
    return alias_detect_second(&self._thisptr, I, Q)

cdef class CN0Estimator:
  """
  Wraps the `libswiftnav` :math:`C / N_0` estimator implementation.

  The estimator state, :libswiftnav:`cn0_est_state_t` is maintained by
  the class instance.

  """

  def __cinit__(self, **kwargs):
    self._thisptr = kwargs
    cn0_est_init(&self._thisptr,
                 self._thisptr.bw,
                 self._thisptr.cn0_0,
                 self._thisptr.cutoff_freq,
                 self._thisptr.loop_freq)

  def update(self, I, Q):
    """
    Wraps the function :libswiftnav:`cn0_est`.

    Parameters
    ----------
    I : float
      The prompt in-phase correlation from the tracking correlator.
    Q : float
      The prompt quadrature correlation from the tracking correlator.

    Returns
    -------
    out : float
      The Carrier-to-Noise Density, :math:`C / N_0`, in dBHz.

    """
    return cn0_est(&self._thisptr, I, Q)


cdef class ChannelMeasurement:

  def __init__(self, **kwargs):
    self._thisptr = kwargs

cdef class NavigationMeasurement:

  def __init__(self, raw_pseudorange, pseudorange, carrier_phase, raw_doppler,
                doppler, sat_pos, sat_vel, snr, lock_time, tot, prn, lock_counter):
    self._thisptr.raw_pseudorange = raw_pseudorange
    self._thisptr.pseudorange = pseudorange
    self._thisptr.carrier_phase = carrier_phase
    self._thisptr.raw_doppler = raw_doppler
    self._thisptr.doppler = doppler
    for i in range(3):
      self._thisptr.sat_pos[i] = sat_pos[i]
      self._thisptr.sat_vel[i] = sat_vel[i]
    self._thisptr.snr = snr
    self._thisptr.lock_time = lock_time
    self._thisptr.tot.tow = tot.tow
    self._thisptr.tot.wn = tot.wn
    self._thisptr.prn = prn
    self._thisptr.lock_counter = lock_counter

  def __rich_cmp__(self, nav_msg, op):
    cdef navigation_measurement_t nav_msg_ = nav_msg._thisptr
    return nav_meas_cmp(<void *>(&self._thisptr), <void *>(&nav_msg_))

# TODO (Buro): Change NavigationMeasurement so that it's
# zero-initialized without having to do this bullshit.
def empty_nav_meas():
  return NavigationMeasurement(0, 0, 0, 0, 0, (0, 0, 0), (0, 0, 0), 0, 0, GpsTime(0, 0), 0, 0)

# TODO (Buro): Remove mallocs, etc. here. Do all this in-place
def _calc_navigation_measurement(chan_meas, nav_meas, nav_time, ephemerides):
  n_channels = len(chan_meas)
  nav_meas = [empty_nav_meas() for n in range(n_channels)]
  cdef channel_measurement_t** chan_meas_ = <channel_measurement_t**>malloc(n_channels*sizeof(channel_measurement_t *))
  cdef navigation_measurement_t** nav_meas_ = <navigation_measurement_t**>malloc(n_channels*sizeof(navigation_measurement_t *))
  cdef ephemeris_t** ephs = <ephemeris_t**>malloc(n_channels*sizeof(ephemeris_t *))
  for n in range(n_channels):
    chan_meas_[n] = &((<ChannelMeasurement ?>chan_meas[n])._thisptr)
    nav_meas_[n] = &((<NavigationMeasurement ?>nav_meas[n])._thisptr)
    ephs[n] = &((<Ephemeris ?>ephemerides[n])._thisptr)
  calc_navigation_measurement_(n_channels, chan_meas_, nav_meas_, nav_time, ephs)
  free(chan_meas_)
  free(nav_meas_)
  free(ephs)
  return nav_meas

# TODO (Buro): Remove mallocs, etc. here. Also, wow, this is awful
def _tdcp_doppler(m_new, m_old):
  n_new = len(m_new)
  n_old = len(m_old)
  n_corrected = min(n_new, n_old)
  m_corrected = [empty_nav_meas() for n in range(n_corrected)]
  cdef navigation_measurement_t** m_new_ = <navigation_measurement_t **>malloc(n_new*sizeof(navigation_measurement_t *))
  cdef navigation_measurement_t** m_old_ = <navigation_measurement_t **>malloc(n_old*sizeof(navigation_measurement_t *))
  cdef navigation_measurement_t* m_corrected_ = <navigation_measurement_t *>malloc(n_corrected*sizeof(navigation_measurement_t *))
  for n in range(n_new):
    m_new_[n] = &((<NavigationMeasurement?>m_new[n])._thisptr)
  for n in range(n_old):
    m_old_[n] = &((<NavigationMeasurement?>m_old[n])._thisptr)
  for n in range(n_corrected):
    m_corrected_[n] = (<NavigationMeasurement?>m_corrected[n])._thisptr
  n_written = tdcp_doppler(n_new, m_new_[0], n_old, m_old_[0], m_corrected_)
  free(m_new_)
  free(m_old_)
  free(m_corrected_)
  return m_corrected[:n_written]
