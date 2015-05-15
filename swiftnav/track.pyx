# Copyright (C) 2012 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

# cython: embedsignature=True

cimport track_c
from nav_msg cimport NavMsg
from nav_msg import NavMsg
from gpstime cimport *
from gpstime_c cimport *
from gpstime import GpsTime
from ephemeris_c cimport ephemeris_t
from track_c cimport channel_measurement_t, navigation_measurement_t
from common cimport *
from libc.stdlib cimport malloc, free

cdef class ChannelMeasurement:
  cdef channel_measurement_t meas

  def __cinit__(self, prn, code_phase, code_freq, carr_phase, carr_freq, TOW_ms, rx_time, snr):
    self.meas.prn = prn
    self.meas.code_phase_chips = code_phase
    self.meas.code_phase_rate = code_freq
    self.meas.carrier_phase = carr_phase
    self.meas.carrier_freq = carr_freq
    self.meas.time_of_week_ms = TOW_ms
    self.meas.receiver_time = rx_time
    self.meas.snr = snr


  property prn:
    def __get__(self):
      return self.meas.prn

  def __repr__(self):
    return '<ChannelMeasurement ' + str((self.meas.prn,
                self.meas.code_phase_chips,
                self.meas.code_phase_rate,
                self.meas.carrier_phase,
                self.meas.carrier_freq,
                self.meas.time_of_week_ms,
                self.meas.receiver_time,
                self.meas.snr)) + '>'

cdef class NavigationMeasurement:
  def __cinit__(self, raw_pseudorange, pseudorange, carrier_phase, raw_doppler,
                doppler, sat_pos, sat_vel, snr, lock_time, tot, prn, lock_counter):
    self.meas.raw_pseudorange = raw_pseudorange
    self.meas.pseudorange = pseudorange
    self.meas.carrier_phase = carrier_phase
    self.meas.raw_doppler = raw_doppler
    self.meas.doppler = doppler
    for i in range(3):
      self.meas.sat_pos[i] = sat_pos[i]
      self.meas.sat_vel[i] = sat_vel[i]
    self.meas.snr = snr
    self.meas.lock_time = lock_time
    self.meas.tot.tow = tot.tow
    self.meas.tot.wn = tot.wn
    self.meas.prn = prn
    self.meas.lock_counter = lock_counter

  property raw_pseudorange:
    def __get__(self):
      return self.meas.raw_pseudorange
    def __set__(self, raw_pseudorange):
      self.meas.raw_pseudorange = raw_pseudorange

  property pseudorange:
    def __get__(self):
      return self.meas.pseudorange
    def __set__(self, pseudorange):
      self.meas.pseudorange = pseudorange

  property carrier_phase:
    def __get__(self):
      return self.meas.carrier_phase
    def __set__(self, carrier_phase):
      self.meas.carrier_phase = carrier_phase

  property raw_doppler:
    def __get__(self):
      return self.meas.raw_doppler
    def __set__(self, raw_doppler):
      self.meas.raw_doppler = raw_doppler

  property doppler:
    def __get__(self):
      return self.meas.doppler
    def __set__(self, doppler):
      self.meas.doppler = doppler

  property sat_pos:
    def __get__(self):
      return (self.meas.sat_pos[0], self.meas.sat_pos[1], self.meas.sat_pos[2])
    def __set__(self, sat_pos):
      for i in range(3):
        self.meas.sat_pos[i] = sat_pos[i]

  property sat_vel:
    def __get__(self):
      return (self.meas.sat_vel[0], self.meas.sat_vel[1], self.meas.sat_vel[2])
    def __set__(self, sat_vel):
      for i in range(3):
        self.meas.sat_vel[i] = sat_vel[i]

  property snr:
    def __get__(self):
      return self.meas.snr
    def __set__(self, snr):
      self.meas.snr = snr

  property lock_time:
    def __get__(self):
      return self.meas.lock_time
    def __set__(self, lock_time):
      self.meas.lock_time = lock_time

  property tot:
    def __get__(self):
      return GpsTime(self.meas.tot.wn, self.meas.tot.tow)
    def __set__(self, tot):
      self.meas.tot.tow = tot.tow
      self.meas.tot.wn = tot.wn

  property prn:
    def __get__(self):
      return self.meas.prn
    def __set__(self, prn):
      self.meas.prn = prn

  property lock_counter:
    def __get__(self):
      return self.meas.lock_counter
    def __set__(self, lock_counter):
      self.meas.lock_counter = lock_counter

  def __repr__(self):
    return '<NavigationMeasurement ' + \
           str((self.meas.tot.tow,
                self.meas.pseudorange,
                self.meas.carrier_phase,
                self.meas.doppler)) + '>'

def calc_navigation_measurement(double t, chan_meas, es):
  n_channels = len(chan_meas)
  nav_meas = [NavigationMeasurement(0, 0, 0, 0, 0, (0,0,0), (0,0,0), 0, 0, 0, 0, 0) for n in range(n_channels)]

  cdef channel_measurement_t** chan_meas_ptrs = <channel_measurement_t**>malloc(n_channels*sizeof(channel_measurement_t*))
  cdef navigation_measurement_t** nav_meas_ptrs = <navigation_measurement_t**>malloc(n_channels*sizeof(navigation_measurement_t*))
  cdef ephemeris_t** es_ptrs = <ephemeris_t**>malloc(n_channels*sizeof(ephemeris_t*))

  for n in range(n_channels):
    chan_meas_ptrs[n] = &((<ChannelMeasurement?>chan_meas[n]).meas)
    nav_meas_ptrs[n] = &((<NavigationMeasurement?>nav_meas[n]).meas)
    es_ptrs[n] = &((<NavMsg?>es[n]).eph)

  track_c.calc_navigation_measurement_(n_channels, chan_meas_ptrs, nav_meas_ptrs, t, es_ptrs)

  free(chan_meas_ptrs)
  free(nav_meas_ptrs)
  free(es_ptrs)

  return nav_meas

def calc_loop_gains(float bw, float zeta, float k, float loop_freq):
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
  cdef float pgain
  cdef float igain
  track_c.calc_loop_gains(bw, zeta, k, loop_freq, &pgain, &igain)
  return (pgain, igain)

def costas_discriminator(complex p):
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
  return track_c.costas_discriminator(p.real, p.imag)

def dll_discriminator(complex e, complex p, complex l):
  """
  Wraps the function :libswiftnav:`dll_discriminator`.

  Parameters
  ----------
  e : complex, :math:`I_E + Q_E j`
    The early correlation. The real component contains the in-phase
    correlation and the imaginary component contains the quadrature
    correlation.
  p : complex, :math:`I_P + Q_P j`
    The prompt correlation.
  l : complex, :math:`I_L + Q_L j`
    The late correlation.

  Returns
  -------
  out : float
    The discriminator value.

  """
  cdef track_c.correlation_t cs[3]
  cs[0].I = e.real
  cs[0].Q = e.imag
  cs[1].I = p.real
  cs[1].Q = p.imag
  cs[2].I = l.real
  cs[2].Q = l.imag

  return track_c.dll_discriminator(cs)

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

  cdef track_c.simple_lf_state_t s

  def __cinit__(self, y0, pgain, igain):
    track_c.simple_lf_init(&self.s, y0, pgain, igain)

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
    return track_c.simple_lf_update(&self.s, error)

  property y:
    """The output variable."""
    def __get__(self):
      return self.s.y

cdef class AidedLoopFilter:
  """
  Wraps the `libswiftnav` aided first-order loop filter implementation.

  It's a PI filter using an extra I aiding error term.

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
  aiding_igain : float
    The aiding integral gain.

  """

  cdef track_c.aided_lf_state_t s

  def __cinit__(self, y0, pgain, igain, aiding_igain):
    track_c.aided_lf_init(&self.s, y0, pgain, igain, aiding_igain)

  def update(self, p_i_error, aiding_error):
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
    return track_c.aided_lf_update(&self.s, p_i_error, aiding_error)

  property y:
    """The output variable."""
    def __get__(self):
      return self.s.y


cdef class SimpleTrackingLoop:
  """
  Wraps the `libswiftnav` simple second-order tracking loop implementation.

  The tracking loop state, :libswiftnav:`simple_tl_state_t` is maintained by
  the class instance.

  For a full description of the loop filter parameters, see
  :libswiftnav:`calc_loop_gains`.

  Parameters
  ----------
  code_params : (float, float, float)
    Code tracking loop parameter tuple, `(bw, zeta, k)`.
  carr_params : (float, float, float)
    Carrier tracking loop parameter tuple, `(bw, zeta, k)`.
  loop_freq : float
    The frequency with which loop updates are performed.

  """

  cdef track_c.simple_tl_state_t s
  cdef float loop_freq
  cdef float code_bw, code_zeta, code_k
  cdef float carr_bw, carr_zeta, carr_k

  def __cinit__(self, code_params, carr_params, loop_freq):
    self.loop_freq = loop_freq
    self.code_bw, self.code_zeta, self.code_k = code_params
    self.carr_bw, self.carr_zeta, self.carr_k = carr_params

  def start(self, code_freq, carr_freq):
    """
    (Re-)initialise the tracking loop.

    Parameters
    ----------
    code_freq : float
      The code phase rate (i.e. frequency).
    carr_freq : float
      The carrier frequency.

    """
    track_c.simple_tl_init(&self.s, self.loop_freq,
                           code_freq, self.code_bw,
                           self.code_zeta, self.code_k,
                           carr_freq, self.carr_bw,
                           self.carr_zeta, self.carr_k)

  def update(self, complex e, complex p, complex l):
    """
    Wraps the function :libswiftnav:`simple_tl_update`.

    Parameters
    ----------
    e : complex, :math:`I_E + Q_E j`
      The early correlation. The real component contains the in-phase
      correlation and the imaginary component contains the quadrature
      correlation.
    p : complex, :math:`I_P + Q_P j`
      The prompt correlation.
    l : complex, :math:`I_L + Q_L j`
      The late correlation.

    Returns
    -------
    out : (float, float)
      The tuple (code_freq, carrier_freq).

    """
    cdef track_c.correlation_t cs[3]
    cs[0].I = e.real
    cs[0].Q = e.imag
    cs[1].I = p.real
    cs[1].Q = p.imag
    cs[2].I = l.real
    cs[2].Q = l.imag
    track_c.simple_tl_update(&self.s, cs)
    return (self.code_freq, self.carr_freq)

  property code_freq:
    """The code phase rate (i.e. frequency)."""
    def __get__(self):
      return self.s.code_freq

  property carr_freq:
    """The carrier frequency."""
    def __get__(self):
      return self.s.carr_freq


cdef class AidedTrackingLoop:
  """
  Wraps the `libswiftnav` simple second-order tracking loop implementation
  and the aided second-order tracking loop implementation.

  The tracking loop state, :libswiftnav:`aided_tl_state_t` is maintained by
  the class instance.

  TODO, add carrier aiding to the code loop.

  For a full description of the loop filter parameters, see
  :libswiftnav:`calc_loop_gains`.

  Parameters
  ----------
  code_params : (float, float, float)
    Code tracking loop parameter tuple, `(bw, zeta, k)`.
  loop_freq : float
    The frequency with which loop updates are performed.

  """

  cdef track_c.aided_tl_state_t s
  cdef float loop_freq
  cdef float code_bw, code_zeta, code_k
  cdef float carr_bw, carr_zeta, carr_k
  cdef float carr_freq_igain

  def __cinit__(self, code_params, carr_params, loop_freq, carr_freq_igain):
    self.loop_freq = loop_freq
    self.code_bw, self.code_zeta, self.code_k = code_params
    self.carr_bw, self.carr_zeta, self.carr_k = carr_params
    self.carr_freq_igain = carr_freq_igain

  def start(self, code_freq, carr_freq):
    """
    (Re-)initialise the tracking loop.

    Parameters
    ----------
    code_freq : float
      The code phase rate (i.e. frequency).
    carr_freq : float
      The carrier frequency.

    """
    track_c.aided_tl_init(&self.s, self.loop_freq,
                           code_freq, self.code_bw,
                           self.code_zeta, self.code_k,
                           carr_freq, self.carr_bw,
                           self.carr_zeta, self.carr_k,
                           self.carr_freq_igain)

  def update(self, complex e, complex p, complex l):
    """
    Wraps the function :libswiftnav:`aided_tl_update`.

    Parameters
    ----------
    e : complex, :math:`I_E + Q_E j`
      The early correlation. The real component contains the in-phase
      correlation and the imaginary component contains the quadrature
      correlation.
    p : complex, :math:`I_P + Q_P j`
      The prompt correlation.
    l : complex, :math:`I_L + Q_L j`
      The late correlation.

    Returns
    -------
    out : (float, float)
      The tuple (code_freq, carrier_freq).

    """
    cdef track_c.correlation_t cs[3]
    cs[0].I = e.real
    cs[0].Q = e.imag
    cs[1].I = p.real
    cs[1].Q = p.imag
    cs[2].I = l.real
    cs[2].Q = l.imag
    track_c.aided_tl_update(&self.s, cs)
    return (self.code_freq, self.carr_freq)

  property code_freq:
    """The code phase rate (i.e. frequency)."""
    def __get__(self):
      return self.s.code_freq

  property carr_freq:
    """The carrier frequency."""
    def __get__(self):
      return self.s.carr_freq


cdef class CompTrackingLoop:
  """
  Wraps the `libswiftnav` code/carrier phase complimentary filter tracking loop
  implementation.

  The tracking loop state, :libswiftnav:`comp_tl_state_t` is maintained by
  the class instance.

  For a full description of the loop filter parameters, see
  :libswiftnav:`calc_loop_gains`.

  Parameters
  ----------
  code_params : (float, float, float)
    Code tracking loop parameter tuple, `(bw, zeta, k)`.
  carr_params : (float, float, float)
    Carrier tracking loop parameter tuple, `(bw, zeta, k)`.
  loop_freq : float
    The frequency with which loop updates are performed.
  tau : float
    The complimentary filter cross-over frequency.
  cpc : float
    The number of carrier cycles per complete code, or equivalently the ratio
    of the carrier frequency to the nominal code frequency.
  sched : int, optional
    The gain scheduling count.

  """

  cdef track_c.comp_tl_state_t s
  cdef float loop_freq, cpc, tau
  cdef float code_bw, code_zeta, code_k
  cdef float carr_bw, carr_zeta, carr_k
  cdef u32 sched

  def __cinit__(self, code_params, carr_params, loop_freq, tau, cpc, sched=0):
    self.loop_freq = loop_freq
    self.tau = tau
    self.cpc = cpc
    self.sched = sched
    self.code_bw, self.code_zeta, self.code_k = code_params
    self.carr_bw, self.carr_zeta, self.carr_k = carr_params

  def start(self, code_freq, carr_freq):
    """
    (Re-)initialise the tracking loop.

    Parameters
    ----------
    code_freq : float
      The code phase rate (i.e. frequency) difference from nominal.
    carr_freq : float
      The carrier frequency difference from nominal, i.e. Doppler shift.

    """
    track_c.comp_tl_init(&self.s, self.loop_freq,
                         code_freq, self.code_bw,
                         self.code_zeta, self.code_k,
                         carr_freq, self.carr_bw,
                         self.carr_zeta, self.carr_k,
                         self.tau, self.cpc, self.sched)

  def update(self, complex e, complex p, complex l):
    """
    Wraps the function :libswiftnav:`comp_tl_update`.

    Parameters
    ----------
    e : complex, :math:`I_E + Q_E j`
      The early correlation. The real component contains the in-phase
      correlation and the imaginary component contains the quadrature
      correlation.
    p : complex, :math:`I_P + Q_P j`
      The prompt correlation.
    l : complex, :math:`I_L + Q_L j`
      The late correlation.

    Returns
    -------
    out : (float, float)
      The tuple (code_freq, carrier_freq).

    """
    cdef track_c.correlation_t cs[3]
    cs[0].I = e.real
    cs[0].Q = e.imag
    cs[1].I = p.real
    cs[1].Q = p.imag
    cs[2].I = l.real
    cs[2].Q = l.imag
    track_c.comp_tl_update(&self.s, cs)
    return (self.code_freq, self.carr_freq)

  property code_freq:
    """The code phase rate (i.e. frequency)."""
    def __get__(self):
      return self.s.code_freq

  property carr_freq:
    """The carrier frequency."""
    def __get__(self):
      return self.s.carr_freq


cdef class CN0Estimator:
  """
  Wraps the `libswiftnav` :math:`C / N_0` estimator implementation.

  The estimator state, :libswiftnav:`cn0_est_state_t` is maintained by
  the class instance.

  Parameters
  ----------
  bw : float
    The loop noise bandwidth in Hz.
  cn0_0 : float
    The initial value of :math:`C / N_0` in dBHz.
  cutoff_freq : float
    The low-pass filter cutoff frequency, :math:`f_c`, in Hz.
  loop_freq : float
    The loop update frequency, :math:`f`, in Hz.

  """

  cdef track_c.cn0_est_state_t s

  def __cinit__(self, bw, cn0_0, cutoff_freq, loop_freq):
    track_c.cn0_est_init(&self.s, bw, cn0_0, cutoff_freq, loop_freq)

  def update(self, I):
    """
    Wraps the function :libswiftnav:`cn0_est`.

    Parameters
    ----------
    I : float
      The prompt in-phase correlation from the tracking correlators.

    Returns
    -------
    out : float
      The Carrier-to-Noise Density, :math:`C / N_0`, in dBHz.

    """
    return track_c.cn0_est(&self.s, I)

