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
  def __cinit__(self, pr, prr, TOT, sat_pos, sat_vel):
    self.meas.pseudorange = pr
    self.meas.pseudorange_rate = prr
    self.meas.TOT = TOT
    for i in range(3):
      self.meas.sat_pos[i] = sat_pos[i]
      self.meas.sat_vel[i] = sat_vel[i]

  def __repr__(self):
    return '<NavigationMeasurement ' + str((self.meas.pseudorange,
                self.meas.pseudorange_rate,
                self.meas.TOT,
                (self.meas.sat_pos[0], self.meas.sat_pos[1], self.meas.sat_pos[2]),
                (self.meas.sat_vel[0], self.meas.sat_vel[1], self.meas.sat_vel[2]))) + '>'

def calc_navigation_measurement(double t, chan_meas, es):
  n_channels = len(chan_meas)
  nav_meas = [NavigationMeasurement(0, 0, 0, (0,0,0), (0,0,0)) for n in range(n_channels)]

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

def calc_loop_gains(double bw, double zeta, double k, double loop_freq):
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
  cdef double pgain
  cdef double igain
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

