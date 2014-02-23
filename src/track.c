/*
 * Copyright (C) 2012 Swift Navigation Inc.
 * Contact: Fergus Noble <fergus@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "constants.h"
#include "prns.h"
#include "track.h"
#include "ephemeris.h"
#include "tropo.h"
#include "coord_system.h"

/** \defgroup track Tracking
 * Functions used in tracking.
 * \{ */

/** \defgroup track_loop Tracking Loops
 * Functions used by the tracking loops.
 * \{ */

/** Calculate coefficients for a 2nd order digital PLL / DLL loop filter.
 *
 * A second order digital PLL consists of a first-order filter and a
 * Numerically Controlled Oscillator (NCO). A linearised model of a digital
 * second order PLL is shown below:
 *
 * \image html 2nd_order_dpll.png Linearised digital second order PLL model.
 *
 * Where \f$K_d\f$ is the discriminator gain and \f$F[z]\f$ and \f$N[z]\f$ are
 * the loop filter and NCO transfer functions. The NCO is essentially an
 * accumulator and hence has a transfer function:
 *
 * \f[
 *   N[z] = \frac{K_0 z^{-1}}{1 - z^{-1}}
 * \f]
 *
 * The first-order loop filter is shown below:
 *
 * \image html 1st_order_loop_filter.png Digital loop filter block diagram.
 *
 * and has transfer function:
 *
 * \f[
 *   F[z] = \frac{(k_p+k_i) - k_p z^{-1}}{1 - z^{-1}}
 * \f]
 *
 * It is useful to be able to calculate the loop filter coefficients, \f$k_p\f$
 * and \f$k_i\f$ by comparison to the parameters used in analog PLL design. By
 * comparison between the digital and analog PLL transfer functions we can show
 * the digital loop filter coefficients are related to the analog PLL natural
 * frequency, \f$\omega_n\f$ and damping ratio, \f$\zeta\f$ by:
 *
 * \f[
 *   k_p = \frac{1}{k} \frac{8 \zeta \omega_n T}
 *         {4 + 4 \zeta \omega_n T + \omega_n^2 T^2}
 * \f]
 * \f[
 *   k_i = \frac{1}{k} \frac{4 \omega_n^2 T^2}
 *         {4 + 4 \zeta \omega_n T + \omega_n^2 T^2}
 * \f]
 *
 * Where \f$T\f$ is the loop update period and the overall loop gain, \f$k\f$,
 * is the product of the NCO and discriminator gains, \f$ k = K_0 K_d \f$. The
 * natural frequency is related to the loop noise bandwidth, \f$B_L\f$ by the
 * following relationship:
 *
 * \f[
 *   \omega_n = \frac{8 \zeta B_L}{4 \zeta^2 + 1}
 * \f]
 *
 * These coefficients are applicable to both the Carrier phase Costas loop and
 * the Code phase DLL.
 *
 * References:
 *  -# Performance analysis of an all-digital BPSK direct-sequence
 *     spread-spectrum IF receiver architecture.
 *     B-Y. Chung, C. Chien, H. Samueli, and R. Jain.
 *     IEEE Journal on Selected Areas in Communications, 11:1096–1107, 1993.
 *
 * \todo This math is all wrong, these slides show the analysis we want:
 *       http://www.compdsp.com/presentations/Jacobsen/abineau_dpll_analysis.pdf
 *
 * \param bw        The loop noise bandwidth, \f$B_L\f$.
 * \param zeta      The damping ratio, \f$\zeta\f$.
 * \param k         The loop gain, \f$k\f$.
 * \param loop_freq The loop update frequency, \f$1/T\f$.
 * \param pgain     Where to store the calculated proportional gain,
 *                  \f$k_p\f$.
 * \param igain     Where to store the calculated integral gain, \f$k_i\f$.
 */
void calc_loop_gains(float bw, float zeta, float k, float loop_freq,
                     float *pgain, float *igain)
{
  /* Find the natural frequency. */
  float omega_n = 8.f*bw*zeta / (4.f*zeta*zeta + 1.f);

  /* Some intermmediate values. */
/*
  float T = 1. / loop_freq;
  float denominator = k*(4 + 4*zeta*omega_n*T + omega_n*omega_n*T*T);

  *pgain = 8*zeta*omega_n*T / denominator;
  *igain = 4*omega_n*omega_n*T*T / denominator;
*/
  *igain = omega_n * omega_n / (k * loop_freq);
  *pgain = 2.f * zeta * omega_n / k;
}

/** Phase discriminator for a Costas loop.
 *
 * \image html costas_loop.png Costas loop block diagram.
 *
 * Implements the \f$\tan^{-1}\f$ Costas loop discriminator.
 *
 * \f[
 *   \varepsilon_k = \tan^{-1} \left(\frac{I_k}{Q_k}\right)
 * \f]
 *
 * References:
 *  -# Understanding GPS: Principles and Applications.
 *     Elliott D. Kaplan. Artech House, 1996.
 *
 * \todo Fix potential divide by zero if Q is zero.
 *
 * \param I The prompt in-phase correlation, \f$I_k\f$.
 * \param Q The prompt quadrature correlation, \f$Q_k\f$.
 * \return The discriminator value, \f$\varepsilon_k\f$.
 */
float costas_discriminator(float I, float Q)
{
  return atanf(Q / I) * (float)(1/(2*M_PI));
}

/** Normalised non-coherent early-minus-late envelope discriminator.
 *
 * Implements the normalised non-coherent early-minus-late envelope DLL
 * discriminator.
 *
 * \f[
 *   \varepsilon_k = \frac{1}{2} \frac{E - L}{E + L}
 * \f]
 *
 * where:
 *
 * \f[
 *   E = \sqrt{I^2_E + Q^2_E}
 * \f]
 * \f[
 *   L = \sqrt{I^2_L + Q^2_L}
 * \f]
 *
 * References:
 *  -# Understanding GPS: Principles and Applications.
 *     Elliott D. Kaplan. Artech House, 1996.
 *
 * \param cs An array [E, P, L] of correlation_t structs for the Early, Prompt
 *           and Late correlations.
 * \return The discriminator value, \f$\varepsilon_k\f$.
 */
float dll_discriminator(correlation_t cs[3])
{
  float early_mag = sqrtf((float)cs[0].I*cs[0].I + (float)cs[0].Q*cs[0].Q);
  float late_mag = sqrtf((float)cs[2].I*cs[2].I + (float)cs[2].Q*cs[2].Q);

  return 0.5f * (early_mag - late_mag) / (early_mag + late_mag);
}

/** Initialise a simple first-order loop filter.
 * The gains can be calculated using calc_loop_gains().
 *
 * \param s The loop filter state struct to initialise.
 * \param y0 The initial value of the output variable, \f$y_0\f$.
 * \param pgain The proportional gain, \f$k_p\f$.
 * \param igain The integral gain, \f$k_i\f$.
 */
void simple_lf_init(simple_lf_state_t *s, float y0,
                    float pgain, float igain)
{
  s->y = y0;
  s->prev_error = 0.f;
  s->pgain = pgain;
  s->igain = igain;
}

/** Update step for the simple first-order loop filter.
 *
 * Implements the first-order loop filter as shown below:
 *
 * \image html 1st_order_loop_filter.png Digital loop filter block diagram.
 *
 * with transfer function:
 *
 * \f[
 *   F[z] = \frac{(k_p+k_i) - k_p z^{-1}}{1 - z^{-1}}
 * \f]
 *
 * \param s The loop filter state struct.
 * \param error The error output from the discriminator, \f$\varepsilon_k\f$.
 * \return The updated output variable, \f$y_k\f$.
 */
float simple_lf_update(simple_lf_state_t *s, float error)
{
  s->y += s->pgain * (error - s->prev_error) + \
          s->igain * error;
  s->prev_error = error;

  return s->y;
}

/** Initialise a simple tracking loop.
 *
 * For a full description of the loop filter parameters, see calc_loop_gains().
 *
 * \param s The tracking loop state struct to initialise.
 * \param code_freq The initial code phase rate (i.e. frequency).
 * \param code_bw The code tracking loop noise bandwidth.
 * \param code_zeta The code tracking loop damping ratio.
 * \param code_k The code tracking loop gain.
 * \param carr_freq The initial carrier frequency.
 * \param carr_bw The carrier tracking loop noise bandwidth.
 * \param carr_zeta The carrier tracking loop damping ratio.
 * \param carr_k The carrier tracking loop gain.
 */
void simple_tl_init(simple_tl_state_t *s, float loop_freq,
                    float code_freq, float code_bw,
                    float code_zeta, float code_k,
                    float carr_freq, float carr_bw,
                    float carr_zeta, float carr_k)
{
  float pgain, igain;

  calc_loop_gains(code_bw, code_zeta, code_k, loop_freq, &pgain, &igain);
  s->code_freq = code_freq;
  simple_lf_init(&(s->code_filt), code_freq, pgain, igain);

  calc_loop_gains(carr_bw, carr_zeta, carr_k, loop_freq, &pgain, &igain);
  s->carr_freq = carr_freq;
  simple_lf_init(&(s->carr_filt), carr_freq, pgain, igain);
}

/** Update step for the simple tracking loop.
 *
 * Implements a basic second-order tracking loop. The code tracking loop is a
 * second-order DLL using dll_discriminator() as its discriminator function.
 * The carrier phase tracking loop is a second-order Costas loop using
 * costas_discriminator().
 *
 * The tracking loop output variables, i.e. code and carrier frequencies can be
 * read out directly from the state struct.
 *
 * \param s The tracking loop state struct.
 * \param cs An array [E, P, L] of correlation_t structs for the Early, Prompt
 *           and Late correlations.
 */
void simple_tl_update(simple_tl_state_t *s, correlation_t cs[3])
{
  float code_error = dll_discriminator(cs);
  s->code_freq = simple_lf_update(&(s->code_filt), -code_error);
  float carr_error = costas_discriminator(cs[1].I, cs[1].Q);
  s->carr_freq = simple_lf_update(&(s->carr_filt), carr_error);
}

/** Initialise a code/carrier phase complimentary filter tracking loop.
 *
 * For a full description of the loop filter parameters, see calc_loop_gains().
 *
 * This filter implements gain scheduling. Before `sched` iterations of the
 * loop filter the behaviour will be identical to the simple loop filter. After
 * `sched` iterations, carrier phase information will be used in the code
 * tracking loop.
 *
 * Note, this filter requires that the code and carrier frequencies are
 * expressed as a difference from the nominal frequncy (e.g. 1.023MHz nominal
 * GPS C/A code phase rate, 4.092MHz IF for the carrier).
 *
 * \param s The tracking loop state struct to initialise.
 * \param code_freq The initial code phase rate (i.e. frequency) difference
 *                  from nominal.
 * \param code_bw The code tracking loop noise bandwidth.
 * \param code_zeta The code tracking loop damping ratio.
 * \param code_k The code tracking loop gain.
 * \param carr_freq The initial carrier frequency difference from nominal, i.e.
 *                  Doppler shift.
 * \param carr_bw The carrier tracking loop noise bandwidth.
 * \param carr_zeta The carrier tracking loop damping ratio.
 * \param carr_k The carrier tracking loop gain.
 * \param tau The complimentary filter cross-over frequency.
 * \param cpc The number of carrier cycles per complete code, or equivalently
 *            the ratio of the carrier frequency to the nominal code frequency.
 * \param sched The gain scheduling count.
 */
void comp_tl_init(comp_tl_state_t *s, float loop_freq,
                  float code_freq, float code_bw,
                  float code_zeta, float code_k,
                  float carr_freq, float carr_bw,
                  float carr_zeta, float carr_k,
                  float tau, float cpc,
                  u32 sched)
{
  float pgain, igain;

  calc_loop_gains(code_bw, code_zeta, code_k, loop_freq, &pgain, &igain);
  s->code_freq = code_freq;
  simple_lf_init(&(s->code_filt), code_freq, pgain, igain);

  calc_loop_gains(carr_bw, carr_zeta, carr_k, loop_freq, &pgain, &igain);
  s->carr_freq = carr_freq;
  simple_lf_init(&(s->carr_filt), carr_freq, pgain, igain);

  s->n = 0;
  s->sched = sched;
  s->carr_to_code = 1.f / cpc;

  s->A = 1.f - (1.f / (loop_freq * tau));
}

/** Update step for a code/carrier phase complimentary filter tracking loop.
 *
 * The tracking loop output variables, i.e. code and carrier frequencies can be
 * read out directly from the state struct.
 *
 * \todo Write proper documentation with math and diagram.
 *
 * \param s The tracking loop state struct.
 * \param cs An array [E, P, L] of correlation_t structs for the Early, Prompt
 *           and Late correlations.
 */
void comp_tl_update(comp_tl_state_t *s, correlation_t cs[3])
{
  float carr_error = costas_discriminator(cs[1].I, cs[1].Q);
  s->carr_freq = simple_lf_update(&(s->carr_filt), carr_error);

  float code_error = dll_discriminator(cs);
  s->code_filt.y = 0.f;
  float code_update = simple_lf_update(&(s->code_filt), -code_error);

  if (s->n > s->sched) {
    s->code_freq = s->A * s->code_freq + \
                   s->A * code_update + \
                   (1.f - s->A)*s->carr_to_code*s->carr_freq;
  } else {
    s->code_freq += code_update;
  }

  s->n++;
}


/** \} */

/** Initialise the \f$ C / N_0 \f$ estimator state.
 *
 * See cn0_est() for a full description.
 *
 * \param s The estimator state struct to initialise.
 * \param bw The loop noise bandwidth in Hz.
 * \param cn0_0 The initial value of \f$ C / N_0 \f$ in dBHz.
 * \param cutoff_freq The low-pass filter cutoff frequency, \f$f_c\f$, in Hz.
 * \param loop_freq The loop update frequency, \f$f\f$, in Hz.
 */
void cn0_est_init(cn0_est_state_t *s, float bw, float cn0_0,
                  float cutoff_freq, float loop_freq)
{
  s->log_bw = 10.f*log10f(bw);
  s->A = cutoff_freq / (loop_freq + cutoff_freq);
  s->I_prev_abs = -1.f;
  s->nsr = powf(10.f, 0.1f*(s->log_bw - cn0_0));
}

/** Estimate the Carrier-to-Noise Density, \f$ C / N_0 \f$ of a tracked signal.
 *
 * Implements a modification of the estimator presented in [1]. In [1] the
 * estimator essentially uses a moving average over the reciprocal of the
 * Signal-to-Noise Ratio (SNR). To reduce memory utilisation a simple IIR
 * low-pass filter is used instead.
 *
 * The noise and signal powers estimates for the \f$k\f$-th observation,
 * \f$\hat P_{N, k}\f$ and \f$\hat P_{S, k}\f$, are calculated as follows:
 *
 * \f[
 *    \hat P_{N, k} = \left( \left| I_k \right| -
 *                    \left| I_{k-1} \right| \right)^2
 * \f]
 * \f[
 *    \hat P_{S, k} = \frac{1}{2} \left( I_k^2 + I_{k-1}^2 \right)
 * \f]
 *
 * Where \f$I_k\f$ is the in-phase output of the prompt correlator for the
 * \f$k\f$-th integration period.
 *
 * The "Noise-to-Signal Ratio" (NSR) is estimated and filtered with a simple
 * low-pass IIR filter:
 *
 * \f[
 *    {NSR}_k = A \frac{\hat P_{N, k}}{\hat P_{S, k}} + (1 - A) {NSR}_{k-1}
 * \f]
 *
 * Where the IIR filter coefficient, \f$A\f$ can be calculated in terms of a
 * cutoff frequency \f$f_c\f$ and the loop update frequency \f$f = 1/T\f$.
 *
 * \f[
 *    A = \frac{f_c}{f_c + f}
 * \f]
 *
 * The filtered NSR value is converted to a \f$ C / N_0 \f$ value and returned.
 *
 * \f[
 *    \left( \frac{C}{N_0} \right)_k =
 *      10 \log_{10} \left( \frac{B_{eqn}}{{NSR}_k} \right)
 * \f]
 *
 * References:
 *    -# "Comparison of Four SNR Estimators for QPSK Modulations",
 *       Norman C. Beaulieu, Andrew S. Toms, and David R. Pauluzzi (2000),
 *       IEEE Communications Letters, Vol. 4, No. 2
 *    -# "Are Carrier-to-Noise Algorithms Equivalent in All Situations?"
 *       Inside GNSS, Jan / Feb 2010.
 *
 * \param s The estimator state struct to initialise.
 * \param I The prompt in-phase correlation from the tracking correlators.
 * \return The Carrier-to-Noise Density, \f$ C / N_0 \f$, in dBHz.
 */
float cn0_est(cn0_est_state_t *s, float I)
{
  float P_n, P_s;

  if (s->I_prev_abs < 0.f) {
    /* This is the first iteration, just update the prev state. */
    s->I_prev_abs = fabs(I);
  } else {
    P_n = fabsf(I) - s->I_prev_abs;
    P_n = P_n*P_n;

    P_s = 0.5f*(I*I + s->I_prev_abs*s->I_prev_abs);

    s->I_prev_abs = fabsf(I);

    s->nsr = s->A * (P_n / P_s) + (1.f - s->A) * s->nsr;
  }

  return s->log_bw - 10.f*log10f(s->nsr);
}

void calc_navigation_measurement(u8 n_channels, channel_measurement_t meas[], navigation_measurement_t nav_meas[], double nav_time, ephemeris_t ephemerides[])
{
  channel_measurement_t* meas_ptrs[n_channels];
  navigation_measurement_t* nav_meas_ptrs[n_channels];
  ephemeris_t* ephemerides_ptrs[n_channels];

  for (u8 i=0; i<n_channels; i++) {
    meas_ptrs[i] = &meas[i];
    nav_meas_ptrs[i] = &nav_meas[i];
    ephemerides_ptrs[i] = &ephemerides[meas[i].prn];
  }

  calc_navigation_measurement_(n_channels, meas_ptrs, nav_meas_ptrs, nav_time, ephemerides_ptrs);
}

void calc_navigation_measurement_(u8 n_channels, channel_measurement_t* meas[], navigation_measurement_t* nav_meas[], double nav_time, ephemeris_t* ephemerides[])
{
  double TOTs[n_channels];
  double min_TOT = DBL_MAX;

  for (u8 i=0; i<n_channels; i++) {
    TOTs[i] = 1e-3 * meas[i]->time_of_week_ms;
    TOTs[i] += meas[i]->code_phase_chips / 1.023e6;
    TOTs[i] += (nav_time - meas[i]->receiver_time) * meas[i]->code_phase_rate / 1.023e6;

    /** \todo Handle GPS time properly here, e.g. week rollover */
    nav_meas[i]->tot.wn = ephemerides[i]->toe.wn;
    nav_meas[i]->tot.tow = TOTs[i];

    if (gpsdifftime(nav_meas[i]->tot, ephemerides[i]->toe) > 3*24*3600)
      nav_meas[i]->tot.wn -= 1;

    if (TOTs[i] < min_TOT)
      min_TOT = TOTs[i];

    nav_meas[i]->raw_doppler = meas[i]->carrier_freq;
    nav_meas[i]->snr = meas[i]->snr;
    nav_meas[i]->prn = meas[i]->prn;

    nav_meas[i]->carrier_phase = meas[i]->carrier_phase;
    nav_meas[i]->carrier_phase += (nav_time - meas[i]->receiver_time) * meas[i]->carrier_freq;
  }

  double clock_err, clock_rate_err;

  for (u8 i=0; i<n_channels; i++) {
    nav_meas[i]->raw_pseudorange = (min_TOT - TOTs[i])*GPS_C + GPS_NOMINAL_RANGE;

    calc_sat_pos(nav_meas[i]->sat_pos, nav_meas[i]->sat_vel, &clock_err, &clock_rate_err, ephemerides[i], nav_meas[i]->tot);

    nav_meas[i]->pseudorange = nav_meas[i]->raw_pseudorange \
                               + clock_err*GPS_C;
    nav_meas[i]->doppler = nav_meas[i]->raw_doppler + clock_rate_err*GPS_L1_HZ;

    nav_meas[i]->tot.tow -= clock_err;
    nav_meas[i]->tot = normalize_gps_time(nav_meas[i]->tot);
  }
}

/** Compare navigation message by PRN.
 * This function is designed to be used together with qsort() etc.
 */
int nav_meas_cmp(const void *a, const void *b)
{
  return (s8)((navigation_measurement_t*)a)->prn
       - (s8)((navigation_measurement_t*)b)->prn;
}

/** Set measurement precise Doppler using time difference of carrier phase.
 * \note The return array `m_tdcp` should have space to contain the number
 * of measurements with common PRNs between `m_new` and `m_old`. Making the
 * array at least `MIN(n_new, n_old)` long will ensure sufficient space.
 *
 * \param n_new Number of measurements in `m_new`
 * \oaram m_new Array of new navigation measurements
 * \param n_old Number of measurements in `m_old`
 * \oaram m_new Array of old navigation measurements, sorted by PRN
 * \param m_tdcp Array in which to store the output measurements
 * \return The number of measurements written to `m_tdcp`
 */
u8 tdcp_doppler(u8 n_new, navigation_measurement_t *m_new,
                u8 n_old, navigation_measurement_t *m_old,
                navigation_measurement_t *m_corrected)
{
  /* Sort m_new, m_old should already be sorted. */
  qsort(m_new, n_new, sizeof(navigation_measurement_t), nav_meas_cmp);

  u8 i, j, n = 0;

  /* Loop over m_new and m_old and check if a PRN is present in both. */
  for (i=0, j=0; i<n_new && j<n_old; i++, j++) {
    if (m_new[i].prn < m_old[j].prn)
      j--;
    else if (m_new[i].prn > m_old[j].prn)
      i--;
    else {
      /* Copy m_new to m_corrected. */
      memcpy(&m_corrected[n], &m_new[i], sizeof(navigation_measurement_t));
      /* Calculate the Doppler correction between raw and corrected. */
      double dopp_corr = m_corrected[n].doppler - m_corrected[n].raw_doppler;
      /* Calculate raw Doppler from time difference of carrier phase. */
      /* TODO: check that using difference of TOTs here is a valid
       * approximation. */
      m_corrected[n].raw_doppler = (m_new[i].carrier_phase - m_old[j].carrier_phase)
                                    / gpsdifftime(m_new[i].tot, m_old[j].tot);
      /* Re-apply the same correction to the raw Doppler to get the corrected Doppler. */
      m_corrected[n].doppler = m_corrected[n].raw_doppler + dopp_corr;
      n++;
    }
  }

  return n;
}

/** \} */

