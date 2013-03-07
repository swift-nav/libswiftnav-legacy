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

#include <math.h>

#include "pvt.h"
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
void calc_loop_gains(double bw, double zeta, double k, double loop_freq,
                     double *pgain, double *igain)
{
  /* Find the natural frequency. */
  double omega_n = bw*8*zeta / (4*zeta*zeta + 1);

  /* Some intermmediate values. */
/*
  double T = 1. / loop_freq;
  double denominator = k*(4 + 4*zeta*omega_n*T + omega_n*omega_n*T*T);

  *pgain = 8*zeta*omega_n*T / denominator;
  *igain = 4*omega_n*omega_n*T*T / denominator;
*/
  *igain = omega_n * omega_n / (k * loop_freq);
  *pgain = 2.0 * zeta * omega_n / k;
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
double costas_discriminator(double I, double Q)
{
  return atan(Q / I) / (2*M_PI);
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
double dll_discriminator(correlation_t cs[3])
{
  double early_mag = sqrt((double)cs[0].I*cs[0].I + (double)cs[0].Q*cs[0].Q);
  double late_mag = sqrt((double)cs[2].I*cs[2].I + (double)cs[2].Q*cs[2].Q);

  return 0.5 * (early_mag - late_mag) / (early_mag + late_mag);
}

/** Initialise a simple first-order loop filter.
 * The gains can be calculated using calc_loop_gains().
 *
 * \param s The loop filter state struct to initialise.
 * \param y0 The initial value of the output variable, \f$y_0\f$.
 * \param pgain The proportional gain, \f$k_p\f$.
 * \param igain The integral gain, \f$k_i\f$.
 */
void simple_lf_init(simple_lf_state_t *s, double y0,
                    double pgain, double igain)
{
  s->y = y0;
  s->prev_error = 0;
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
double simple_lf_update(simple_lf_state_t *s, double error)
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
void simple_tl_init(simple_tl_state_t *s, double loop_freq,
                    double code_freq, double code_bw,
                    double code_zeta, double code_k,
                    double carr_freq, double carr_bw,
                    double carr_zeta, double carr_k)
{
  double pgain, igain;

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
  double code_error = dll_discriminator(cs);
  s->code_freq = simple_lf_update(&(s->code_filt), -code_error);
  double carr_error = costas_discriminator(cs[1].I, cs[1].Q);
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
 * \param s The tracking loop state struct to initialise.
 * \param code_freq The initial code phase rate (i.e. frequency).
 * \param code_bw The code tracking loop noise bandwidth.
 * \param code_zeta The code tracking loop damping ratio.
 * \param code_k The code tracking loop gain.
 * \param carr_freq The initial carrier frequency.
 * \param carr_bw The carrier tracking loop noise bandwidth.
 * \param carr_zeta The carrier tracking loop damping ratio.
 * \param carr_k The carrier tracking loop gain.
 * \param tau The complimentary filter cross-over frequency.
 * \param sched The gain scheduling count.
 */
void comp_tl_init(comp_tl_state_t *s, double loop_freq,
                    double code_freq, double code_bw,
                    double code_zeta, double code_k,
                    double carr_freq, double carr_bw,
                    double carr_zeta, double carr_k,
                    double tau, u32 sched)
{
  double pgain, igain;

  calc_loop_gains(code_bw, code_zeta, code_k, loop_freq, &pgain, &igain);
  s->code_freq = code_freq;
  simple_lf_init(&(s->code_filt), code_freq, pgain, igain);

  calc_loop_gains(carr_bw, carr_zeta, carr_k, loop_freq, &pgain, &igain);
  s->carr_freq = carr_freq;
  simple_lf_init(&(s->carr_filt), carr_freq, pgain, igain);

  s->n = 0;
  s->sched = sched;

  s->A = 1.0 - (1.0 / (loop_freq * tau));
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
  double carr_error = costas_discriminator(cs[1].I, cs[1].Q);
  s->carr_freq = simple_lf_update(&(s->carr_filt), carr_error);

  double code_error = dll_discriminator(cs);
  s->code_filt.y = 0;
  double code_update = simple_lf_update(&(s->code_filt), -code_error);

  if (s->n > s->sched) {
    s->code_freq = s->A * s->code_freq + \
                   s->A * code_update + \
                   (1.0 - s->A)*
                     1.023e6*(1 + (s->carr_freq-4.092e6)/1.57542e9);
  } else {
    s->code_freq += code_update;
  }

  s->n++;
}


/** \} */

void calc_navigation_measurement(u8 n_channels, channel_measurement_t meas[], navigation_measurement_t nav_meas[], double nav_time, ephemeris_t ephemerides[])
{
  channel_measurement_t* meas_ptrs[n_channels];
  navigation_measurement_t* nav_meas_ptrs[n_channels];
  ephemeris_t* ephemerides_ptrs[n_channels];

  for (u8 i=0; i<n_channels; i++) {
    meas_ptrs[i] = &meas[i];
    nav_meas_ptrs[i] = &nav_meas[i];
    ephemerides_ptrs[i] = &ephemerides[i];
  }

  calc_navigation_measurement_(n_channels, meas_ptrs, nav_meas_ptrs, nav_time, ephemerides_ptrs);
}

void calc_navigation_measurement_(u8 n_channels, channel_measurement_t* meas[], navigation_measurement_t* nav_meas[], double nav_time, ephemeris_t* ephemerides[])
{
  double TOTs[n_channels];
  double mean_TOT = 0;

  for (u8 i=0; i<n_channels; i++) {
    TOTs[i] = 1e-3 * meas[i]->time_of_week_ms;
    TOTs[i] += meas[i]->code_phase_chips / 1.023e6;
    TOTs[i] += (nav_time - meas[i]->receiver_time) * meas[i]->code_phase_rate / 1.023e6;

    nav_meas[i]->TOT = TOTs[i];
    mean_TOT += TOTs[i];
    nav_meas[i]->pseudorange_rate = NAV_C * -meas[i]->carrier_freq / GPS_L1_HZ;
  }

  mean_TOT = mean_TOT/n_channels;

  double clock_err, clock_rate_err;

  double az, el;

  const double WPR_llh[3] = {D2R*37.038350, D2R*-122.141812, 376.7};
  double WPR_ecef[3];
  wgsllh2ecef(WPR_llh, WPR_ecef);

  for (u8 i=0; i<n_channels; i++) {
    nav_meas[i]->pseudorange = (mean_TOT - TOTs[i])*NAV_C + NOMINAL_RANGE;

    calc_sat_pos(nav_meas[i]->sat_pos, nav_meas[i]->sat_vel, &clock_err, &clock_rate_err, ephemerides[i], TOTs[i]);
    wgsecef2azel(nav_meas[i]->sat_pos, WPR_ecef, &az, &el);

    /*nav_meas[i]->pseudorange -= tropo_correction(el);*/
    nav_meas[i]->pseudorange += clock_err*NAV_C;
    nav_meas[i]->pseudorange_rate -= clock_rate_err*NAV_C;
  }
}

/** \} */

