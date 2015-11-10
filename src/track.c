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
 * A second order PLL consists of a phase detector, first-order filter and a
 * Numerically Controlled Oscillator (NCO). A linearised model of a
 * second order PLL is shown below:
 *
 * \image html 2nd_order_dpll.png Linearised second order PLL model.
 *
 * Where \f$K_d\f$ is the discriminator gain and \f$F(s)\f$ and \f$N(s)\f$ are
 * the loop filter and VCO transfer functions. The VCO is essentially an
 * integrator and hence has a transfer function:
 *
 * \f[
 *   N(s) = \frac{K_0}{s}
 * \f]
 *
 * The first-order loop filter is of the form:
 * \f[
 *   F(s) = \frac{s\tau_2 + 1}{s\tau_1}
 * \f]
 * Where \f$\tau_1\f$ and \f$\tau_2\f$ are time constants.
 *
 * The closed loop transfer function for the above PLL is:
 * \f[
 *   H(s) = \frac{k_0\tau_2 s + k_0}
 *               {s^2 + \frac{k_0\tau_2}{\tau_1} s + \frac{k_0}{\tau_1}}
 * \f]
 *
 * Comparing the denominator above to the standard form for a second order
 * system \f$s^2 + 2\omega_n\zeta s + \omega_n^2\f$, we find find the
 * critical frequency and damping ration:
 * \f[
 *   \omega_n = \sqrt{\frac{k_0}{\tau_1}} \\
 *   \zeta = \frac{\omega_n\tau_2}{2}
 * \f]
 *
 * To design our digital implementation of the loop filter, we apply the
 * bilinear transfrom to the loop filter transfer function above
 * \f[
 *   H(z) = \left[\frac{k_0\tau_2 s + k_0}
 *               {s^2 + \frac{k_0\tau_2}{\tau_1} s + \frac{k_0}{\tau_1}}
 *          \right]_{s \leftarrow \frac{2(z-1)}{T(z + 1)}} \\
 *        = \frac{\frac{2\tau_2+T}{2\tau_1} + \frac{T-2\tau_2}{2\tau_1}z^{-1}}
 *               {1 - z^{-1}}
 * \f]
 *
 * The denominator coefficients \f$a_n\f$ are undependent of the loop
 * characteristics.  The numerator coefficients \f$b_n\f$ given by
 * \f[
 *   b_0 = \frac{2\tau_2 + T}{2\tau_1} \\
 *   b_1 = \frac{T - 2\tau_2}{2\tau_1}
 * \f]
 * where
 * \f[
 *   \tau_1 = \frac{k_0}{\omega_n^2} \\
 *   \tau_2 = \frac{2\zeta}{\omega_n}
 * \f]
 *
 * \param bw        The loop noise bandwidth, \f$B_L\f$.
 * \param zeta      The damping ratio, \f$\zeta\f$.
 * \param k         The loop gain, \f$k\f$.
 * \param loop_freq The loop update frequency, \f$1/T\f$.
 * \param b0        First filter coefficient, \f$b_0\f$.
 * \param b1        Second filter coefficient, \f$b_1\f$.
 */
void calc_loop_gains(float bw, float zeta, float k, float loop_freq,
                     float *b0, float *b1)
{
  float T = 1/loop_freq;
  /* Find the natural frequency. */
  float omega_n = 8.f*bw*zeta / (4.f*zeta*zeta + 1.f);

  /* Some intermmediate values. */
  float tau_1 = k / (omega_n * omega_n);
  float tau_2 = 2 * zeta / omega_n;

  *b0 = (2*tau_2 + T) / (2*tau_1);
  *b1 = (T - 2*tau_2) / (2*tau_1);
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
 * \param I The prompt in-phase correlation, \f$I_k\f$.
 * \param Q The prompt quadrature correlation, \f$Q_k\f$.
 * \return The discriminator value, \f$\varepsilon_k\f$.
 */
float costas_discriminator(float I, float Q)
{
  if (I == 0) {
    // Technically, it should be +/- 0.25, but then we'd have to keep track
    //  of the previous sign do it right, so it's simple enough to just return
    //  the average of 0.25 and -0.25 in the face of that ambiguity, so zero.
    return 0;
  }
  return atanf(Q / I) * (float)(1/(2*M_PI));
}

/** Frequency discriminator for a FLL, used to aid the PLL.
 *
 * Implements the \f$atan2\f$ frequency discriminator.
 *
 * Strictly speaking this should have a 1 / dt factor, but then we
 * would later multiply the gain by dt again, and it would cancel out.
 * So why do those unnecessary operations at all?
 *
 * \f[
 *    dot_k = abs(I_k * I_{k-1}) + abs(Q_k * Q_{k-1})\\
 *    cross_k = I_{k-1} * Q_k - I_k * Q_{k-1}\\
 *    \varepsilon_k = atan2 \left(cross_k, dot_k\right) / \pi
 * \f]
 *
 * References:
 *  -# Understanding GPS: Principles and Applications.
 *     Elliott D. Kaplan. Artech House, 1996.
 *
 * \param I The prompt in-phase correlation, \f$I_k\f$.
 * \param Q The prompt quadrature correlation, \f$Q_k\f$.
 * \param prev_I The prompt in-phase correlation, \f$I_{k-1}\f$.
 * \param prev_Q The prompt quadrature correlation, \f$Q_{k-1}\f$.
 * \return The discriminator value, \f$\varepsilon_k\f$.
 */
float frequency_discriminator(float I, float Q, float prev_I, float prev_Q)
{
  float dot = fabsf(I * prev_I) + fabsf(Q * prev_Q);
  float cross = prev_I * Q - I * prev_Q;
  return atan2f(cross, dot) / ((float) M_PI);
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

/** Initialize an integral aided loop filter.
 *
 * This initializes a feedback loop with a PI component, plus an extra independent I term.
 *
 * \param s The loop filter state struct to initialize.
 * \param y0 The initial value of the output variable, \f$y_0\f$.
 * \param b0 The proportional gain of the PI error term, \f$b_0\f$.
 * \param b1 The integral gain of the PI error term, \f$b_1\f$.
 * \param aiding_igain The integral gain of the aiding error term, \f$k_{ia}\f$.
 */
void aided_lf_init(aided_lf_state_t *s, float y0,
                   float b0, float b1,
                   float aiding_igain)
{
  s->y = y0;
  s->prev_error = 0.f;
  s->b0 = b0;
  s->b1 = b1;
  s->aiding_igain = aiding_igain;
}

/** Update step for the integral aided loop filter.
 *
 * This updates a feedback loop with a PI component plus an extra independent I term.
 *
 * \param s The loop filter state struct.
 * \param p_i_error The error output from the discriminator used in both P and I terms.
 * \param aiding_error The error output from the discriminator use just in an I term.
 * \return The updated output variable.
 */
float aided_lf_update(aided_lf_state_t *s, float p_i_error, float aiding_error)
{
  s->y += (s->b0 * p_i_error) + (s->b1 * s->prev_error) +
          s->aiding_igain * aiding_error;
  s->prev_error = p_i_error;

  return s->y;
}

/** Initialise a simple first-order loop filter.
 * The gains can be calculated using calc_loop_gains().
 *
 * \param s The loop filter state struct to initialise.
 * \param y0 The initial value of the output variable, \f$y_0\f$.
 * \param b0 First filter coefficient, \f$b_0\f$.
 * \param b1 Second filter coefficient, \f$b_1\f$.
 */
void simple_lf_init(simple_lf_state_t *s, float y0,
                    float b0, float b1)
{
  s->y = y0;
  s->prev_error = 0.f;
  s->b0 = b0;
  s->b1 = b1;
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
 *   F[z] = \frac{b_0 + b_1 z^{-1}}{1 - z^{-1}}
 * \f]
 *
 * \param s The loop filter state struct.
 * \param error The error output from the discriminator, \f$x_n\f$.
 * \return The updated output variable, \f$y_n\f$.
 */
float simple_lf_update(simple_lf_state_t *s, float error)
{
  s->y += (s->b0 * error) + (s->b1 * s->prev_error);
  s->prev_error = error;

  return s->y;
}

/** Initialise an aided tracking loop.
 *
 * For a full description of the loop filter parameters, see calc_loop_gains().
 *
 * \param s The tracking loop state struct to initialise.
 * \param loop_freq The loop update rate.
 * \param code_freq The initial code phase rate (i.e. frequency).
 * \param code_bw The code tracking loop noise bandwidth.
 * \param code_zeta The code tracking loop damping ratio.
 * \param code_k The code tracking loop gain.
 * \param carr_to_code Ratio of carrier to code frequencies,
 *                     or 0 to disable carrier aiding.
 * \param carr_freq The initial carrier frequency.
 * \param carr_bw The carrier tracking loop noise bandwidth.
 * \param carr_zeta The carrier tracking loop damping ratio.
 * \param carr_k The carrier tracking loop gain.
 * \param carr_freq_b1 The integral gain of the aiding error term, \f$k_{ia}\f$.
 */
void aided_tl_init(aided_tl_state_t *s, float loop_freq,
                   float code_freq,
                   float code_bw, float code_zeta, float code_k,
                   float carr_to_code,
                   float carr_freq,
                   float carr_bw, float carr_zeta, float carr_k,
                   float carr_freq_b1)
{
  float b0, b1;

  s->carr_freq = carr_freq;
  s->prev_I = 1.0f; // This works, but is it a really good way to do it?
  s->prev_Q = 0.0f;
  calc_loop_gains(carr_bw, carr_zeta, carr_k, loop_freq, &b0, &b1);
  aided_lf_init(&(s->carr_filt), carr_freq, b0, b1, carr_freq_b1);

  calc_loop_gains(code_bw, code_zeta, code_k, loop_freq, &b0, &b1);
  s->code_freq = code_freq;
  s->carr_to_code = carr_to_code;
  /* If using carrier aiding, initialize code_freq in code loop filter
     to zero to avoid double-counting. */
  simple_lf_init(&(s->code_filt), carr_to_code ? 0 : code_freq, b0, b1);
}

void aided_tl_retune(aided_tl_state_t *s, float loop_freq,
                     float code_bw, float code_zeta, float code_k,
                     float carr_to_code,
                     float carr_bw, float carr_zeta, float carr_k,
                     float carr_freq_b1)
{
  calc_loop_gains(carr_bw, carr_zeta, carr_k, loop_freq,
                  &s->carr_filt.b0, &s->carr_filt.b1);
  s->carr_filt.aiding_igain = carr_freq_b1;

  calc_loop_gains(code_bw, code_zeta, code_k, loop_freq,
                  &s->code_filt.b0, &s->code_filt.b1);
  s->carr_to_code = carr_to_code;
}

/** Update step for the aided tracking loop.
 *
 * Implements a basic second-order tracking loop. The code tracking loop is a
 * second-order DLL using dll_discriminator() as its discriminator function.
 * The carrier phase tracking loop is a second-order Costas loop using
 * costas_discriminator(), aided by a frequency discriminator using
 * frequency_discriminator().
 *
 * The tracking loop output variables, i.e. code and carrier frequencies can be
 * read out directly from the state struct.
 *
 * \param s The tracking loop state struct.
 * \param cs An array [E, P, L] of correlation_t structs for the Early, Prompt
 *           and Late correlations.
 */
void aided_tl_update(aided_tl_state_t *s, correlation_t cs[3])
{
  /* Carrier loop */
  float carr_error = costas_discriminator(cs[1].I, cs[1].Q);
  float freq_error = 0;
  if (s->carr_filt.aiding_igain != 0) { /* Don't waste cycles if not using FLL */
    freq_error = frequency_discriminator(cs[1].I, cs[1].Q, s->prev_I, s->prev_Q);
    s->prev_I = cs[1].I;
    s->prev_Q = cs[1].Q;
  }
  s->carr_freq = aided_lf_update(&(s->carr_filt), carr_error, freq_error);

  /* Code loop */
  float code_error = dll_discriminator(cs);
  s->code_freq = simple_lf_update(&(s->code_filt), -code_error);
  if (s->carr_to_code) /* Optional carrier aiding of code loop */
    s->code_freq += s->carr_freq / s->carr_to_code;
}

/** Initialise a simple tracking loop.
 *
 * For a full description of the loop filter parameters, see calc_loop_gains().
 *
 * \param s The tracking loop state struct to initialise.
 * \param loop_freq The loop update frequency, \f$1/T\f$.
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
  float b0, b1;

  calc_loop_gains(code_bw, code_zeta, code_k, loop_freq, &b0, &b1);
  s->code_freq = code_freq;
  simple_lf_init(&(s->code_filt), code_freq, b0, b1);

  calc_loop_gains(carr_bw, carr_zeta, carr_k, loop_freq, &b0, &b1);
  s->carr_freq = carr_freq;
  simple_lf_init(&(s->carr_filt), carr_freq, b0, b1);
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
 * \param loop_freq The loop update frequency, \f$1/T\f$.
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
  float b0, b1;

  calc_loop_gains(code_bw, code_zeta, code_k, loop_freq, &b0, &b1);
  s->code_freq = code_freq;
  simple_lf_init(&(s->code_filt), code_freq, b0, b1);

  calc_loop_gains(carr_bw, carr_zeta, carr_k, loop_freq, &b0, &b1);
  s->carr_freq = carr_freq;
  simple_lf_init(&(s->carr_filt), carr_freq, b0, b1);

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

/** Initialise alias lock detection structure.
 *
 * \param a         The alias lock detect struct.
 * \param acc_len   The number of points to accumulate before calculating error.
 * \param time_diff Time difference between calls to alias_detect_first and
 *                  alias_detect_full.
 */
void alias_detect_init(alias_detect_t *a, u32 acc_len, float time_diff)
{
  memset(a, 0, sizeof(*a));
  a->acc_len = acc_len;
  a->dt = time_diff;
}

/** Load first I/Q sample into alias lock detect structure.
 *
 * \param a The alias lock detect struct.
 * \param I The prompt in-phase correlation.
 * \param Q The prompt quadrature-phase correlation.
 */
void alias_detect_first(alias_detect_t *a, float I, float Q)
{
  a->first_I = I;
  a->first_Q = Q;
}

/** Load second I/Q sample into alias lock detect structure and return
 *  frequency error if ready.
 *
 * \param a The alias lock detect struct.
 * \param I The prompt in-phase correlation.
 * \param Q The prompt quadrature-phase correlation.
 * \returns Calculated frequency error or zero.
 */
float alias_detect_second(alias_detect_t *a, float I, float Q)
{
  a->dot += (I * a->first_I + Q * a->first_Q) / a->acc_len;
  a->cross += (a->first_I * Q - I * a->first_Q) / a->acc_len;
  a->fl_count++;
  if (a->fl_count == a->acc_len) {
    float err = atan2f(a->cross, a->dot) / ((float)M_PI * a->dt);
    a->fl_count = 0;
    a->cross = 0;
    a->dot = 0;
    return err;
  }
  return 0;
}

/** Initialise the lock detector state.
 * \param l
 * \param k1 LPF coefficient.
 * \param k2 I Scale factor.
 * \param lp Pessimistic count threshold.
 * \param lo Optimistic count threshold
 */
void lock_detect_init(lock_detect_t *l, float k1, float k2, u16 lp, u16 lo)
{
  memset(l, 0, sizeof(*l)); /* sets l->lpf[iq].y = 0 */
  lock_detect_reinit(l, k1, k2, lp, lo);
}

/** Update the lock detector parameters, preserving internal state.
 * \param l
 * \param k1 LPF coefficient.
 * \param k2 I Scale factor.
 * \param lp Pessimistic count threshold.
 * \param lo Optimistic count threshold
 */
void lock_detect_reinit(lock_detect_t *l, float k1, float k2, u16 lp, u16 lo)
{
  /* Adjust LPF coefficient */
  l->lpfi.k1 = k1;
  l->lpfq.k1 = k1;

  l->k2 = k2;
  l->lp = lp;
  l->lo = lo;
}

static float lock_detect_lpf_update(struct loop_detect_lpf *lpf, float x)
{
  lpf->y += lpf->k1 * (x - lpf->y);
  return lpf->y;
}

/** Update the lock detector with new prompt correlations.
 * \param l
 * \param I In-phase prompt correlation.
 * \param Q Quadrature prompt correlation.
 * \param DT Integration time
 *
 * References:
 *  -# Understanding GPS: Principles and Applications, 2nd Edition
 *     Section 5.11.2, pp 233-235
 *     Elliott D. Kaplan. Artech House, 1996.
 */
void lock_detect_update(lock_detect_t *l, float I, float Q, float DT)
{
  float a, b;
  /* Calculated low-pass filtered prompt correlations */
  a = lock_detect_lpf_update(&l->lpfi, fabs(I) / DT) / l->k2;
  b = lock_detect_lpf_update(&l->lpfq, fabs(Q) / DT);

  if (a > b) {
    /* In-phase > quadrature, looks like we're locked */
    l->outo = true;
    l->pcount2 = 0;
    /* Wait before raising the pessimistic indicator */
    if (l->pcount1 > l->lp) {
      l->outp = true;
    } else {
      l->pcount1++;
    }
  } else {
    /* In-phase < quadrature, looks like we're not locked */
    l->outp = false;
    l->pcount1 = 0;
    /* Wait before lowering the optimistic indicator */
    if (l->pcount2 > l->lo) {
      l->outo = false;
    } else {
      l->pcount2++;
    }
  }
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
  float Tw0 = (2*M_PI*cutoff_freq) / loop_freq;
  s->b = Tw0 / (Tw0 + 2);
  s->a = (Tw0 - 2) / (Tw0 + 2);

  s->log_bw = 10.f*log10f(bw);
  s->I_prev_abs = -1.f;
  s->Q_prev_abs = -1.f;
  s->nsr = powf(10.f, 0.1f*(s->log_bw - cn0_0));
}

/** Estimate the Carrier-to-Noise Density, \f$ C / N_0 \f$ of a tracked signal.
 *
 * Implements a modification of the estimator presented in [1]. In [1] the
 * estimator essentially uses a moving average over the reciprocal of the
 * Signal-to-Noise Ratio (SNR). To reduce memory utilisation a simple IIR
 * low-pass filter is used instead.
 *
 * The noise and signal power estimates for the \f$k\f$-th observation,
 * \f$\hat P_{N, k}\f$ and \f$\hat P_{S, k}\f$, are calculated as follows:
 *
 * \f[
 *    \hat P_{N, k} = \left( \left| Q_k \right| -
 *                    \left| Q_{k-1} \right| \right)^2
 * \f]
 * \f[
 *    \hat P_{S, k} = \frac{1}{2} \left( I_k^2 + I_{k-1}^2 \right)
 * \f]
 *
 * Where \f$I_k\f$ is the in-phase output of the prompt correlator for the
 * \f$k\f$-th integration period.
 *
 * The "Noise-to-Signal Ratio" (NSR) is estimated and filtered with a first
 * order low-pass IIR filter with transfer function:
 *
 * \f[
 *    F(s) = \frac{\omega_0}{s + \omega_0}
 * \f]
 * The bilinear transform is applied to obtain a digital equivalent
 * \f[
 *    F(z) = \frac{b + bz^{-1}}{1 + az^{-1}}
 * \f]
 * where \f$ b = \frac{T\omega_0}{T\omega_0 + 2} \f$
 * and \f$ a = \frac{T\omega_0 - 2}{T\omega_0 + 2} \f$.
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
 * \param I The prompt in-phase correlation from the tracking correlator.
 * \param Q The prompt quadrature correlation from the tracking correlator.
 * \return The Carrier-to-Noise Density, \f$ C / N_0 \f$, in dBHz.
 */
float cn0_est(cn0_est_state_t *s, float I, float Q)
{
  float P_n, P_s;

  if (s->I_prev_abs < 0.f) {
    /* This is the first iteration, just update the prev state. */
    s->I_prev_abs = fabsf(I);
    s->Q_prev_abs = fabsf(Q);
  } else {
    P_n = fabsf(Q) - s->Q_prev_abs;
    P_n = P_n*P_n;

    P_s = 0.5f*(I*I + s->I_prev_abs*s->I_prev_abs);

    s->I_prev_abs = fabsf(I);
    s->Q_prev_abs = fabsf(Q);

    float tmp = s->b * P_n / P_s;
    s->nsr = tmp + s->xn - s->a * s->nsr;
    s->xn = tmp;
  }

  return s->log_bw - 10.f*log10f(s->nsr);
}

void calc_navigation_measurement(u8 n_channels, channel_measurement_t meas[],
                                 navigation_measurement_t nav_meas[], double nav_time,
                                 ephemeris_t ephemerides[])
{
  channel_measurement_t* meas_ptrs[n_channels];
  navigation_measurement_t* nav_meas_ptrs[n_channels];
  ephemeris_t* ephemerides_ptrs[n_channels];

  for (u8 i=0; i<n_channels; i++) {
    meas_ptrs[i] = &meas[i];
    nav_meas_ptrs[i] = &nav_meas[i];
    ephemerides_ptrs[i] = &ephemerides[meas[i].sid.sat];
  }

  calc_navigation_measurement_(n_channels, meas_ptrs, nav_meas_ptrs, nav_time, ephemerides_ptrs);
}

void calc_navigation_measurement_(u8 n_channels, channel_measurement_t* meas[], navigation_measurement_t* nav_meas[], double nav_time, ephemeris_t* ephemerides[])
{
  double TOTs[n_channels];
  double min_TOF = -DBL_MAX;
  double clock_err[n_channels], clock_rate_err[n_channels];

  for (u8 i=0; i<n_channels; i++) {
    TOTs[i] = 1e-3 * meas[i]->time_of_week_ms;
    TOTs[i] += meas[i]->code_phase_chips / 1.023e6;
    TOTs[i] += (nav_time - meas[i]->receiver_time) * meas[i]->code_phase_rate / 1.023e6;

    /** \todo Maybe keep track of week number in tracking channel
        state or derive it from system time. */
    nav_meas[i]->tot.tow = TOTs[i];
    gps_time_match_weeks(&nav_meas[i]->tot, &ephemerides[i]->toe);

    nav_meas[i]->raw_doppler = meas[i]->carrier_freq;
    nav_meas[i]->snr = meas[i]->snr;
    nav_meas[i]->sid.sat = meas[i]->sid.sat;

    nav_meas[i]->carrier_phase = meas[i]->carrier_phase;
    nav_meas[i]->carrier_phase += (nav_time - meas[i]->receiver_time) * meas[i]->carrier_freq;

    nav_meas[i]->lock_counter = meas[i]->lock_counter;

    /* calc sat clock error */
    calc_sat_state(ephemerides[i], nav_meas[i]->tot,
                   nav_meas[i]->sat_pos, nav_meas[i]->sat_vel,
                   &clock_err[i], &clock_rate_err[i]);

    /* remove clock error to put all tots within the same time window */
    if ((TOTs[i] + clock_err[i]) > min_TOF)
      min_TOF = TOTs[i];
  }

  for (u8 i=0; i<n_channels; i++) {
    nav_meas[i]->raw_pseudorange = (min_TOF - TOTs[i])*GPS_C + GPS_NOMINAL_RANGE;

    nav_meas[i]->pseudorange = nav_meas[i]->raw_pseudorange \
                               + clock_err[i]*GPS_C;
    nav_meas[i]->doppler = nav_meas[i]->raw_doppler + clock_rate_err[i]*GPS_L1_HZ;

    nav_meas[i]->tot.tow -= clock_err[i];
    nav_meas[i]->tot = normalize_gps_time(nav_meas[i]->tot);
  }
}

/** Compare navigation message by PRN.
 * This function is designed to be used together with qsort() etc.
 */
int nav_meas_cmp(const void *a, const void *b)
{
  return sid_compare(((navigation_measurement_t*)a)->sid,
                     ((navigation_measurement_t*)b)->sid);
}

/** Set measurement precise Doppler using time difference of carrier phase.
 * \note The return array `m_tdcp` should have space to contain the number
 * of measurements with common PRNs between `m_new` and `m_old`. Making the
 * array at least `MIN(n_new, n_old)` long will ensure sufficient space.
 *
 * \param n_new Number of measurements in `m_new`
 * \param m_new Array of new navigation measurements
 * \param n_old Number of measurements in `m_old`
 * \param m_old Array of old navigation measurements, sorted by PRN
 * \param m_corrected Array in which to store the output measurements
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
    if (sid_compare(m_new[i].sid, m_old[j].sid) < 0)
      j--;
    else if (sid_compare(m_new[i].sid, m_old[j].sid) > 0)
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

