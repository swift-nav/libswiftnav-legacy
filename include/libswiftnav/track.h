/*
 * Copyright (C) 2012, 2016 Swift Navigation Inc.
 * Contact: Fergus Noble <fergus@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef LIBSWIFTNAV_TRACK_H
#define LIBSWIFTNAV_TRACK_H

#include <libswiftnav/common.h>
#include <libswiftnav/ephemeris.h>
#include <libswiftnav/signal.h>

/** \addtogroup track
 * \{ */

/** \addtogroup track_loop
 * \{ */

/** State structure for the I-aided loop filter.
 * Should be initialised with aided_lf_init().
 */
typedef struct {
  float b0;            /**< Filter coefficient. */
  float b1;            /**< Filter coefficient. */
  float aiding_igain;  /**< Aiding integral gain. */
  float prev_error;    /**< Previous error. */
  float y;             /**< Output variable. */
} aided_lf_state_t;

/** State structure for the simple loop filter.
 * Should be initialised with simple_lf_init().
 */
typedef struct {
  float b0;         /**< Filter coefficient. */
  float b1;         /**< Filter coefficient. */
  float prev_error; /**< Previous error. */
  float y;          /**< Output variable. */
} simple_lf_state_t;

typedef struct {
  float carr_freq;             /**< Code frequency. */
  aided_lf_state_t carr_filt ; /**< Carrier loop filter state. */
  float code_freq;             /**< Carrier frequency. */
  simple_lf_state_t code_filt; /**< Code loop filter state. */
  float prev_I, prev_Q;        /**< Previous timestep's in-phase and
				    quadrature integration, for FLL. */
  float carr_to_code;          /**< Ratio of code to carrier freqs, or
				    zero to disable carrier aiding */
} aided_tl_state_t;

/**
 * Tracking loop control object
 *
 * This object is made according to Elliott D.Kaplan and Chrostopher J.Hegarty
 * book with bilinear transform integrators.
 */
typedef struct {
  float T;

  /* First order frequency filter */
  float freq_omega_0;          /**< Frequency error multiplier: omega_0f */
  float freq_acc;              /**< Frequency error accumulator */

  /* Second order phase filter */
  float phase_omega_02;        /**< Phase error multiplier 1: omega_n^2 */
  float phase_omega_0a;        /**< Phase error multiplier 2: omega_n*a */
  float phase_acc2;            /**< Phase error accumulator 1 */
  float phase_acc;             /**< Phase error accumulator 2 */

  /* Doppler frequency */
  float carr_freq;             /**< Carrier frequency output. */

  float code_freq;             /**< Code frequency error. */
  simple_lf_state_t code_filt; /**< Code loop filter state. */
  float prev_I, prev_Q;        /**< Previous timestep's in-phase and
                                    quadrature integration, for FLL. */
  float carr_to_code;          /**< Ratio of code to carrier freqs, or
                                    zero to disable carrier aiding */
} aided_tl_state_new_t;

/** State structure for a simple tracking loop.
 * Should be initialised with simple_tl_init().
 */
typedef struct {
  float code_freq;             /**< Code phase rate (i.e. frequency). */
  float carr_freq;             /**< Carrier frequency. */
  simple_lf_state_t code_filt; /**< Code loop filter state. */
  simple_lf_state_t carr_filt; /**< Carrier loop filter state. */
} simple_tl_state_t;

/** State structure for a code/carrier phase complimentary filter tracking
 * loop.
 * Should be initialised with comp_tl_init().
 */
typedef struct {
  float code_freq;             /**< Code phase rate (i.e. frequency). */
  float carr_freq;             /**< Carrier frequency. */
  simple_lf_state_t code_filt; /**< Code loop filter state. */
  simple_lf_state_t carr_filt; /**< Carrier loop filter state. */
  u32 sched;                   /**< Gain scheduling count. */
  u32 n;                       /**< Iteration counter. */
  float A;                     /**< Complementary filter crossover gain. */
  float carr_to_code;          /**< Scale factor from carrier to code. */
} comp_tl_state_t;

typedef struct {
  u8 acc_len;     /**< Accumulation length. */
  float dt;       /**< Time difference between sample points. */
  float dot;        /**< Accumulated dot products. */
  float cross;      /**< Accumulated cross products. */
  u8 fl_count;    /**< Currently accumulated point count. */
  float first_I;    /**< First in-phase sample. */
  float first_Q;    /**< First quadrature-phase sample. */
} alias_detect_t;

struct loop_detect_lpf {
  float k1;                    /**< Filter coefficient. */
  float y;                     /**< Output and state variable. */
};

/** State structure for basic lock detector with optimistic and pessimistic
 *  indicators.
 */
typedef struct {
  struct loop_detect_lpf lpfi; /**< I path LPF state. */
  struct loop_detect_lpf lpfq; /**< Q path LPF state. */
  float k2;                    /**< I Scale factor. */
  u16 lo, lp;                  /**< Optimistic and pessimistic count threshold. */
  u16 pcount1, pcount2;        /**< Counter state variables. */
  bool outo, outp;             /**< Optimistic and pessimistic indicator. */
} lock_detect_t;

/** \} */

/** Structure representing a complex valued correlation. */
typedef struct {
  float I; /**< In-phase correlation. */
  float Q; /**< Quadrature correlation. */
} correlation_t;

/** State structure for the \f$ C / N_0 \f$ estimator.
 * Should be initialized with cn0_est_init().
 */
typedef struct {
  float log_bw;     /**< Noise bandwidth in dBHz. */
  float I_prev_abs; /**< Abs. value of the previous in-phase correlation. */
  float Q_prev_abs; /**< Abs. value of the previous quadrature correlation. */
  float cn0;        /**< Signal to noise ratio in dB/Hz. */
} cn0_est_state_t;

/** State structure for first order low-pass filter.
 *
 * \see lp1_filter_init
 * \see lp1_filter_update
 */
typedef struct {
  float b;          /**< IIR filter coeff. */
  float a;          /**< IIR filter coeff. */
  float xn;         /**< Last pre-filter sample. */
  float yn;         /**< Last post-filter sample. */
} lp1_filter_t;

/**
 * Second order Butterworth filter object.
 *
 * Structure for filtering CN0 values using 2nd order Butterworth filter.
 *
 * \see bw2_filter_init
 * \see bw2_filter_update
 */
typedef struct {
  float b;          /**< IIR filter coeff. */
  float a2;         /**< IIR filter coeff. */
  float a3;         /**< IIR filter coeff. */
  float yn;         /**< Last post-filter sample. */
  float yn_prev;    /**< Previous post-filter sample. */
  float xn;         /**< Last pre-filter sample. */
  float xn_prev;    /**< Previous pre-filter sample. */
} bw2_filter_t;

/** This struct holds the state of a tracking channel at a given receiver time epoch.
 *
 * The struct contains the information necessary to calculate the pseudorange,
 * carrier phase and Doppler information needed for a PVT solution but is
 * formatted closer to the natural outputs from the tracking channels.
 *
 * \see calc_navigation_measurement()
 */
typedef struct {
  gnss_signal_t sid;       /**< Satellite signal. */
  double code_phase_chips; /**< The code-phase in chips at `receiver_time`. */
  double code_phase_rate;  /**< Code phase rate in chips/s. */
  double carrier_phase;    /**< Carrier phase in cycles. */
  double carrier_freq;     /**< Carrier frequency in Hz. */
  u32 time_of_week_ms;     /**< Number of milliseconds since the start of the
                                GPS week corresponding to the last code rollover.  */
  double rec_time_delta;   /**< Difference between receiver clock time at which
                                this measurement is valid and reference time
                                (seconds). */
  double snr;              /**< Signal to noise ratio. */
  u16 lock_counter;        /**< This number is changed each time the tracking
                                channel is re-locked or a cycle slip is
                                detected and the carrier phase integer
                                ambiguity is reset.  If this number changes it
                                is an indication you should reset integer
                                ambiguity resolution for this channel. */
} channel_measurement_t;

typedef struct {
  double raw_pseudorange;
  double pseudorange;
  double raw_carrier_phase;
  double carrier_phase;
  double raw_doppler;
  double doppler;
  double sat_pos[3];
  double sat_vel[3];
  double snr;
  double lock_time;
  gps_time_t tot;
  gnss_signal_t sid;
  u16 lock_counter;
} navigation_measurement_t;

/** \} */

void calc_loop_gains(float bw, float zeta, float k, float loop_freq,
                     float *b0, float *b1);
float costas_discriminator(float I, float Q);
float frequency_discriminator(float I, float Q, float prev_I, float prev_Q);
float dll_discriminator(correlation_t cs[3]);

void aided_lf_init(aided_lf_state_t *s, float y0,
                   float pgain, float igain,
                   float aiding_igain);
float aided_lf_update(aided_lf_state_t *s, float p_i_error, float aiding_error);

void simple_lf_init(simple_lf_state_t *s, float y0,
                    float pgain, float igain);
float simple_lf_update(simple_lf_state_t *s, float error);


void simple_tl_init(simple_tl_state_t *s, float loop_freq,
                    float code_freq, float code_bw,
                    float code_zeta, float code_k,
                    float carr_freq, float carr_bw,
                    float carr_zeta, float carr_k);
void simple_tl_update(simple_tl_state_t *s, correlation_t cs[3]);

void aided_tl_init(aided_tl_state_t *s, float loop_freq,
                   float code_freq,
                   float code_bw, float code_zeta, float code_k,
                   float carr_to_code,
                   float carr_freq,
                   float carr_bw, float carr_zeta, float carr_k,
                   float carr_freq_b1);

void aided_tl_retune(aided_tl_state_t *s, float loop_freq,
                     float code_bw, float code_zeta, float code_k,
                     float carr_to_code,
                     float carr_bw, float carr_zeta, float carr_k,
                     float carr_freq_b1);

void aided_tl_update(aided_tl_state_t *s, correlation_t cs[3]);

void aided_tl_init_new(aided_tl_state_new_t *s, float loop_freq,
                   float code_freq,
                   float code_bw, float code_zeta, float code_k,
                   float carr_to_code,
                   float carr_freq,
                   float carr_bw, float carr_zeta, float carr_k,
                   float carr_freq_b1);

void aided_tl_retune_new(aided_tl_state_new_t *s, float loop_freq,
                     float code_bw, float code_zeta, float code_k,
                     float carr_to_code,
                     float carr_bw, float carr_zeta, float carr_k,
                     float carr_freq_b1);

void aided_tl_update_new(aided_tl_state_new_t *s, correlation_t cs[3]);

void comp_tl_init(comp_tl_state_t *s, float loop_freq,
                    float code_freq, float code_bw,
                    float code_zeta, float code_k,
                    float carr_freq, float carr_bw,
                    float carr_zeta, float carr_k,
                    float tau, float cpc, u32 sched);
void comp_tl_update(comp_tl_state_t *s, correlation_t cs[3]);

void alias_detect_init(alias_detect_t *a, u32 acc_len, float time_diff);
void alias_detect_reinit(alias_detect_t *a, u32 acc_len, float time_diff);
void alias_detect_first(alias_detect_t *a, float I, float Q);
float alias_detect_second(alias_detect_t *a, float I, float Q);

void lock_detect_init(lock_detect_t *l, float k1, float k2, u16 lp, u16 lo);
void lock_detect_reinit(lock_detect_t *l, float k1, float k2, u16 lp, u16 lo);
void lock_detect_update(lock_detect_t *l, float I, float Q, float DT);

void cn0_est_bl_init(cn0_est_state_t *s,
                     float bw, float cn0_0, float f_s, float f_i);
float cn0_est_bl_update(cn0_est_state_t *s, float I, float Q);
void cn0_est_snv_init(cn0_est_state_t *s,
                      float cn0_0, float bw, float f_s, float f_i);
float cn0_est_snv_update(cn0_est_state_t *s, float I, float Q);

void lp1_filter_init(lp1_filter_t *f, float initial,
                     float cutoff_freq, float loop_freq);
float lp1_filter_update(lp1_filter_t *f, float value);

void bw2_filter_init(bw2_filter_t *f, float initial,
                     float cutoff_freq, float loop_freq);
float bw2_filter_update(bw2_filter_t *f, float value);

s8 calc_navigation_measurement(u8 n_channels, const channel_measurement_t *meas[],
                               navigation_measurement_t *nav_meas[], gps_time_t *rec_time,
                               const ephemeris_t* e[]);

int nav_meas_cmp(const void *a, const void *b);
u8 tdcp_doppler(u8 n_new, navigation_measurement_t *m_new,
                u8 n_old, navigation_measurement_t *m_old,
                navigation_measurement_t *m_corrected, double dt);

#endif /* LIBSWIFTNAV_TRACK_H */
