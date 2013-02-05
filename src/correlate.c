/*
 * Copyright (C) 2013 Swift Navigation Inc.
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

#ifdef __SSSE3__
#include <tmmintrin.h>
#endif

#include "correlate.h"

/** \defgroup corr Correlation
 * Correlators used for tracking.
 * \{ */

#ifndef __SSSE3__

void track_correlate(s8* samples, s8* code,
                     double* init_code_phase, double code_step, double* init_carr_phase, double carr_step,
                     double* I_E, double* Q_E, double* I_P, double* Q_P, double* I_L, double* Q_L, u32* num_samples)
{
  double code_phase = *init_code_phase;
  double carr_phase = *init_carr_phase;

  double carr_sin = sin(carr_phase);
  double carr_cos = cos(carr_phase);
  double sin_delta = sin(carr_step);
  double cos_delta = cos(carr_step);

  *I_E = *Q_E = *I_P = *Q_P = *I_L = *Q_L = 0;

  double code_E, code_P, code_L;
  double baseband_Q, baseband_I;

  *num_samples = (int)ceil((1023.0 - code_phase) / code_step);

  for (u32 i=0; i<*num_samples; i++) {
    /*code_E = get_chip(code, (int)ceil(code_phase-0.5));*/
    /*code_P = get_chip(code, (int)ceil(code_phase));*/
    /*code_L = get_chip(code, (int)ceil(code_phase+0.5));*/
    /*code_E = code[(int)ceil(code_phase-0.5)];*/
    /*code_P = code[(int)ceil(code_phase)];*/
    /*code_L = code[(int)ceil(code_phase+0.5)];*/
    code_E = code[(int)(code_phase+0.5)];
    code_P = code[(int)(code_phase+1.0)];
    code_L = code[(int)(code_phase+1.5)];

    baseband_Q = carr_cos * samples[i];
    baseband_I = carr_sin * samples[i];

    double carr_sin_ = carr_sin*cos_delta + carr_cos*sin_delta;
    double carr_cos_ = carr_cos*cos_delta - carr_sin*sin_delta;
    double i_mag = (3.0 - carr_sin_*carr_sin_ - carr_cos_*carr_cos_) / 2.0;
    carr_sin = carr_sin_ * i_mag;
    carr_cos = carr_cos_ * i_mag;

    *I_E += code_E * baseband_I;
    *Q_E += code_E * baseband_Q;
    *I_P += code_P * baseband_I;
    *Q_P += code_P * baseband_Q;
    *I_L += code_L * baseband_I;
    *Q_L += code_L * baseband_Q;

    code_phase += code_step;
    carr_phase += carr_step;
  }
  *init_code_phase = code_phase - 1023;
  *init_carr_phase = fmod(carr_phase, 2*M_PI);
}

#else

void track_correlate(s8* samples, s8* code,
                     double* init_code_phase, double code_step, double* init_carr_phase, double carr_step,
                     double* I_E, double* Q_E, double* I_P, double* Q_P, double* I_L, double* Q_L, u32* num_samples)
{
  double code_phase = *init_code_phase;

  float carr_sin = sin(*init_carr_phase);
  float carr_cos = cos(*init_carr_phase);
  float sin_delta = sin(carr_step);
  float cos_delta = cos(carr_step);

  *num_samples = (int)ceil((1023.0 - code_phase) / code_step);

  __m128 IE_QE_IP_QP;
  __m128 CE_CE_CP_CP;
  __m128 IL_QL_X_X;
  __m128 CL_CL_X_X;
  __m128 S_C_S_C;
  __m128 BI_BQ_BI_BQ;
  __m128 dC_dS_dS_dC;
  __m128 a1, a2, a3;

  IE_QE_IP_QP = _mm_set_ps(0, 0, 0, 0);
  IL_QL_X_X = _mm_set_ps(0, 0, 0, 0);
  S_C_S_C = _mm_set_ps(carr_sin, carr_cos, carr_sin, carr_cos);
  dC_dS_dS_dC = _mm_set_ps(cos_delta, sin_delta, sin_delta, cos_delta);

  for (u32 i=0; i<*num_samples; i++) {
    CE_CE_CP_CP = _mm_set_ps(code[(int)(code_phase+0.5)],
                             code[(int)(code_phase+0.5)],
                             code[(int)(code_phase+1.0)],
                             code[(int)(code_phase+1.0)]);
    CL_CL_X_X = _mm_set_ps(code[(int)(code_phase+1.5)],
                           code[(int)(code_phase+1.5)],
                           0, 0);

    /* Load sample and multiply by sin/cos carrier to mix down to baseband. */
    a1 = _mm_set1_ps((float)samples[i]); // S, S, S, S
    BI_BQ_BI_BQ = _mm_mul_ps(a1, S_C_S_C);

    /* Update carrier sin/cos values by multiplying by the constant rotation
     * matrix corresponding to carr_step. */
    a1 = _mm_mul_ps(S_C_S_C, dC_dS_dS_dC); // SdC, CdS, SdS, CdC
    a2 = _mm_shuffle_ps(a1, a1, _MM_SHUFFLE(3, 0, 3, 0)); // SdC_CdC_SdC_CdC
    a3 = _mm_shuffle_ps(a1, a1, _MM_SHUFFLE(2, 1, 2, 1)); // CdS_SdS_CdS_SdS
    // S = SdC + CdS, C = CdC - SdS
    S_C_S_C = _mm_addsub_ps(a2, a3); // C_S_C_S

    /* Multiply code and baseband signal. */
    a1 = _mm_mul_ps(CE_CE_CP_CP, BI_BQ_BI_BQ);
    a2 = _mm_mul_ps(CL_CL_X_X, BI_BQ_BI_BQ);

    /* Increment accumulators. */
    IE_QE_IP_QP = _mm_add_ps(IE_QE_IP_QP, a1);
    IL_QL_X_X   = _mm_add_ps(IL_QL_X_X, a2);

    code_phase += code_step;
  }
  *init_code_phase = code_phase - 1023;
  *init_carr_phase = fmod(*init_carr_phase + *num_samples*carr_step, 2*M_PI);

  float res[8];
  _mm_storeu_ps(res, IE_QE_IP_QP);
  _mm_storeu_ps(res+4, IL_QL_X_X);

  *I_E = res[3];
  *Q_E = res[2];
  *I_P = res[1];
  *Q_P = res[0];
  *I_L = res[7];
  *Q_L = res[6];
}

#endif /* !__SSSE3__ */

/** \} */


