/*
 * Copyright (C) 2014 Swift Navigation Inc.
 * Contact: Ian Horn <ian@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

/*
 * This is a Bierman-Thornton kalman filter implementation, as described in:
 *  [1] Gibbs, Bruce P. "Advanced Kalman Filtering, Least-Squares, and Modeling."
 *      John C. Wiley & Sons, Inc., 2011.
 */

#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <cblas.h>
#include <clapack.h>
#include <math.h>

#include "logging.h"
#include "linear_algebra.h"
#include "constants.h"
#include "track.h"
#include "almanac.h"
#include "gpstime.h"
#include "baseline.h"
#include "amb_kf.h"

/** Measure the integer ambiguity just from the code and carrier measurements.
 * The expectation value of carrier + code / lambda is
 * integer ambiguity + bias. Currently, pseudorange bias can get up to 10s of
 * wavelengths for minutes at a time, so averaging carrier + code isn't
 * sufficient for determining the ambiguity. Regardless of bias, this is an
 * important measurement. It is especially useful as a simple initialization
 * of the float filter.
 *
 * \param carrier A carrier phase measurement in units of wavelengths.
 * \param code    A code measurement in the same units as GPS_L1_LAMBDA_NO_VAC.
 * \return An estimate of the integer ambiguity. Its expectation value is the
 *         integer ambiguity plus carrier and code bias.
 */
double simple_amb_measurement(double carrier, double code)
{
  return carrier + code / GPS_L1_LAMBDA_NO_VAC;
}

/** \defgroup amb_kf Float Ambiguity Resolution
 * Preliminary integer ambiguity estimation with a Kalman Filter.
 * \{ */

/** In place updating of the state cov and k vec using a scalar observation
 * This is from section 10.2.1 of Gibbs [1], with some extra logic for handling
 *    singular matrices, dictating that zeros from cov_D dominate.
 *
 */
static void incorporate_scalar_measurement(u32 state_dim, double *h, double R,
                                           double *U, double *D, double *k)
{
  DEBUG_ENTRY();

  if (DEBUG) {
    VEC_PRINTF(h, state_dim);
    printf("R = %.16f", R);
    if (abs(R) == 0) {
      printf(" \t (R == 0 exactly)\n");
    }
    else {
      printf("\n");
    }
  }

  double f[state_dim]; /*  f = U^T * h. */
  memcpy(f, h, state_dim * sizeof(double));
  cblas_dtrmv(CblasRowMajor, CblasUpper, CblasTrans, CblasUnit,
              /* ^ CBLAS_ORDER, CBLAS_UPLO, CBLAS_TRANSPOSE transA, CBLAS_DIAG. */
              state_dim, U,     /* int N, double *A. */
              state_dim, f, 1); /*  int lda, double *X, int incX. */


  double g[state_dim]; /*  g = diag(D) * f. */
  double alpha = R;    /*  alpha = f * g + R = f^T * diag(D) * f + R. */
  for (u32 i=0; i<state_dim; i++) {
    g[i] = D[i] * f[i];
    alpha += f[i] * g[i];
  }
  if (DEBUG) {
    VEC_PRINTF(f, state_dim);
    VEC_PRINTF(g, state_dim);
    printf("alpha = %.16f", alpha);
    if (abs(alpha) == 0) {
      printf(" \t (alpha == 0 exactly)\n");
    }
    else {
      printf("\n");
    }
  }

  double gamma[state_dim];
  double U_bar[state_dim * state_dim];
  double D_bar[state_dim];

  memset(gamma, 0,             state_dim * sizeof(double));
  memset(U_bar, 0, state_dim * state_dim * sizeof(double));
  memset(D_bar, 0,             state_dim * sizeof(double));
  memset(k,     0,             state_dim * sizeof(double));

  gamma[0] = R + g[0] * f[0];
  if (D[0] == 0 || R == 0) {
    /*  This is just an expansion of the other branch with the proper
     *  0 `div` 0 definitions. */
    D_bar[0] = 0;
  }
  else {
    D_bar[0] = D[0] * R / gamma[0];
  }
  k[0] = g[0];
  U_bar[0] = 1;
  if (DEBUG) {
    printf("gamma[0] = %f\n", gamma[0]);
    printf("D_bar[0] = %f\n", D_bar[0]);
    VEC_PRINTF(k, state_dim);
    printf("U_bar[:,0] = {");
    for (u32 i=0; i < state_dim; i++) {
      printf("%f, ", U_bar[i*state_dim]);
    }
    printf("}\n");
  }
  for (u32 j=1; j<state_dim; j++) {
    gamma[j] = gamma[j-1] + g[j] * f[j];
    if (D[j] == 0 || gamma[j-1] == 0) {
      /* This is just an expansion of the other branch with the proper
       * 0 `div` 0 definitions. */
      D_bar[j] = 0;
    }
    else {
      D_bar[j] = D[j] * gamma[j-1] / gamma[j];
    }
    double f_over_gamma = f[j] / gamma[j-1];
    for (u32 i=0; i<=j; i++) {
      if (k[i] == 0) {
      /* This is just an expansion of the other branch with the proper
       * 0 `div` 0 definitions. */
        U_bar[i*state_dim + j] = U[i*state_dim + j];
      }
      else {
        /*  U_bar[:,j] = U[:,j] - f[j]/gamma[j-1] * k. */
        U_bar[i*state_dim + j] = U[i*state_dim + j] - f_over_gamma * k[i];
      }
      k[i] += g[j] * U[i*state_dim + j]; /*  k = k + g[j] * U[:,j]. */
    }
    if (DEBUG) {
      printf("gamma[%"PRIu32"] = %f\n", j, gamma[j]);
      printf("D_bar[%"PRIu32"] = %f\n", j, D_bar[j]);
      VEC_PRINTF(k, state_dim);
      printf("U_bar[:,%"PRIu32"] = {", j);
      for (u32 i=0; i < state_dim; i++) {
        printf("%f, ", U_bar[i*state_dim + j]);
      }
      printf("}\n");
    }
  }
  for (u32 i=0; i<state_dim; i++) {
    k[i] /= alpha;
  }
  memcpy(U, U_bar, state_dim * state_dim * sizeof(double));
  memcpy(D, D_bar,             state_dim * sizeof(double));
  if (DEBUG) {
    MAT_PRINTF(U, state_dim, state_dim);
    VEC_PRINTF(D, state_dim);
  }

  DEBUG_EXIT();
}

/** In place updating of the state mean and covariances to use the (decorrelated) observations
 * This is directly from section 10.2.1 of Gibbs [1]
 */
static void incorporate_obs(nkf_t *kf, double *decor_obs)
{
  DEBUG_ENTRY();

  for (u32 i=0; i<kf->obs_dim; i++) {
    double *h = &kf->decor_obs_mtx[kf->state_dim * i]; /* vector of length kf->state_dim. */
    double R = kf->decor_obs_cov[i]; /* scalar. */
    double k[kf->state_dim]; /*  vector of length kf->state_dim. */

    /* updates cov and sets k. */
    incorporate_scalar_measurement(kf->state_dim, h, R, kf->state_cov_U, kf->state_cov_D, &k[0]);

    double predicted_obs = 0;
    /* TODO take advantage of sparsity of h. */
    for (u32 j=0; j<kf->state_dim; j++) {
      predicted_obs += h[j] * kf->state_mean[j];
    }
    double obs_minus_predicted_obs = decor_obs[i] - predicted_obs;

    for (u32 j=0; j<kf->state_dim; j++) {
      kf->state_mean[j] += k[j] * obs_minus_predicted_obs; /* uses k to update mean. */
    }
  }

  DEBUG_EXIT();
}

/*  Turns (phi, rho) into Q_tilde * (phi, rho). */
static void make_residual_measurements(nkf_t *kf, double *measurements, double *resid_measurements)
{
  u8 constraint_dim = CLAMP_DIFF(kf->state_dim, 3);
  cblas_dgemv (CblasRowMajor, CblasNoTrans, /* Order, TransA. */
               constraint_dim, kf->state_dim, /*  M, N. */
               1, kf->null_basis_Q, kf->state_dim, /*  alpha A, lda. */
               measurements, 1, /*  X, incX. */
               0, resid_measurements, 1); /*  beta, Y, incY. */
  for (u8 i=0; i< kf->state_dim; i++) {
    resid_measurements[i+constraint_dim] =
      simple_amb_measurement(measurements[i],
                             measurements[i+kf->state_dim]);
  }
}

static void diffuse_state(nkf_t *kf)
{
  for (u8 i=0; i< kf->state_dim; i++) {
    /* TODO make this a tunable parameter defined at the right time. */
    kf->state_cov_D[i] += kf->amb_drift_var;
  }
}


/** In place updating of the state mean and covariance. Modifies measurements.
 */
void nkf_update(nkf_t *kf, double *measurements)
{
  DEBUG_ENTRY();

  double resid_measurements[kf->obs_dim];
  make_residual_measurements(kf, measurements, resid_measurements);

  /* Replaces residual measurements by their decorrelated version. */
  cblas_dtrmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasUnit, /*  Order, Uplo, TransA, Diag. */
              kf->obs_dim, kf->decor_mtx, kf->obs_dim, /*  N, A, lda. */
              resid_measurements, 1); /*  X, incX. */

  /*  predict_forward(kf);. */
  diffuse_state(kf);
  incorporate_obs(kf, resid_measurements);

  if (DEBUG) {
    MAT_PRINTF(kf->state_cov_U, kf->state_dim, kf->state_dim);
    VEC_PRINTF(kf->state_cov_D, kf->state_dim);
    VEC_PRINTF(kf->state_mean, kf->state_dim);
  }

  DEBUG_EXIT();
}

/*  Presumes that the first alm entry is the reference sat. */
void assign_de_mtx(u8 num_sats, const sdiff_t *sats_with_ref_first,
                   const double ref_ecef[3], double *DE)
{
  DEBUG_ENTRY();

  if (DEBUG) {
    printf("num_sats = %u\nsdiff prns&positions = {\n", num_sats);
    for (u8 i=0; i < num_sats; i++) {
      printf("i = %u, prn = %u, \tpos = [%f, \t%f, \t%f]\n",
             i,
             sats_with_ref_first[i].prn,
             sats_with_ref_first[i].sat_pos[0],
             sats_with_ref_first[i].sat_pos[1],
             sats_with_ref_first[i].sat_pos[2]);
    }
    printf("}\nref_ecef = {%f, \t%f, \t%f}\n", ref_ecef[0], ref_ecef[1], ref_ecef[2]);
  }

  assert(num_sats > 1);
  u8 de_length = num_sats - 1;

  if (num_sats <= 1) {
    log_debug("not enough sats\n");
    DEBUG_EXIT();
    return;
  }

  memset(DE, 0, de_length * 3 * sizeof(double));
  double e0[3];
  double x0 = sats_with_ref_first[0].sat_pos[0] - ref_ecef[0];
  double y0 = sats_with_ref_first[0].sat_pos[1] - ref_ecef[1];
  double z0 = sats_with_ref_first[0].sat_pos[2] - ref_ecef[2];
  double norm0 = sqrt(x0*x0 + y0*y0 + z0*z0);
  e0[0] = x0 / norm0;
  e0[1] = y0 / norm0;
  e0[2] = z0 / norm0;
  for (u8 i=1; i<num_sats; i++) {
    double x = sats_with_ref_first[i].sat_pos[0] - ref_ecef[0];
    double y = sats_with_ref_first[i].sat_pos[1] - ref_ecef[1];
    double z = sats_with_ref_first[i].sat_pos[2] - ref_ecef[2];
    double norm = sqrt(x*x + y*y + z*z);
    DE[3*(i-1)] = x / norm - e0[0];
    DE[3*(i-1) + 1] = y / norm - e0[1];
    DE[3*(i-1) + 2] = z / norm - e0[2];
  }
  if (DEBUG) {
    MAT_PRINTF(DE, (num_sats-1), 3);
  }
  DEBUG_EXIT();
}

void least_squares_solve_b(nkf_t *kf, const sdiff_t *sdiffs_with_ref_first,
                      const double *dd_measurements, const double ref_ecef[3],
                      double b[3])
{
  return least_squares_solve_b_external_ambs(kf->state_dim, kf->state_mean,
      sdiffs_with_ref_first, dd_measurements, ref_ecef, b);
}

/* Initializes the ambiguity means and variances.
 * Note that the covariance is  in UDU form, and U starts as identity. */
static void initialize_state(nkf_t *kf, double *dd_measurements, double init_var)
{
  u8 num_dds = kf->state_dim;
  for (u32 i=0; i<num_dds; i++) {
    kf->state_mean[i] =
      simple_amb_measurement(dd_measurements[i],
                             dd_measurements[i + num_dds]);
    /*  Sigma begins as a diagonal. */
    kf->state_cov_D[i] = init_var;
  }
  matrix_eye(num_dds, kf->state_cov_U);
}

static void QR_part1(integer m, integer n, double *A, double *tau)
{
  double w[1];
  integer lwork = -1;
  integer info;
  integer jpvt[3];
  memset(jpvt, 0, 3 * sizeof(integer));
  dgeqp3_(&m, &n,
          A, &m,
          jpvt,
          tau,
          w, &lwork, &info);
  lwork = round(w[0]);
  double work[lwork];
  dgeqp3_(&m, &n,
          A, &m,
          jpvt,
          tau,
          work, &lwork, &info); /* set A = QR(A). */
}

static void QR_part2(integer m, integer n, double *A, double *tau)
{
  double w[1];
  integer lwork = -1;
  integer info;
  dorgqr_(&m, &m, &n,
          A, &m,
          tau,
          w, &lwork, &info);
  lwork = round(w[0]);
  double work[lwork];
  dorgqr_(&m, &m, &n,
          A, &m,
          tau,
          work, &lwork, &info);
}

void assign_phase_obs_null_basis(u8 num_dds, double *DE_mtx, double *q)
{
  /* DE is num_sats-1 by 3, need to transpose it to column major. */
  double A[num_dds * num_dds];
  for (u8 i=0; i < num_dds; i++) {
    for (u8 j=0; j<3; j++) {
      A[j*num_dds + i] = DE_mtx[i*3 + j]; /* set A = Transpose(DE_mtx). */
    }
  }
  integer m = num_dds;
  integer n = 3;
  double tau[3];
  QR_part1(m, n, A, tau);
  QR_part2(m, n, A, tau);
  memcpy(q, &A[3*num_dds], CLAMP_DIFF(num_dds, 3) * num_dds * sizeof(double));
}

/* TODO this could be made more efficient, if it matters. */
static void assign_dd_obs_cov(u8 num_dds, double phase_var, double code_var,
                              double *dd_obs_cov)
{
  for (u8 i = 0; i < num_dds; i++) {
    for (u8 j = 0; j < num_dds; j++) {
      if (i == j) {
        dd_obs_cov[i*2*num_dds + j] = 2 * phase_var;
        dd_obs_cov[(i+num_dds)*2*num_dds + num_dds + j] = 2 * code_var;
      }
      else {
        dd_obs_cov[i*2*num_dds + j] = phase_var;
        dd_obs_cov[(i+num_dds)*2*num_dds + num_dds + j] = code_var;
      }
    }
    memset(&dd_obs_cov[i*2*num_dds + num_dds], 0, num_dds * sizeof(double));
    memset(&dd_obs_cov[(i+num_dds)*2*num_dds], 0, num_dds * sizeof(double));
  }
}

/* TODO make this more efficient (e.g. via pages 3/6.2-3/2014 of ian's notebook). */
static void assign_residual_obs_cov(u8 num_dds, double phase_var, double code_var,
                                    double *q, double *r_cov)
{
  double dd_obs_cov[4 * num_dds * num_dds];
  assign_dd_obs_cov(num_dds, phase_var, code_var, dd_obs_cov);
  integer nullspace_dim = CLAMP_DIFF(num_dds, 3);
  integer dd_dim = 2*num_dds;
  integer res_dim = num_dds + nullspace_dim;
  double q_tilde[res_dim * dd_dim];
  memset(q_tilde, 0, res_dim * dd_dim * sizeof(double));

  for (u8 i=0; i<nullspace_dim; i++) {
    memcpy(&q_tilde[i*dd_dim], &q[i*num_dds], num_dds * sizeof(double));
  }
  for (u8 i=0; i<num_dds; i++) {
    q_tilde[(i+nullspace_dim)*dd_dim + i] = 1;
    q_tilde[(i+nullspace_dim)*dd_dim + i+num_dds] = 1 / GPS_L1_LAMBDA_NO_VAC;
  }

  /* TODO make more efficient via the structure of q_tilde, and its relation to
   * the I + 1*1^T structure of the obs cov mtx. */
  double QC[res_dim * dd_dim];
  cblas_dsymm(CblasRowMajor, CblasRight, CblasUpper, /* CBLAS_ORDER, CBLAS_SIDE, CBLAS_UPLO. */
              res_dim, dd_dim,                       /* int M, int N. */
              1, dd_obs_cov, dd_dim,                 /* double alpha, double *A, int lda. */
              q_tilde, dd_dim,                       /* double *B, int ldb. */
              0, QC, dd_dim);                        /* double beta, double *C, int ldc. */

  /* TODO make more efficient via the structure of q_tilde, and its relation to
   * the I + 1*1^T structure of the obs cov mtx. */
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, /*  CBLAS_ORDER, CBLAS_TRANSPOSE transA, cBLAS_TRANSPOSE transB. */
              res_dim, res_dim, dd_dim,                /* int M, int N, int K. */
              1, QC, dd_dim,                           /* double alpha, double *A, int lda. */
              q_tilde, dd_dim,                         /* double *B, int ldb. */
              0, r_cov, res_dim);                      /* beta, double *C, int ldc. */
}

/*  In place inversion of U. */
static void invert_U(u8 res_dim, double *U)
{
  char uplo = 'U'; /* upper triangular. */
  char diag = 'U'; /* unit triangular. */
  integer dim = res_dim;
  integer lda = res_dim;
  integer info;
  dtrtri_(&uplo, &diag,
          &dim, U, &lda, &info);
}

static void assign_simple_sig(u8 num_dds, double var, double *simple_cov)
{
  for (u8 i = 0; i < num_dds; i++) {
    for (u8 j = 0; j < num_dds; j++) {
      if (i == j) {
        simple_cov[i*num_dds + j] = 2 * var;
      }
      else {
        simple_cov[i*num_dds + j] = var;
      }
    }
  }
}

/* From get_kf_matrices:
 * U^-1 * y == y'
 *           = U^-1 * H * x
 *          == H' * x
 *
 * H = ( Q )
 *     ( I )
 * where Q's rows form a basis for the left null space for DE
 */
static void assign_H_prime(u8 res_dim, u8 constraint_dim, u8 num_dds,
                           double *Q, double *U_inv, double *H_prime)
{
  /*  set the H_prime variable to equal H. */
  memcpy(H_prime, Q, constraint_dim * num_dds * sizeof(double));
  matrix_eye(num_dds, &H_prime[constraint_dim * num_dds]);

  /*  multiply H_prime by U_inv to make it the actual H_prime. */
  cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasUnit,
              /* ^ CBLAS_ORDER, CBLAS_SIDE, CBLAS_UPLO, CBLAS_TRANSPOSE, CBLAS_DIAG. */
              res_dim, num_dds,  /* M, N. */
              1, U_inv, res_dim, /* alpha, A, lda. */
              H_prime, num_dds); /* B, ldb. */
}

/* REQUIRES num_sdiffs > 0 */
/* y = H * x
 * Var[y] = Sig = U * D * U^T
 * ==> Var[U^-1 * y] = D
 * U^-1 * y == y'
 *           = U^-1 * H * x
 *          == H' * x
 *
 * H = ( Q )
 *     ( I )
 * where Q's rows form a basis for the left null space for DE
 *
 * Sig = Q~ * Sig_v * Q~^T
 *   where    Q~ = ( Q   0        )
 *                 ( I  -I/lambda )
 *     and Sig_v = ( D*D^T * var_phi   0               )
 *                 ( 0                 D*D^T * var_rho )
 *     and D*D^T = 1*1^T + I
 *
 * This function constructs D, U^-1, and H'
 */
static void get_kf_matrices(u8 num_sdiffs, sdiff_t *sdiffs_with_ref_first,
                            double ref_ecef[3],
                            double phase_var, double code_var,
                            double *null_basis_Q,
                            double *U_inv, double *D,
                            double *H_prime)
{
  assert (num_sdiffs > 0);

  u8 num_dds = num_sdiffs - 1;
  u8 constraint_dim = CLAMP_DIFF(num_dds, 3);
  u8 res_dim = num_dds + constraint_dim;

  double Sig[res_dim * res_dim];

  /* assign Sig and H. */
  if (constraint_dim > 0) {
    double DE[num_dds * 3];
    assign_de_mtx(num_sdiffs, sdiffs_with_ref_first, ref_ecef, DE);
    assign_phase_obs_null_basis(num_dds, DE, null_basis_Q);
    assign_residual_obs_cov(num_dds, phase_var, code_var, null_basis_Q, Sig);
    /* TODO U seems to have that fancy blockwise structure we love so much. Use it. */
    matrix_udu(res_dim, Sig, U_inv, D); /* U_inv holds U after this. */
    invert_U(res_dim, U_inv);
    /* TODO this also has fancy structure. */
    assign_H_prime(res_dim, constraint_dim, num_dds, null_basis_Q, U_inv, H_prime);
  }
  else {
    assign_simple_sig(num_dds,
                      phase_var + code_var / (GPS_L1_LAMBDA_NO_VAC * GPS_L1_LAMBDA_NO_VAC),
                      Sig);
    matrix_udu(res_dim, Sig, U_inv, D); /* U_inv holds U after this. */
    invert_U(res_dim, U_inv);

    /* H = I in this case, so H' = U^-1 * H = U^-1. */
    memcpy(H_prime, U_inv, num_dds * num_dds * sizeof(double));
  }

}


/* REQUIRES num_sats > 1 */
void set_nkf(nkf_t *kf, double amb_drift_var, double phase_var, double code_var, double amb_init_var,
            u8 num_sdiffs, sdiff_t *sdiffs_with_ref_first, double *dd_measurements, double ref_ecef[3])
{
  DEBUG_ENTRY();

  kf->amb_drift_var = amb_drift_var;
  set_nkf_matrices(kf, phase_var, code_var, num_sdiffs, sdiffs_with_ref_first, ref_ecef);
  /* Given plain old measurements, initialize the state. */
  initialize_state(kf, dd_measurements, amb_init_var);

  DEBUG_EXIT();
}

void set_nkf_matrices(nkf_t *kf, double phase_var, double code_var,
                     u8 num_sdiffs, sdiff_t *sdiffs_with_ref_first, double ref_ecef[3])
{
  assert(num_sdiffs > 1);

  u32 num_diffs = num_sdiffs - 1;
  kf->state_dim = num_sdiffs - 1;
  u32 constraint_dim = CLAMP_DIFF(num_diffs, 3);
  kf->obs_dim = num_diffs + constraint_dim;

  get_kf_matrices(num_sdiffs, sdiffs_with_ref_first,
                  ref_ecef,
                  phase_var, code_var,
                  kf->null_basis_Q,
                  kf->decor_mtx, kf->decor_obs_cov,
                  kf->decor_obs_mtx);
}

s32 find_index_of_element_in_u8s(const u32 num_elements, const u8 x, const u8 *list)
{
  for (u32 i=0; i<num_elements; i++) {
    if (x == list[i]) {
      return i;
    }
  }
  return -1;
}

/* REQUIRES num_sats > 1 */
void rebase_mean_N(double *mean, const u8 num_sats, const u8 *old_prns, const u8 *new_prns)
{
  assert(num_sats > 1);
  u8 state_dim = num_sats - 1;

  u8 old_ref = old_prns[0];
  u8 new_ref = new_prns[0];

  double new_mean[state_dim];
  s32 index_of_new_ref_in_old = find_index_of_element_in_u8s(num_sats, new_ref, &old_prns[1]);
  double val_for_new_ref_in_old_basis = mean[index_of_new_ref_in_old];
  for (u8 i=0; i<state_dim; i++) {
    u8 new_prn = new_prns[1+i];
    if (new_prn == old_ref) {
      new_mean[i] = - val_for_new_ref_in_old_basis;
    }
    else {
      s32 index_of_this_sat_in_old_basis = find_index_of_element_in_u8s(num_sats, new_prn, &old_prns[1]);
      new_mean[i] = mean[index_of_this_sat_in_old_basis] - val_for_new_ref_in_old_basis;
    }
  }
  memcpy(mean, new_mean, (state_dim) * sizeof(double));
}

/* REQUIRES num_sats > 1 */
static void assign_state_rebase_mtx(const u8 num_sats, const u8 *old_prns,
                                    const u8 *new_prns, double *rebase_mtx)
{
  assert(num_sats > 1);
  u8 state_dim = num_sats - 1;

  memset(rebase_mtx, 0, state_dim * state_dim * sizeof(double));
  u8 old_ref = old_prns[0];
  u8 new_ref = new_prns[0];

  s32 index_of_new_ref_in_old = find_index_of_element_in_u8s(state_dim, new_ref, &old_prns[1]);
  s32 index_of_old_ref_in_new = find_index_of_element_in_u8s(state_dim, old_ref, &new_prns[1]);
  for (u8 i=0; i<state_dim; i++) {
    rebase_mtx[i*state_dim + index_of_new_ref_in_old] = -1;
    if (i != (u8) index_of_old_ref_in_new) {
      s32 index_of_this_sat_in_old_basis = find_index_of_element_in_u8s(state_dim, new_prns[i+1], &old_prns[1]);
      rebase_mtx[i*state_dim + index_of_this_sat_in_old_basis] = 1;
    }
  }
}

/* REQUIRES num_sats > 1 */
void rebase_covariance_sigma(double *state_cov, const u8 num_sats, const u8 *old_prns, const u8 *new_prns)
{
  assert(num_sats > 1);
  u8 state_dim = num_sats - 1;

  double rebase_mtx[state_dim * state_dim];
  assign_state_rebase_mtx(num_sats, old_prns, new_prns, rebase_mtx);

  double intermediate_cov[state_dim * state_dim];
  /* TODO make more efficient via structure of rebase_mtx. */
  cblas_dsymm(CblasRowMajor, CblasRight, CblasUpper, /* CBLAS_ORDER, CBLAS_SIDE, CBLAS_UPLO. */
              state_dim, state_dim,                  /* int M, int N. */
              1, state_cov, state_dim,               /* double alpha, double *A, int lda. */
              rebase_mtx, state_dim,                 /* double *B, int ldb. */
              0, intermediate_cov, state_dim);       /* double beta, double *C, int ldc. */

  /* TODO make more efficient via the structure of rebase_mtx. */
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, /* CBLAS_ORDER, CBLAS_TRANSPOSE transA, cBLAS_TRANSPOSE transB. */
              state_dim, state_dim, state_dim,         /* int M, int N, int K. */
              1, intermediate_cov, state_dim,          /* double alpha, double *A, int lda. */
              rebase_mtx, state_dim,                   /* double *B, int ldb. */
              0, state_cov, state_dim);                /* beta, double *C, int ldc. */
}

/* REQUIRES num_sats > 1 */
void rebase_covariance_udu(double *state_cov_U, double *state_cov_D, u8 num_sats, u8 *old_prns, u8 *new_prns)
{
  assert(num_sats > 1);
  u8 state_dim = num_sats - 1;

  double state_cov[state_dim * state_dim];
  matrix_reconstruct_udu(state_dim, state_cov_U, state_cov_D, state_cov);
  rebase_covariance_sigma(state_cov, num_sats, old_prns, new_prns);
  matrix_udu(state_dim, state_cov, state_cov_U, state_cov_D);
}


/* REQUIRES num_sats > 1 */
void rebase_nkf(nkf_t *kf, u8 num_sats, u8 *old_prns, u8 *new_prns)
{
  assert(num_sats > 1);
  rebase_mean_N(kf->state_mean, num_sats, old_prns, new_prns);
  rebase_covariance_udu(kf->state_cov_U, kf->state_cov_D, num_sats, old_prns, new_prns);
}

void nkf_state_projection(nkf_t *kf,
                                    u8 num_old_non_ref_sats,
                                    u8 num_new_non_ref_sats,
                                    u8 *ndx_of_new_sat_in_old)
{
  u8 old_state_dim = num_old_non_ref_sats;
  double old_cov[old_state_dim * old_state_dim];
  matrix_reconstruct_udu(old_state_dim, kf->state_cov_U, kf->state_cov_D, old_cov);

  u8 new_state_dim = num_new_non_ref_sats;
  double new_cov[new_state_dim * new_state_dim];
  double new_mean[new_state_dim];

  for (u8 i=0; i<num_new_non_ref_sats; i++) {
    u8 ndxi = ndx_of_new_sat_in_old[i];
    new_mean[i] = kf->state_mean[ndxi];
    for (u8 j=0; j<num_new_non_ref_sats; j++) {
      u8 ndxj = ndx_of_new_sat_in_old[j];
      new_cov[i*new_state_dim + j] = old_cov[ndxi*old_state_dim + ndxj];
    }
  }

  /* Put it all back into the kf. */
  memcpy(kf->state_mean, new_mean, new_state_dim * sizeof(double));
  matrix_udu(new_state_dim, new_cov, kf->state_cov_U, kf->state_cov_D);
  /* NOTE: IT DOESN'T UPDATE THE OBSERVATION OR TRANSITION MATRICES, JUST THE STATE. */
}

/** Add new sats to the Kalman Filter
 * Given some space Z = X x Y, and some state mean/cov on X,
 * we construct a state mean/cov on Z by doing the inclusion of X
 * in Z and making initial estimates for the state of the Y elements
 * of Z. Here, X is the space of old satellite DDs, and Y is new sats.
 *
 * \param kf                    The KF struct to be updated
 * \param num_old_non_ref_sats  The old number of DDs (dimension of X)
 * \param num_new_non_ref_sats  The new number of DDs (dimension of Z)
 * \param ndx_of_old_sat_in_new The indices of the old sats in the new full
 *                              space. For example, if the first (ndx=0) element
 *                              of the state vector in X was for PRN 3 and
 *                              PRN 3 is the 4th (ndx=3) element of the state
 *                              vector in Z, then ndx_of_old_sat_in_new[0] = 3.
 * \param init_amb_est          Estimates of the ambiguities (in Z) used to
 *                              initialize the ambiguities of the new sats (Y).
 *                              Has length num_new_non_ref_sats.
 * \param int_init_var          Each element of init_amb_est should have
 *                              variance equal to int_init_var.
 */
void nkf_state_inclusion(nkf_t *kf,
                         u8 num_old_non_ref_sats,
                         u8 num_new_non_ref_sats,
                         u8 *ndx_of_old_sat_in_new,
                         double *init_amb_est,
                         double int_init_var)
{
  u8 old_state_dim = num_old_non_ref_sats;
  double old_cov[old_state_dim * old_state_dim];
  matrix_reconstruct_udu(old_state_dim, kf->state_cov_U, kf->state_cov_D, old_cov);

  u8 new_state_dim = num_new_non_ref_sats;
  double new_cov[new_state_dim * new_state_dim];
  memset(new_cov, 0, new_state_dim * new_state_dim * sizeof(double));
  double new_mean[new_state_dim];
  /* Initialize the ambiguity means/vars, including estimates for new sats. */
  memcpy(new_mean, init_amb_est, new_state_dim * sizeof(double));
  for (u8 i=0; i<num_new_non_ref_sats; i++) {
    new_cov[i*new_state_dim + i] = int_init_var;
  }

  /* Overwrite the ambiguity means/covars for the sats we were already tracking
   */
  for (u8 i=0; i<num_old_non_ref_sats; i++) {
    u8 ndxi = ndx_of_old_sat_in_new[i];
    new_mean[ndxi] = kf->state_mean[i];
    for (u8 j=0;j<num_old_non_ref_sats; j++) {
      u8 ndxj = ndx_of_old_sat_in_new[j];
      new_cov[ndxi*new_state_dim + ndxj] = old_cov[i*old_state_dim + j];
    }
  }
  matrix_udu(new_state_dim, new_cov, kf->state_cov_U, kf->state_cov_D);
  memcpy(kf->state_mean, new_mean, new_state_dim * sizeof(double));
}

/** \} */
