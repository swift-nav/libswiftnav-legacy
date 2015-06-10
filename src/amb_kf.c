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
 * It has been modified to be more robust against outliers.
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
#include "filter_utils.h"
#include "amb_kf.h"


/** \defgroup amb_kf Float Ambiguity Resolution
 * Preliminary integer ambiguity estimation with a Kalman Filter.
 * \{ */

/** Calculation of vectors needed for the innovation scaling.
 * We compute two vectors needed to make the Bierman update,
 * as well as the variance of the innovation.
 *
 * \param state_dim The dimension of the KF state.
 * \param h         A row of the observation matrix.
 * \param R         The variance of the observation.
 * \param U         KF state covariance U (not packed form)
 * \param D         KF state covariance D (stored as a vector)
 * \param f         U^T * h.
 * \param g         diag(D) * f.
 * \return          alpha = f * g + R, the innovation variance.
 */
double compute_innovation_terms(u32 state_dim, const double *h,
                                double R, const double *U,
                                const double *D, double *f, double *g)
{
  memcpy(f, h, state_dim * sizeof(double));
  /*  f = U^T * h. */
  cblas_dtrmv(CblasRowMajor, CblasUpper, CblasTrans, CblasUnit,
              /* ^ CBLAS_ORDER, CBLAS_UPLO, CBLAS_TRANSPOSE transA, CBLAS_DIAG. */
              state_dim, U,     /* int N, double *A. */
              state_dim, f, 1); /*  int lda, double *X, int incX. */


  /*  g = diag(D) * f.
      alpha = f * g + R = f^T * diag(D) * f + R. */
  double alpha = R;
  for (u32 i=0; i<state_dim; i++) {
    g[i] = D[i] * f[i];
    alpha += f[i] * g[i];
  }
  return alpha;
}


/** In place updating of the state cov and k vec using a scalar observation
 * This is from section 10.2.1 of Gibbs [1], with some extra logic for handling
 * singular matrices, dictating that zeros from cov_D dominate in a particular
 * potential 0 / 0.
 * We also make it more robust, by multiplying k by k_scalar <=  1.
 *
 * \param kf        The KF to update
 * \param R         The measurement variance
 * \param f         U^T * h
 * \param g         diag(D) * f
 * \param alpha     The innovation variance
 * \param k_scalar  A scalar to multiply the Kalman gain by (softens outliers).
 *                  Must be between 0 and 1 inclusive.
 * \param innov     The difference between the actual and predicted observation
 */
void update_kf_state(nkf_t *kf, double R, const double *f, const double *g,
                     double alpha, double k_scalar,
                     double innov)
{
  DEBUG_ENTRY();
  if (kf->state_dim == 0) {
    return;
  }
  /* If we are scaling the update by 0, we aren't updating at all,
   * so return early. */
  if (k_scalar == 0) {
    return;
  }
  assert(k_scalar <= 1);
  assert(k_scalar >= 0);
  u32 state_dim = kf->state_dim;
  double *U = kf->state_cov_U;
  double *D = kf->state_cov_D;
  double k[state_dim];

  double U_bar[state_dim * state_dim];
  double D_bar[state_dim];

  memset(U_bar, 0, state_dim * state_dim * sizeof(double));
  memset(D_bar, 0,             state_dim * sizeof(double));
  memset(k,     0,             state_dim * sizeof(double));

  /* K is inversely proportional to alpha, so we scale alpha to scale K.
   * Solving for an R that would give the properly scaled alpha and thus the
   * correct K, we get the following: */
  R += alpha * (1 - k_scalar) / k_scalar;
  double gamma = R  + g[0] * f[0];
  if (D[0] == 0 || R == 0) {
    /*  This is just an expansion of the other branch with the proper
     *  0 `div` 0 definitions. */
    D_bar[0] = 0;
  }
  else {
    D_bar[0] = D[0] * R / gamma;
  }
  k[0] = g[0];
  U_bar[0] = 1;
  if (DEBUG) {
    printf("gamma[0] = %f\n", gamma);
    printf("D_bar[0] = %f\n", D_bar[0]);
    VEC_PRINTF(k, state_dim);
    printf("U_bar[:,0] = {");
    for (u32 i=0; i < state_dim; i++) {
      printf("%f, ", U_bar[i*state_dim]);
    }
    printf("}\n");
  }
  for (u32 j=1; j<state_dim; j++) {
    double gamma_prev = gamma;
    gamma += g[j] * f[j];
    if (D[j] == 0 || gamma_prev == 0) {
      /* This is just an expansion of the other branch with the proper
       * 0 `div` 0 definitions. */
      D_bar[j] = 0;
    }
    else {
      D_bar[j] = D[j] * gamma_prev / gamma;
    }
    double f_over_gamma = f[j] / gamma_prev;
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
      printf("gamma[%"PRIu32"] = %f\n", j, gamma);
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

  /* Update the KF mean, scaled by some heuristic term for robustness */
  for (u32 j=0; j<kf->state_dim; j++) {
      kf->state_mean[j] += k[j] * k_scalar * innov;
  }
  if (DEBUG) {
    MAT_PRINTF(U, state_dim, state_dim);
    VEC_PRINTF(D, state_dim);
  }

  DEBUG_EXIT();
}

/** Get the weighted sum of squared innovations.
 * It's a normalized error metric. High is bad.
 * More precisely, we square the difference between the predicted observations
 * and the actual observations, divide them by their variances (HPH' + R),
 * but pretending they're independent, and then sum them.
 * It's: sum_i (z - H * x)_i / (H * P * H' + R)_ii.
 * If HPH'+R was actually diagonal and the innovations were independent
 * these would be chi square with kf->obs_dim degrees of freedom.
 *
 * If either dimension (state or observation) is zero, we define the output
 * to be zero.
 *
 * \param kf        Kalman filter struct
 * \param decor_obs Decorrelated observation vector
 * \return          The weighted SOS.
 */
double get_sos_innov(const nkf_t *kf, const double *decor_obs)
{
  assert(kf != NULL);
  assert(decor_obs != NULL);
  if (kf->state_dim == 0 || kf->obs_dim == 0) {
    return 0;
  }

  double predicted_obs[kf->obs_dim];
  matrix_multiply(kf->obs_dim, kf->state_dim, 1,
                  kf->decor_obs_mtx, kf->state_mean, predicted_obs);
  double hu[kf->obs_dim * kf->state_dim];
  /* TODO use the fact that U is unit triangular to save a ton of time */
  matrix_multiply(kf->obs_dim, kf->state_dim, kf->state_dim,
                  kf->decor_obs_mtx, kf->state_cov_U, hu);
  /* (H * U * D * U^T * H^T)_ii = (HU * D * HU^T)_ii
   *                            = Sum_kl (HU_ik * D_kl * HU^T_li)
   *                            = Sum_kl (HU_ik * D_kl * HU_il)
   *                            = Sum_k (HU_ik * D_kk * HU_ik) */
  double sos = 0;
  for (u8 i=0; i < kf->obs_dim; i++) {
    double hph_r_ii = kf->decor_obs_cov[i];
    for (u8 k=0; k < kf->state_dim; k++) {
      hph_r_ii += hu[i * kf->state_dim + k] *
                  hu[i * kf->state_dim + k] *
                  kf->state_cov_D[k];
    }
    sos += (predicted_obs[i] - decor_obs[i]) *
           (predicted_obs[i] - decor_obs[i]) /
           hph_r_ii;
  }
  return sos;
}

/** Compute a scale factor to soften outliers, updating an outlier filter.
 * We want to have some form of outlier detection that allows them
 * (especially near the edge of the classifier) to influence the filter.
 * This way, we soften the influence of any outliers, while eventually letting
 * change-points through.
 *
 * As of version 0.16 of piksi_firmware, this moving average in log space to
 * work best, only surpassed by a moving median. Both seem optimal with a
 * timescale of 7 observations, as per detecting change points, but if we
 * reduce the number of cycle slips in future releases, we may want to increase
 * it.
 *
 * \param kf        The Kalman filter.
 * \param decor_obs Decorrelated observation vector.
 * \param k_scalar  A scalar to multiply Kalman gain by.
 * \return          Whether we think the observation was an outlier (true=bad).
 */
bool outlier_check(nkf_t *kf, const double *decor_obs, double *k_scalar)
{
  assert (kf != NULL);
  assert (decor_obs != NULL);
  assert (k_scalar != NULL);
  if (kf->state_dim == 0 || kf->obs_dim == 0) {
    *k_scalar = 1;
    return false;
  }

  double sos = get_sos_innov(kf, decor_obs) / kf->obs_dim;
  double l_sos = log(MAX(1e-10,sos));
  double new_weight = 1.0f / KF_SOS_TIMESCALE;
  *k_scalar =  MIN(1, SOS_SWITCH * exp(kf->l_sos_avg - l_sos));
  /* l_sos_avg is a simple weighted average of its previous value and the
   * current l_sos:
   * kf->l_sos_avg = kf->l_sos_avg * (1 - new_weight) + l_sos * new_weight.
   * This can be simplified to the following: */
  kf->l_sos_avg += (l_sos - kf->l_sos_avg) * new_weight;
  return (*k_scalar < 1);
}

/** In place measurement update of the state mean and covariances.
 * A modification of the update in section 10.2.1 of Gibbs [1].
 * Changes are: -Scaling the kalman gain to deal with outliers.
 *              -Some minor tweaks to handle singular matrices.
 *
 * \param kf        Kalman filter to be updated.
 * \param decor_obs Decorrelated observation vector to be incorporated.
 * \return          Whether we think the observations were bad. (true = bad).
 */
static bool incorporate_obs(nkf_t *kf, double *decor_obs)
{
  DEBUG_ENTRY();

  double k_scalar;
  bool is_outlier = outlier_check(kf, decor_obs, &k_scalar);

  for (u32 i=0; i<kf->obs_dim; i++) {
    double *h = &kf->decor_obs_mtx[kf->state_dim * i]; /* vector of length kf->state_dim. */
    double R = kf->decor_obs_cov[i]; /* scalar. */

    double f[kf->state_dim];
    double g[kf->state_dim];

    double alpha = compute_innovation_terms(kf->state_dim, h, R,
                                            kf->state_cov_U, kf->state_cov_D,
                                            f, g);
    double predicted_obs = 0;
    /* TODO take advantage of sparsity of h. */
    for (u32 j=0; j<kf->state_dim; j++) {
      predicted_obs += h[j] * kf->state_mean[j];
    }
    double obs_minus_predicted_obs = decor_obs[i] - predicted_obs;

    /* updates kf state. */
    update_kf_state(kf, R, f, g, alpha, k_scalar, obs_minus_predicted_obs);
  }
  DEBUG_EXIT();
  return is_outlier;
}

/*  Turns (phi, rho) into Q_tilde * (phi, rho). */
static void make_residual_measurements(const nkf_t *kf, const double *measurements, double *resid_measurements)
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

/** The prediction step of the KF.
 * Since we're just doing parameter estimation, where we allow the parameters
 * to drift a bit (cycle slips and biases), we only update the covariances.
 *
 * To be really strict, since changes are equally likely on all channels,
 * This should be += (I + 1 * 1^T)*var instead of I*var, but it's
 * unlikely to be significant.
 *
 * \param kf The KF to be updated.
 */
static void diffuse_state(nkf_t *kf)
{
  double cov[kf->state_dim * kf->state_dim];
  matrix_reconstruct_udu(kf->state_dim, kf->state_cov_U, kf->state_cov_D, cov);
  for (u8 i=0; i< kf->state_dim; i++) {
    /* TODO make this a tunable parameter defined at the right time. */
    cov[i*kf->state_dim + i] += kf->amb_drift_var;
  }
  matrix_udu(kf->state_dim, cov, kf->state_cov_U, kf->state_cov_D);
}

/** In place updating of the KF state mean and covariance.
 * Does both the prediction and measurement updates to the KF.
 *
 * \param kf            The KF to update
 * \param measurements  The observations. The first (kf->state_dim) elements are
 *                      carrier phases, and the next (kf->state_dim) are
 *                      pseudoranges.
 * \return              Whether the KF thought the measurement was a bad
 *                      measurement. (true = bad)
 */
bool nkf_update(nkf_t *kf, const double *measurements)
{
  DEBUG_ENTRY();

  double resid_measurements[kf->obs_dim];
  make_residual_measurements(kf, measurements, resid_measurements);

  /* Replaces residual measurements by their decorrelated version. */
  cblas_dtrmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasUnit, /*  Order, Uplo, TransA, Diag. */
              kf->obs_dim, kf->decor_mtx, kf->obs_dim, /*  N, A, lda. */
              resid_measurements, 1); /*  X, incX. */

  /*  Prediction update */
  diffuse_state(kf);
  /* Measurement update */
  bool is_bad_measurement = incorporate_obs(kf, resid_measurements);

  if (DEBUG) {
    MAT_PRINTF(kf->state_cov_U, kf->state_dim, kf->state_dim);
    VEC_PRINTF(kf->state_cov_D, kf->state_dim);
    VEC_PRINTF(kf->state_mean, kf->state_dim);
  }
  DEBUG_EXIT();
  return is_bad_measurement;
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
  kf->l_sos_avg = 1;
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

/* Currently this function is only ever used to find the new position of a prn
 * in a new basis (permuted list) of prns, so all callsites assert the result
 * is not -1.
 */
s32 find_index_of_element_in_u8s(const u32 num_elements, const u8 x, const u8 *list)
{
  for (u32 i=0; i<num_elements; i++) {
    if (x == list[i]) {
      return i;
    }
  }
  assert(false);
}

/* REQUIRES num_sats > 1 */
void rebase_mean_N(double *mean, const u8 num_sats, const u8 *old_prns, const u8 *new_prns)
{
  assert(num_sats > 1);
  u8 state_dim = num_sats - 1;

  u8 old_ref = old_prns[0];
  u8 new_ref = new_prns[0];

  double new_mean[state_dim];
  s32 index_of_new_ref_in_old = find_index_of_element_in_u8s(num_sats-1, new_ref, &old_prns[1]);
  double val_for_new_ref_in_old_basis = mean[index_of_new_ref_in_old];
  for (u8 i=0; i<state_dim; i++) {
    u8 new_prn = new_prns[1+i];
    if (new_prn == old_ref) {
      new_mean[i] = - val_for_new_ref_in_old_basis;
    }
    else {
      s32 index_of_this_sat_in_old_basis = find_index_of_element_in_u8s(num_sats-1, new_prn, &old_prns[1]);
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

  s32 index_of_new_ref_in_old = find_index_of_element_in_u8s(num_sats-1, new_ref, &old_prns[1]);
  s32 index_of_old_ref_in_new = find_index_of_element_in_u8s(num_sats-1, old_ref, &new_prns[1]);
  for (u8 i=0; i<state_dim; i++) {
    rebase_mtx[i*state_dim + index_of_new_ref_in_old] = -1;
    if (i != (u8) index_of_old_ref_in_new) {
      s32 index_of_this_sat_in_old_basis = find_index_of_element_in_u8s(num_sats-1, new_prns[i+1], &old_prns[1]);
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
