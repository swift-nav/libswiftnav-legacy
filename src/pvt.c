/*
 * Copyright (C) 2010 Swift Navigation Inc.
 * Contact: Henry Hallam <henry@swift-nav.com>
 *          Matt Peddie <peddie@alum.mit.edu>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "constants.h"
#include "linear_algebra.h"
#include "coord_system.h"
#include "track.h"

#include "pvt.h"

static double vel_solve(double rx_vel[],
                        const u8 n_used,
                        const navigation_measurement_t const nav_meas[n_used],
                        const double const G[n_used][4],
                        const double const X[4][n_used])
{
  /* Velocity Solution
   *
   * G and X matrices already exist from the position
   * solution loop through valid measurements.  Here we form satellite
   * velocity and pseudorange rate vectors -- it's the same
   * prediction-error least-squares thing, but we do only one step.
  */

  double tempvX[n_used];
  double pdot_pred;

  for (u8 j = 0; j < n_used; j++) {
    /* Calculate predicted pseudorange rates from the satellite velocity
     * and the geometry matix G which contains normalised line-of-sight
     * vectors to the satellites.
     */
    pdot_pred = -vector_dot(3, G[j], nav_meas[j].sat_vel);

    /* The residual is due to the user's motion. */
    tempvX[j] = -nav_meas[j].doppler * GPS_C / GPS_L1_HZ - pdot_pred;
  }

  /* Use X to map our pseudorange rate residuals onto the Jacobian update.
   *
   *   rx_vel[j] = X[j] . tempvX[j]
   */
  matrix_multiply(4, n_used, 1, (double *) X, (double *) tempvX, (double *) rx_vel);

  /* Return just the receiver clock bias. */
  return rx_vel[3];
}

void compute_dops(const double const H[4][4],
                  const double const pos_ecef[3],
                  dops_t *dops)
{
  double H_pos_diag[3];
  double H_ned[3];

  dops->gdop = dops->pdop = dops->tdop = dops->hdop = dops->vdop = 0;

  /* PDOP is the norm of the position elements of tr(H) */
  for (u8 i=0; i<3; i++) {
    dops->pdop += H[i][i];
    /* Also get the trace of H position states for use in HDOP/VDOP
     * calculations.
     */
    H_pos_diag[i] = H[i][i];
  }
  dops->pdop = sqrt(dops->pdop);

  /* TDOP is like PDOP but for the time state. */
  dops->tdop = sqrt(H[3][3]);

  /* Calculate the GDOP -- ||tr(H)|| = sqrt(PDOP^2 + TDOP^2) */
  dops->gdop = sqrt(dops->tdop*dops->tdop + dops->pdop*dops->pdop);

  /* HDOP and VDOP are Horizontal and Vertical; we need to rotate the
   * PDOP into NED frame and then take the separate components.
   */
  wgsecef2ned(H_pos_diag, pos_ecef, H_ned);
  dops->vdop = sqrt(H_ned[2]*H_ned[2]);
  dops->hdop = sqrt(H_ned[0]*H_ned[0] + H_ned[1]*H_ned[1]);
}


/* This function is the key to GPS solution, so it's commented
 * liberally.  It does a single step of a multi-dimensional
 * Newton-Raphson solution for the variables X, Y, Z (in ECEF) plus
 * the clock offset for each receiver used to make pseudorange
 * measurements.  The steps involved are roughly the following:
 *
 *     1. Account for the Earth's rotation during transmission
 *
 *     2. Estimate the ECEF position for each satellite measured using
 *     the downloaded ephemeris
 *
 *     3. Compute the Jacobian of pseudorange versus estimated state.
 *     There's no explicit differentiation; it's done symbolically
 *     first and just coded as a "line of sight" vector.
 *
 *     4. Get the inverse of the Jacobian times its transpose.  This
 *     matrix is normalized to one, but it tells us the direction we
 *     must move the state estimate during this step.
 *
 *     5. Multiply this inverse matrix (H) by the transpose of the
 *     Jacobian (to yield X).  This maps the direction of our state
 *     error into a direction of pseudorange error.
 *
 *     6. Multiply this matrix (X) by the error between the estimated
 *     (ephemeris) position and the measured pseudoranges.  This
 *     yields a vector of corrections to our state estimate.  We apply
 *     these to our current estimate and recurse to the next step.
 *
 *     7. If our corrections are very small, we've arrived at a good
 *     enough solution.  Solve for the receiver's velocity (with
 *     vel_solve) and do some bookkeeping to pass the solution back
 *     out.
 */
static double pvt_solve(double rx_state[],
                        const u8 n_used,
                        const navigation_measurement_t const nav_meas[n_used],
                        double H[4][4])
{
  double p_pred[n_used];

  /* Vector of prediction errors */
  double omp[n_used];

  /* G is a geometry matrix tells us how our pseudoranges relate to
   * our state estimates -- it's the Jacobian of d(p_i)/d(x_j) where
   * x_j are x, y, z, Δt. */
  double G[n_used][4];
  double Gtrans[4][n_used];
  double GtG[4][4];

  /* H is the square of the Jacobian matrix; it tells us the shape of
     our error (or, if you prefer, the direction in which we need to
     move to get a better solution) in terms of the receiver state. */

  /* X is just H * Gtrans -- it maps our pseudoranges onto our
   * Jacobian update */
  double X[4][n_used];

  double tempv[3];
  double los[3];
  double xk_new[3];
  double tempd;
  double correction[4];

  for (u8 j=0; j<4; j++) {
    correction[j] = 0.0;
  }

  for (u8 j = 0; j < n_used; j++) {
    /* The satellite positions need to be corrected for Earth's rotation during
     * the signal time of flight. */
    /* TODO: Explain more about how this corrects for the Sagnac effect. */

    /* Magnitude of range vector converted into an approximate time in secs. */
    vector_subtract(3, rx_state, nav_meas[j].sat_pos, tempv);
    double tau = vector_norm(3, tempv) / GPS_C;

    /* Rotation of Earth during time of flight in radians. */
    double wEtau = GPS_OMEGAE_DOT * tau;

    /* Apply linearlised rotation about Z-axis which will adjust for the
     * satellite's position at time t-tau. Note the rotation is through
     * -wEtau because it is the ECEF frame that is rotating with the Earth and
     * hence in the ECEF frame free falling bodies appear to rotate in the
     * opposite direction.
     *
     * Making a small angle approximation here leads to less than 1mm error in
     * the satellite position. */
    xk_new[0] = nav_meas[j].sat_pos[0] + wEtau * nav_meas[j].sat_pos[1];
    xk_new[1] = nav_meas[j].sat_pos[1] - wEtau * nav_meas[j].sat_pos[0];
    xk_new[2] = nav_meas[j].sat_pos[2];

    /* Line of sight vector. */
    vector_subtract(3, xk_new, rx_state, los);

    /* Predicted range from satellite position and estimated Rx position. */
    p_pred[j] = vector_norm(3, los);

    /* omp means "observed minus predicted" range -- this is E, the
     * prediction error vector (or innovation vector in Kalman/LS
     * filtering terms).
     */
    omp[j] = nav_meas[j].pseudorange - p_pred[j];

    /* Construct a geometry matrix.  Each row (satellite) is
     * independently normalized into a unit vector. */
    for (u8 i=0; i<3; i++) {
      G[j][i] = -los[i] / p_pred[j];
    }

    /* Set time covariance to 1. */
    G[j][3] = 1;

  } /* End of channel loop. */

  /* Solve for position corrections using batch least-squares.  When
   * all-at-once least-squares estimation for a nonlinear problem is
   * mixed with numerical iteration (not time-series recursion, but
   * iteration on a single set of measurements), it's basically
   * Newton's method.  There's a reasonably clear explanation of this
   * in Wikipedia's article on GPS.
   */

  /* Gt := G^{T} */
  matrix_transpose(n_used, 4, (double *) G, (double *) Gtrans);
  /* GtG := G^{T} G */
  matrix_multiply(4, n_used, 4, (double *) Gtrans, (double *) G, (double *) GtG);
  /* H \elem \mathbb{R}^{4 \times 4} := GtG^{-1} */
  matrix_inverse(4, (const double *) GtG, (double *) H);
  /* X := H * G^{T} */
  matrix_multiply(4, 4, n_used, (double *) H, (double *) Gtrans, (double *) X);
  /* correction := X * E (= X * omp) */
  matrix_multiply(4, n_used, 1, (double *) X, (double *) omp, (double *) correction);

  /* Increment ecef estimate by the new corrections */
  for (u8 i=0; i<3; i++) {
    rx_state[i] += correction[i];
  }

  /* Set the Δt estimates according to this solution */
  for (u8 i=3; i<4; i++) {
    rx_state[i] = correction[i];
  }

  /* Look at the magnintude of the correction to see if
   * the solution has converged yet.
   */
  tempd = vector_norm(3, correction);
  if(tempd > 0.001) {
    /* The solution has not converged, return a negative value to
     * indicate that we should continue iterating.
     */
    return -tempd;
  }

  /* The solution has converged! */

  /* Perform the velocity solution. */
  vel_solve(&rx_state[4], n_used, nav_meas, (const double (*)[4]) G, (const double (*)[n_used]) X);

  return tempd;
}

u8 filter_solution(gnss_solution* soln, dops_t* dops)
{
  if (dops->pdop > 50.0)
    /* PDOP is too high to yeild a good solution. */
    return 1;

  if (soln->pos_llh[2] < -1e3 || soln->pos_llh[2] > 1e8)
    /* Altitude is unreasonable. */
    return 2;

  /* NOTE: The following condition is required to comply with US export
   * regulations. It must not be removed. Any modification to this condition
   * is strictly not approved by Swift Navigation Inc. */

  if (soln->pos_llh[2] >= 0.3048*60000 &&
      vector_norm(3, soln->vel_ecef) >= 0.514444444*1000)
    /* Altitude is greater than 60000' and velocity is greater than 1000kts. */
    return 3;

  return 0;
}

u8 calc_PVT(const u8 n_used,
            const navigation_measurement_t const nav_meas[n_used],
            gnss_solution *soln,
            dops_t *dops)
{
  /* Initial state is the center of the Earth with zero velocity and zero
   * clock error, if we have some a priori position estimate we could use
   * that here to speed convergence a little on the first iteration.
   *
   *  rx_state format:
   *    pos[3], clock error, vel[3], intermediate freq error
   */
  static double rx_state[8];

  double H[4][4];

  soln->valid = 0;

  soln->n_used = n_used; // Keep track of number of working channels

  /* reset state to zero !? */
  for(u8 i=4; i<8; i++) {
    rx_state[i] = 0;
  }

  double update;
  u8 iters;
  /* Newton-Raphson iteration. */
  for (iters=0; iters<PVT_MAX_ITERATIONS; iters++) {
    if ((update = pvt_solve(rx_state, n_used, nav_meas, H)) > 0) {
      break;
    }
  }

  /* Compute various dilution of precision metrics. */
  compute_dops((const double(*)[4])H, rx_state, dops);
  soln->err_cov[6] = dops->gdop;

  /* Populate error covariances according to layout in definition
   * of gnss_solution struct.
   */
  soln->err_cov[0] = H[0][0];
  soln->err_cov[1] = H[0][1];
  soln->err_cov[2] = H[0][2];
  soln->err_cov[3] = H[1][1];
  soln->err_cov[4] = H[1][2];
  soln->err_cov[5] = H[2][2];

  if (iters >= PVT_MAX_ITERATIONS) {
    printf("Position solution not available after %d iterations, giving up.\n", iters);
    /* Reset state if solution fails */
    rx_state[0] = 0;
    rx_state[1] = 0;
    rx_state[2] = 0;
    return -1;
  }

  /* Save as x, y, z. */
  for (u8 i=0; i<3; i++) {
    soln->pos_ecef[i] = rx_state[i];
    soln->vel_ecef[i] = rx_state[4+i];
  }

  wgsecef2ned(soln->vel_ecef, soln->pos_ecef, soln->vel_ned);

  /* Convert to lat, lon, hgt. */
  wgsecef2llh(rx_state, soln->pos_llh);

  soln->clock_offset = rx_state[3] / GPS_C;
  soln->clock_bias = rx_state[7] / GPS_C;

  /* Time at receiver is TOT plus time of flight. Time of flight is eqaul to
   * the pseudorange minus the clock bias. */
  soln->time = nav_meas[0].tot;
  soln->time.tow += nav_meas[0].pseudorange / GPS_C;
  /* Subtract clock offset. */
  soln->time.tow -= rx_state[3] / GPS_C;
  soln->time = normalize_gps_time(soln->time);

  u8 ret;
  if ((ret = filter_solution(soln, dops))) {
    memset(soln, 0, sizeof(*soln));
    memset(dops, 0, sizeof(*dops));
    /* Reset state if solution fails */
    rx_state[0] = 0;
    rx_state[1] = 0;
    rx_state[2] = 0;
    printf("Solution filtered: %d\n", ret);
    return -1;
  }

  soln->valid = 1;
  return 0;
}

