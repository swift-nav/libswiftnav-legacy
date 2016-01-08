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


#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include <libswiftnav/constants.h>
#include <libswiftnav/logging.h>
#include <libswiftnav/linear_algebra.h>
#include <libswiftnav/coord_system.h>
#include <libswiftnav/track.h>
#include <libswiftnav/pvt.h>

#include <libswiftnav/plover/pvt.h>

static void compute_dops(const double H[4][4],
                         const double pos_ecef[3],
                         dops_t *dops)
{
  /* PDOP is the norm of the position elements of tr(H) */
  double pdop_sq = H[0][0] + H[1][1] + H[2][2];
  dops->pdop = sqrt(pdop_sq);

  /* TDOP is like PDOP but for the time state. */
  dops->tdop = sqrt(H[3][3]);

  /* Calculate the GDOP -- ||tr(H)|| = sqrt(PDOP^2 + TDOP^2) */
  dops->gdop = sqrt(pdop_sq + H[3][3]);

  /* HDOP and VDOP are Horizontal and Vertical.  We could rotate H
   * into NED frame and then take the separate components, but a more
   * computationally efficient approach is to find the vector in the
   * ECEF frame that represents the Down unit vector, and project it
   * through H.  That gives us VDOP^2, then we find HDOP from the
   * relation PDOP^2 = HDOP^2 + VDOP^2. */
  double M[3][3];
  ecef2ned_matrix(pos_ecef, M);
  double down_ecef[4] = {M[2][0], M[2][1], M[2][2], 0};
  double tmp[3];
  matrix_multiply(3, 4, 1, (double *)H, down_ecef, tmp);
  double vdop_sq = vector_dot(3, down_ecef, tmp);
  dops->vdop = sqrt(vdop_sq);
  dops->hdop = sqrt(pdop_sq - vdop_sq);
}


/** This function is the key to GPS solution, so it's commented
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

/* ! see pvt.plv ! */

static u8 filter_solution(gnss_solution* soln, dops_t* dops)
{
  if (dops->pdop > 50.0)
    /* PDOP is too high to yield a good solution. */
    return 1;

  if (soln->pos_llh[2] < -1e3 || soln->pos_llh[2] > 1e6)
    /* Altitude is unreasonable. */
    return 2;

  /* NOTE: The following condition is required to comply with US export
   * regulations. It must not be removed. Any modification to this condition
   * is strictly not approved by Swift Navigation Inc. */

  if (vector_norm(3, soln->vel_ecef) >= 0.514444444*1000)
    /* Velocity is greater than 1000kts. */
    return 3;

  return 0;
}

/** Checks pvt_iter residuals.
 *
 * \param n_used length of omp
 * \param omp residual vector calculated by pvt_solve
 * \param rx_state reference to pvt solver state allocated in calc_PVT
 * \param residual If not null, used to output double value of residual
 *
 * \return residual < PVT_RESIDUAL_THRESHOLD
 */
static bool residual_test(u8 n_used, double omp[n_used],
                          const double rx_state[],
                          double *residual)
{
  /* Very liberal threshold. Typical range 20 - 120 */
#define PVT_RESIDUAL_THRESHOLD 3000

  /* Need to add clock offset to observed-minus-predicted calculated by last
   * iteration of pvt_solve before computing residual. */
  for (int i = 0; i < n_used; i++) {
    omp[i] -= rx_state[3];
  }
  double norm = vector_norm(n_used, omp);
  if (residual) {
    *residual = norm;
  }
  return norm < PVT_RESIDUAL_THRESHOLD;
}

/** Iterates pvt_solve until it converges or PVT_MAX_ITERATIONS is reached.
 *
 * \return
 *   - `0`: solution converged
 *   - `-1`: solution failed to converge
 *
 *  Results stored in rx_state, omp, H
 */
static s8 pvt_iter(double rx_state[],
                   const u8 n_used,
                   const navigation_measurement_t *nav_meas[n_used],
                   double omp[n_used],
                   double H[4][4])
{
  /* Reset state to zero */
  for(u8 i=4; i<8; i++) {
    rx_state[i] = 0;
  }

  u8 iters;
  /* Newton-Raphson iteration. */
  for (iters=0; iters<PVT_MAX_ITERATIONS; iters++) {
    if (calc_pvt(rx_state, n_used, (const navigation_measurement_t * *)nav_meas, omp, (double *)H) > 0) {
      break;
    }
  }

  if (iters >= PVT_MAX_ITERATIONS) {
    /* Reset state if solution fails */
    rx_state[0] = 0;
    rx_state[1] = 0;
    rx_state[2] = 0;
    return -1;
  }

  return 0;
}

/** See pvt_solve_raim() for parameter meanings.
 *
 * \return
 *   - `1`: repaired solution, using one fewer observation
 *          returns sid of removed measurement if removed_sid ptr is passed
 *
 *   - `-1`: no reasonable solution possible
 */
static s8 pvt_repair(double rx_state[],
                     const u8 n_used,
                     const navigation_measurement_t nav_meas[n_used],
                     double omp[n_used],
                     double H[4][4],
                     gnss_signal_t *removed_sid)
{
  /* Try solving with n-1 navigation measurements. */
  s8 one_less = n_used - 1;
  s8 bad_sat = -1;
  u8 num_passing = 0;

  const navigation_measurement_t *nav_meas_subset[n_used];
  for (s8 i = 0; i < n_used; i++) {
    nav_meas_subset[i] = &nav_meas[i];
  }

  /* Carefully ordered.
   * Permutes nav measurements so that each one is excluded from one test.
   */
  for (s8 drop = one_less; drop >= 0; drop--) {
    /* Swaps the last omitted value with the one at index `drop'.
     * On first iteration, does nothing (omits last nav_meas) */
    const navigation_measurement_t *temp;
    temp = nav_meas_subset[drop];
    nav_meas_subset[drop] = nav_meas_subset[one_less];
    nav_meas_subset[one_less] = temp;

    s8 flag = pvt_iter(rx_state, n_used - 1, nav_meas_subset, omp, H);

    if (flag == -1) {
      /* Didn't converge. */
      /* TODO(dsk) this may be unnecessary; use continue instead. */
      return -1;
    }

    if (residual_test(n_used-1, omp, rx_state, 0)) {
      num_passing++;
      bad_sat = drop;
    }
  }

  if (num_passing == 1) {
    /* Repair is possible by omitting bad_sat. Recalculate that solution. */
    for (s8 i = 0; i < n_used; i++) {
      nav_meas_subset[i] = &nav_meas[i];
    }
    nav_meas_subset[bad_sat] = nav_meas_subset[one_less];
    s8 flag = pvt_iter(rx_state, n_used - 1, nav_meas_subset, omp, H);
    assert(flag == 0);
    if (removed_sid) {
      *removed_sid = nav_meas[bad_sat].sid;
    }
    return 1;
  } else {
    return -1;
  }
}

/** Calculate pvt solution, perform RAIM check, attempt to repair if needed.
 *
 * See calc_PVT for parameter meanings.
 * \param rx_state reference to pvt solver state allocated in calc_PVT
 * \param n_used number of measurments
 * \param nav_meas array of measurements
 * \param disable_raim passing True will omit raim check/repair functionality
 * \param H see pvt_solve
 * \param removed_sid if not null and repair occurs, returns dropped sid
 * \param residual if not null, return double value of residual
 *
 * \return Non-negative values indicate success; see below
 *         For negative values, refer to pvt_err_msg().
 * Return values:
 *    `2`: solution ok, but raim check was not used
 *        (exactly 4 measurements, or explicitly disabled)
 *
 *    `1`: repaired solution, using one fewer observation
 *        returns sid of removed measurement if removed_sid ptr is passed
 *
 *    `0`: initial solution ok
 *
 *   - `-1`: repair failed
 *   - `-2`: not enough satellites to attempt repair
 *   - `-3`: pvt_iter didn't converge
 *
 *  Results stored in rx_state, H
 */
static s8 pvt_solve_raim(double rx_state[],
                         const u8 n_used,
                         const navigation_measurement_t nav_meas[n_used],
                         bool disable_raim,
                         double H[4][4],
                         gnss_signal_t *removed_sid,
                         double residual)
{
  double omp[n_used];

  assert(n_used <= MAX_CHANNELS);

  const navigation_measurement_t *nav_meas_ptrs[n_used];
  for (s8 i = 0; i < n_used; i++) {
    nav_meas_ptrs[i] = &nav_meas[i];
  }

  s8 flag = pvt_iter(rx_state, n_used, nav_meas_ptrs, omp, H);

  if (flag == -1) {
    /* Iteration didn't converge. Don't attempt to repair; too CPU intensive. */
    return -3;
  }
  if (flag >= 0 && (disable_raim || residual_test(n_used, omp, rx_state, &residual))) {
    /* Solution ok, or raim check disabled. */
    if (disable_raim || n_used == 4) {
      /* Residual test couldn't have detected an error. */
      return 2;
    }
    return 0;
  } else {
    if (n_used < 6) {
      /* Not enough measurements to repair.
       * 6 are needed because a 4 dimensional system is exactly constrained,
       * so the bad measurement can't be detected.
       */
      return -2;
    }
    return pvt_repair(rx_state, n_used, nav_meas, omp, H, removed_sid);
  }
}

/** Error strings for calc_PVT() negative (failure) return codes.
 *  e.g. `pvt_err_msg[-ret - 1]`
 *    where `ret` is the return value of calc_PVT(). */
const char *pvt_err_msg[] = {
  "PDOP too high",
  "Altitude unreasonable",
  "Velocity >= 1000 kts",
  "RAIM repair attempted, failed",
  "RAIM repair impossible (not enough measurements)",
  "Took too long to converge",
  "Not enough measurements for solution (< 4)",
};

/** Try to calculate a single point gps solution
 *
 * \param n_used number of measurments
 * \param nav_meas array of measurements
 * \param disable_raim passing True will omit raim check/repair functionality
 * \param soln output solution struct
 * \param dops output doppler information
 * \return Non-negative values indicate a valid solution.
 *   -  `2`: Solution converged but RAIM unavailable or disabled
 *   -  `1`: Solution converged, failed RAIM but was successfully repaired
 *   -  `0`: Solution converged and verified by RAIM
 *   - `-1`: PDOP is too high to yield a good solution.
 *   - `-2`: Altitude is unreasonable.
 *   - `-3`: Velocity is greater than or equal to 1000 kts.
 *   - `-4`: RAIM check failed and repair was unsuccessful
 *   - `-5`: RAIM check failed and repair was impossible (not enough measurements)
 *   - `-6`: pvt_iter didn't converge
 *   - `-7`: < 4 measurements
 */
s8 calc_PVT(const u8 n_used,
            const navigation_measurement_t nav_meas[n_used],
            bool disable_raim,
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

  if (n_used < 4) {
    return -7;
  }

  soln->valid = 0;
  soln->n_used = n_used; // Keep track of number of working channels

  gnss_signal_t removed_sid;
  s8 raim_flag = pvt_solve_raim(rx_state, n_used, nav_meas, disable_raim,
                                H, &removed_sid, 0);

  if (raim_flag < 0) {
    /* Didn't converge or least squares integrity check failed. */
    return raim_flag - 3;
  }

  /* Initial solution failed, but repair was successful. */
  if (raim_flag == 1) {
    soln->n_used--;
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
  normalize_gps_time(&soln->time);

  u8 ret;
  if ((ret = filter_solution(soln, dops))) {
    memset(soln, 0, sizeof(*soln));
    /* Reset position elements of state if solution fails. */
    rx_state[0] = 0;
    rx_state[1] = 0;
    rx_state[2] = 0;
    return -ret;
  }

  soln->valid = 1;

  return raim_flag;
}
