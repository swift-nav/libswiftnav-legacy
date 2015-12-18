/*
 * Copyright (C) 2010 Swift Navigation Inc.
 * Contact: Henry Hallam <henry@swift-nav.com>
 *          Fergus Noble <fergus@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "logging.h"
#include "linear_algebra.h"
#include "constants.h"
#include "ephemeris.h"

/** \defgroup ephemeris Ephemeris
 * Functions and calculations related to the GPS ephemeris.
 * \{ */

/** Calculate satellite position, velocity from xyz ephemeris.
 *
 * References:
 *   -# WAAS Specification FAA-E-2892b 4.4.11
 *
 * \param e Ephemeris struct
 * \param pos Array into which to write calculated satellite position [m]
 * \param vel Array into which to write calculated satellite velocity [m/s]
 *
 * \return  0 on success,
 *         -1 if ephemeris is older (or newer) than 4 hours
 */
static s8 calc_sat_state_xyz(const ephemeris_t *e, gps_time_t t,
                             double pos[3], double vel[3],
                             double *clock_err, double *clock_rate_err)
{
  /*
   *Ephemeris did not get time-stammped when it was received.
   */
  if (e->toe.wn == 0)
    return -1;

  if (!e->valid || !e->healthy)
    return -1;

  const ephemeris_xyz_t *ex = &e->xyz;
  u8 ndays = t.tow / DAY_SECS;
  double tod = t.tow - DAY_SECS * ndays;
  double dt = tod - ex->toa;

  if (dt > DAY_SECS/2)
    dt -= DAY_SECS;
  else if (dt < -DAY_SECS/2)
    dt += DAY_SECS;

  vel[0] = ex->rate[0];
  vel[1] = ex->rate[1];
  vel[2] = ex->rate[2];

  pos[0] = ex->pos[0] + ex->rate[0] * dt +
           0.5 * ex->acc[0] * pow(dt, 2);
  pos[1] = ex->pos[1] + ex->rate[1] * dt +
           0.5 * ex->acc[1] * pow(dt, 2);
  pos[2] = ex->pos[2] + ex->rate[2] * dt +
           0.5 * ex->acc[2] * pow(dt, 2);

  memcpy(clock_err, &(ex->a_gf0), sizeof(double));
  memcpy(clock_rate_err, &(ex->a_gf1), sizeof(double));

  return 0;
}

/** Calculate satellite position, velocity and clock offset from ephemeris.
 *
 * References:
 *   -# IS-GPS-200D, Section 20.3.3.3.3.1 and Table 20-IV
 *
 * \param t GPS time at which to calculate the satellite state
 * \param ephemeris Ephemeris struct
 * \param pos Array into which to write calculated satellite position [m]
 * \param vel Array into which to write calculated satellite velocity [m/s]
 * \param clock_err Pointer to where to store the calculated satellite clock
 *                  error [s]
 * \param clock_rate_err Pointer to where to store the calculated satellite
 *                       clock error [s/s]
 *
 * \return  0 on success,
 *         -1 if ephemeris is older (or newer) than 4 hours
 */
static s8 calc_sat_state_kepler(const ephemeris_t *e,
                                gps_time_t t,
                                double pos[3], double vel[3],
                                double *clock_err, double *clock_rate_err)
{
  const ephemeris_kepler_t *k = &e->kepler;

  /* Calculate satellite clock terms */

  /* Seconds from clock data reference time (toc) */
  double dt = gpsdifftime(t, k->toc);
  *clock_err = k->af0 + dt * (k->af1 + dt * k->af2) - k->tgd;
  *clock_rate_err = k->af1 + 2.0 * dt * k->af2;

  /* Seconds from the time from ephemeris reference epoch (toe) */
  dt = gpsdifftime(t, e->toe);

  /* If dt is greater than fit_interval hours our ephemeris isn't valid. */
  if (fabs(dt) > ((u32)e->fit_interval)*60*60) {
    log_error("Using ephemeris outside validity period, dt = %+.0f", dt);
    return -1;
  }

  /* Calculate position per IS-GPS-200D p 97 Table 20-IV */

  /* Semi-major axis in meters. */
  double a = k->sqrta * k->sqrta;
  /* Corrected mean motion in radians/sec. */
  double ma_dot = sqrt(GPS_GM / (a * a * a)) + k->dn;
  /* Corrected mean anomaly in radians. */
  double ma = k->m0 + ma_dot * dt;

  /* Iteratively solve for the Eccentric Anomaly
   * (from Keith Alter and David Johnston) */
  double ea = ma; /* Starting value for E. */
  double ea_old;
  double temp;
  double ecc = k->ecc;
  u8 count = 0;

  /* TODO: Implement convergence test using integer difference of doubles,
   * http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm */
  /* TODO: Bound number of iterations. */
  do {
    ea_old = ea;
    temp = 1.0 - ecc * cos(ea_old);
    ea = ea + (ma - ea_old + ecc * sin(ea_old)) / temp;
    count++;
    if (count > 5)
      break;
  } while (fabs(ea - ea_old) > 1.0E-14);

  double ea_dot = ma_dot / temp;

  /* Relativistic correction term. */
  double einstein = GPS_F * ecc * k->sqrta * sin(ea);
  *clock_err += einstein;

  /* Begin calc for True Anomaly and Argument of Latitude */
  double temp2 = sqrt(1.0 - ecc * ecc);
  /* Argument of Latitude = True Anomaly + Argument of Perigee. */
  double al = atan2(temp2 * sin(ea), cos(ea) - ecc) + k->w;
  double al_dot = temp2 * ea_dot / temp;

  /* Calculate corrected argument of latitude based on position. */
  double cal = al + k->cus * sin(2.0 * al) + k->cuc * cos(2.0 * al);
  double cal_dot = al_dot * (1.0 + 2.0 * (k->cus * cos(2.0 * al)
                                          - k->cuc * sin(2.0 * al)));

  /* Calculate corrected radius based on argument of latitude. */
  double r = a * temp + k->crc * cos(2.0 * al) + k->crs * sin(2.0 * al);
  double r_dot = a * ecc * sin(ea) * ea_dot
                 + 2.0 * al_dot * (k->crs * cos(2.0 * al)
                                   - k->crc * sin(2.0 * al));

  /* Calculate inclination based on argument of latitude. */
  double inc = k->inc + k->inc_dot * dt + k->cic * cos(2.0 * al)
               + k->cis * sin(2.0 * al);
  double inc_dot = k->inc_dot
                   + 2.0 * al_dot * (k->cis * cos(2.0 * al)
                                     - k->cic * sin(2.0 * al));

  /* Calculate position and velocity in orbital plane. */
  double x = r * cos(cal);
  double y = r * sin(cal);
  double x_dot = r_dot * cos(cal) - y * cal_dot;
  double y_dot = r_dot * sin(cal) + x * cal_dot;

  /* Corrected longitude of ascenting node. */
  double om_dot = k->omegadot - GPS_OMEGAE_DOT;
  double om = k->omega0 + dt * om_dot - GPS_OMEGAE_DOT * e->toe.tow;

  /* Compute the satellite's position in Earth-Centered Earth-Fixed
   * coordiates. */
  pos[0] = x * cos(om) - y * cos(inc) * sin(om);
  pos[1] = x * sin(om) + y * cos(inc) * cos(om);
  pos[2] = y * sin(inc);

  /* Compute the satellite's velocity in Earth-Centered Earth-Fixed
   * coordiates. */
  temp = y_dot * cos(inc) - y * sin(inc) * inc_dot;
  vel[0] = -om_dot * pos[1] + x_dot * cos(om) - temp * sin(om);
  vel[1] = om_dot * pos[0] + x_dot * sin(om) + temp * cos(om);
  vel[2] = y * cos(inc) * inc_dot + y_dot * sin(inc);

  return 0;
}

/** Calculate satellite position, velocity and clock offset from ephemeris.
 *
 * Dispatch to internal function for Kepler/XYZ ephemeris depending on
 * constellation.
 *
 * \param t GPS time at which to calculate the satellite state
 * \param e Ephemeris struct
 * \param pos Array into which to write calculated satellite position [m]
 * \param vel Array into which to write calculated satellite velocity [m/s]
 * \param clock_err Pointer to where to store the calculated satellite clock
 *                  error [s]
 * \param clock_rate_err Pointer to where to store the calculated satellite
 *                       clock error [s/s]
 *
 * \return  0 on success,
 *         -1 if ephemeris is older (or newer) than 4 hours
 */
s8 calc_sat_state(const ephemeris_t *e, gps_time_t t,
                  double pos[3], double vel[3],
                  double *clock_err, double *clock_rate_err)
{
  assert(pos != NULL);
  assert(vel != NULL);
  assert(clock_err != NULL);
  assert(clock_rate_err != NULL);
  assert(e != NULL);

  switch (e->sid.constellation) {
  case CONSTELLATION_GPS:
    return calc_sat_state_kepler(e, t, pos, vel, clock_err, clock_rate_err);
  case CONSTELLATION_SBAS:
    return calc_sat_state_xyz(e, t, pos, vel, clock_err, clock_rate_err);
  default:
    assert("unsupported constellation");
    return -1;
  }
}

/** Is this ephemeris usable?
 *
 * \todo This should actually be more than just the "valid" flag.
 *       When we write an is_usable() function, lets use that instead
 *       of just es[prn].valid.
 *
 * \param eph Ephemeris struct
 * \param t
 * \return 1 if the ephemeris is valid and not too old.
 *         0 otherwise.
 */
u8 ephemeris_good(const ephemeris_t *eph, gps_time_t t)
{
  /* Seconds from the time from ephemeris reference epoch (toe) */
  double dt = gpsdifftime(t, eph->toe);

  /* TODO: this doesn't exclude ephemerides older than a week so could be made
   * better. */
  return (eph->valid && fabs(dt) < ((u32)eph->fit_interval)*60*60);
}

/** Convert a GPS URA index into a value.
*
* |param index URA index.
* |return the URA in meters.
*/
float decode_ura_index(const u8 index) {
  static float values[16] = {
    [0]  = 2.0f,
    [1]  = 2.8f,
    [2]  = 4.0f,
    [3]  = 5.7f,
    [4]  = 8.0f,
    [5]  = 11.3f,
    [6]  = 16.0f,
    [7]  = 32.0f,
    [8]  = 64.0f,
    [9]  = 128.0f,
    [10] = 256.0f,
    [11] = 512.0f,
    [12] = 1024.0f,
    [13] = 2048.0f,
    [14] = 4096.0f,
    [15] = 6144.0f,
  };

  return values[index];
}

/** Calculate the GPS ephemeris curve fit interval.
*
* |param fit_interval_flag The curve fit interval flag. 0 is 4 hours, 1 is >4 hours.
* |param iodc The IODC value.
* |return the curve fit interval in hours.
*/
u8 decode_fit_interval(u8 fit_interval_flag, u16 iodc) {
  u8 fit_interval = 4;

  if (fit_interval_flag) {
    fit_interval = 6;

    if ((iodc >= 240) && (iodc <= 247)) {
      fit_interval = 8;
    } else if (((iodc >= 248) && (iodc <= 255)) || (iodc == 496)) {
      fit_interval = 14;
    } else if ((iodc >= 497) && (iodc <= 503)) {
      fit_interval = 26;
    } else if ((iodc >= 504) && (iodc <= 510)) {
      fit_interval = 50;
    } else if ((iodc == 511) || ((iodc >= 752) && (iodc <= 756))) {
      fit_interval = 74;
    } else if (iodc == 757) {
      fit_interval = 98;
    }
  }

  return fit_interval;
}

/** Decode ephemeris from L1 C/A GPS navigation message frames.
 *
 * \note This function does not check for parity errors. You should check the
 *       subframes for parity errors before calling this function.
 *
 * References:
 *   -# IS-GPS-200D, Section 20.3.2 and Figure 20-1
 *
 * \param frame_words Array containing words 3 through 10 of subframes
 *                    1, 2 and 3. Word is in the 30 LSBs of the u32.
 * \param e Pointer to an ephemeris struct to fill in.
 */
void decode_ephemeris(u32 frame_words[3][8], ephemeris_t *e)
{
  assert(frame_words != NULL);
  assert(e != NULL);
  assert(e->sid.constellation == CONSTELLATION_GPS);
  ephemeris_kepler_t *k = &e->kepler;

  char buf[SID_STR_LEN_MAX];
  sid_to_string(buf, sizeof(buf), e->sid);

  /* These unions facilitate signed/unsigned conversion and sign extension. */
  union {
    s8 s8;
    u8 u8;
  } onebyte;

  union
  {
    s16 s16;
    u16 u16;
  } twobyte;

  union
  {
    s32 s32;
    u32 u32;
  } fourbyte;


  /* Subframe 1: WN, URA, SV health, T_GD, IODC, t_oc, a_f2, a_f1, a_f0 */

  /* GPS week number (mod 1024): Word 3, bits 1-10 */
  e->toe.wn = frame_words[0][3-3] >> (30-10) & 0x3FF;
  e->toe.wn += GPS_WEEK_CYCLE * 1024;
  k->toc.wn = e->toe.wn;

  /* URA: Word 3, bits 13-16 */
  /* Value of 15 is unhealthy */
  u8 ura_index = frame_words[0][3-3] >> (30-16) & 0xF;
  e->ura = decode_ura_index(ura_index);
  log_debug("URA for %s = index %d, value %.1f", buf, ura_index, e->ura);

  /* NAV data and signal health bits: Word 3, bits 17-22 */
  // TODO: Check only L1 C/A health
  u8 health_bits = frame_words[0][3-3] >> (30-22) & 0x3F;
  log_debug("Health bits for %s = 0x%02x", buf, health_bits);
  e->healthy = (health_bits == 0x00) && (ura_index < 15);
  if (!e->healthy) {
    log_warn("Latest ephemeris for %s is unhealthy. Ignoring satellite.", buf);
  }

  /* t_gd: Word 7, bits 17-24 */
  onebyte.u8 = frame_words[0][7-3] >> (30-24) & 0xFF;
  k->tgd = onebyte.s8 * pow(2,-31);

  /* iodc: Word 3, bits 23-24 and word 8, bits 1-8 */
  k->iodc = ((frame_words[0][3-3] >> (30-24) & 0x3) << 8)
          | (frame_words[0][8-3] >> (30-8) & 0xFF);

  /* t_oc: Word 8, bits 8-24 */
  k->toc.tow = (frame_words[0][8-3] >> (30-24) & 0xFFFF) * 16;

  /* a_f2: Word 9, bits 1-8 */
  onebyte.u8 = frame_words[0][9-3] >> (30-8) & 0xFF;
  k->af2 = onebyte.s8 * pow(2,-55);

  /* a_f1: Word 9, bits 9-24 */
  twobyte.u16 = frame_words[0][9-3] >> (30-24) & 0xFFFF;
  k->af1 = twobyte.s16 * pow(2,-43);

  /* a_f0: Word 10, bits 1-22 */
  fourbyte.u32 = frame_words[0][10-3] >> (30-22) & 0x3FFFFF;
  /* Shift to the left for sign extension */
  fourbyte.u32 <<= 10;
  /* Carry the sign bit back down and reduce to signed 22 bit value */
  fourbyte.s32 >>= 10;
  k->af0 = fourbyte.s32 * pow(2,-31);

  /* Subframe 2: IODE, crs, dn, m0, cuc, ecc, cus, sqrta, toe, fit_interval */

  /* iode: Word 3, bits 1-8 */
  u8 iode_sf2 = frame_words[1][3-3] >> (30-8) & 0xFF;

  /* crs: Word 3, bits 9-24 */
  twobyte.u16 = frame_words[1][3-3] >> (30-24) & 0xFFFF;
  k->crs = twobyte.s16 * pow(2,-5);

  /* dn: Word 4, bits 1-16 */
  twobyte.u16 = frame_words[1][4-3] >> (30-16) & 0xFFFF;
  k->dn = twobyte.s16 * pow(2,-43) * GPS_PI;

  /* m0: Word 4, bits 17-24 and word 5, bits 1-24 */
  fourbyte.u32 = ((frame_words[1][4-3] >> (30-24) & 0xFF) << 24)
               | (frame_words[1][5-3] >> (30-24) & 0xFFFFFF);
  k->m0 = fourbyte.s32 * pow(2,-31) * GPS_PI;

  /* cuc: Word 6, bits 1-16 */
  twobyte.u16 = frame_words[1][6-3] >> (30-16) & 0xFFFF;
  k->cuc = twobyte.s16 * pow(2,-29);

  /* ecc: Word 6, bits 17-24 and word 7, bits 1-24 */
  fourbyte.u32 = ((frame_words[1][6-3] >> (30-24) & 0xFF) << 24)
               | (frame_words[1][7-3] >> (30-24) & 0xFFFFFF);
  k->ecc = fourbyte.u32 * pow(2,-33);

  /* cus: Word 8, bits 1-16 */
  twobyte.u16 = frame_words[1][8-3] >> (30-16) & 0xFFFF;
  k->cus = twobyte.s16 * pow(2,-29);

  /* sqrta: Word 8, bits 17-24 and word 9, bits 1-24 */
  fourbyte.u32 = ((frame_words[1][8-3] >> (30-24) & 0xFF) << 24)
               | (frame_words[1][9-3] >> (30-24) & 0xFFFFFF);
  k->sqrta = fourbyte.u32 * pow(2,-19);

  /* t_oe: Word 10, bits 1-16 */
  e->toe.tow = (frame_words[1][10-3] >> (30-16) & 0xFFFF) * 16;

  /* fit_interval_flag: Word 10, bit 17 */
  u8 fit_interval_flag = frame_words[1][10-3] >> (30-17) & 0x1;
  e->fit_interval = decode_fit_interval(fit_interval_flag, k->iodc);
  log_debug("Fit interval for %s = %d", buf, e->fit_interval);

  /* Subframe 3: cic, omega0, cis, inc, crc, w, omegadot, IODE, inc_dot */

  /* cic: Word 3, bits 1-16 */
  twobyte.u16 = frame_words[2][3-3] >> (30-16) & 0xFFFF;
  k->cic = twobyte.s16 * pow(2,-29);

  /* omega0: Word 3, bits 17-24 and word 4, bits 1-24 */
  fourbyte.u32 = ((frame_words[2][3-3] >> (30-24) & 0xFF) << 24)
               | (frame_words[2][4-3] >> (30-24) & 0xFFFFFF);
  k->omega0 = fourbyte.s32 * pow(2,-31) * GPS_PI;

  /* cis: Word 5, bits 1-16 */
  twobyte.u16 = frame_words[2][5-3] >> (30-16) & 0xFFFF;
  k->cis = twobyte.s16 * pow(2,-29);

  /* inc (i0): Word 5, bits 17-24 and word 6, bits 1-24 */
  fourbyte.u32 = ((frame_words[2][5-3] >> (30-24) & 0xFF) << 24)
               | (frame_words[2][6-3] >> (30-24) & 0xFFFFFF);
  k->inc = fourbyte.s32 * pow(2,-31) * GPS_PI;

  /* crc: Word 7, bits 1-16 */
  twobyte.u16 = frame_words[2][7-3] >> (30-16) & 0xFFFF;
  k->crc = twobyte.s16 * pow(2,-5);

  /* w (omega): Word 7, bits 17-24 and word 8, bits 1-24 */
  fourbyte.u32 = ((frame_words[2][7-3] >> (30-24) & 0xFF) << 24)
               | (frame_words[2][8-3] >> (30-24) & 0xFFFFFF);
  k->w = fourbyte.s32 * pow(2,-31) * GPS_PI;

  /* Omega_dot: Word 9, bits 1-24 */
  fourbyte.u32 = frame_words[2][9-3] >> (30-24) & 0xFFFFFF;
  /* Shift left for sign extension */
  fourbyte.u32 <<= 8;
  /* Sign extend it */
  fourbyte.s32 >>= 8;
  k->omegadot = fourbyte.s32 * pow(2,-43) * GPS_PI;

  /* iode: Word 10, bits 1-8 */
  k->iode = frame_words[2][10-3] >> (30-8) & 0xFF;

  /* inc_dot (IDOT): Word 10, bits 9-22 */
  twobyte.u16 = frame_words[2][10-3] >> (30-22) & 0x3FFF;
  twobyte.u16 <<= 2;
  /* Sign extend */
  twobyte.s16 >>= 2;
  k->inc_dot = twobyte.s16 * pow(2,-43) * GPS_PI;

  /* Both IODEs and IODC (8 LSBs) must match */
  log_debug("Check ephemeris with %s. IODC = 0x%03x IODE = 0x%02x and 0x%02x.", buf, k->iodc, iode_sf2, k->iode);
  e->valid = (iode_sf2 == k->iode) && (k->iode == (k->iodc & 0xFF));
  if (!e->valid) {
    log_warn("Latest ephemeris for %s had IODC/IODE mismatch. Ignoring ephemeris.", buf);
  }
}

static bool ephemeris_xyz_equal(const ephemeris_xyz_t *a,
                                const ephemeris_xyz_t *b)
{
  return (a->iod == b->iod) &&
         (a->toa == b->toa) &&
         (a->a_gf0 == b->a_gf0) &&
         (a->a_gf1 == b->a_gf1) &&
         (memcmp(a->pos, b->pos, sizeof(a->pos)) == 0) &&
         (memcmp(a->rate, b->rate, sizeof(a->rate)) == 0) &&
         (memcmp(a->acc, b->acc, sizeof(a->acc)) == 0);
}

static bool ephemeris_kepler_equal(const ephemeris_kepler_t *a,
                                   const ephemeris_kepler_t *b)
{
  return (a->iodc == b->iodc) &&
         (a->iode == b->iode) &&
         (a->tgd == b->tgd) &&
         (a->crs == b->crs) &&
         (a->crc == b->crc) &&
         (a->cuc == b->cuc) &&
         (a->cus == b->cus) &&
         (a->cic == b->cic) &&
         (a->cis == b->cis) &&
         (a->dn == b->dn) &&
         (a->m0 == b->m0) &&
         (a->ecc == b->ecc) &&
         (a->sqrta == b->sqrta) &&
         (a->omega0 == b->omega0) &&
         (a->omegadot == b->omegadot) &&
         (a->w == b->w) &&
         (a->inc == b->inc) &&
         (a->inc_dot == b->inc_dot) &&
         (a->af0 == b->af0) &&
         (a->af1 == b->af1) &&
         (a->af2 == b->af2) &&
         (a->toc.wn == b->toc.wn) &&
         (a->toc.tow == b->toc.tow);
}

bool ephemeris_equal(const ephemeris_t *a, const ephemeris_t *b)
{
  if (!sid_is_equal(a->sid, b->sid) ||
      (a->ura != b->ura) ||
      (a->fit_interval != b->fit_interval) ||
      (a->valid != b->valid) ||
      (a->healthy != b->healthy) ||
      (a->toe.wn != b->toe.wn) ||
      (a->toe.tow != b->toe.tow))
    return false;

  switch (a->sid.constellation) {
  case CONSTELLATION_GPS:
    return ephemeris_kepler_equal(&a->kepler, &b->kepler);
  case CONSTELLATION_SBAS:
    return ephemeris_xyz_equal(&a->xyz, &b->xyz);
  default:
    assert("unsupported constellation");
    return false;
  }
}

/** \} */
