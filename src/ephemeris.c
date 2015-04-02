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

#include "logging.h"
#include "linear_algebra.h"
#include "constants.h"
#include "ephemeris.h"

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
s8 calc_sat_state(const ephemeris_t *ephemeris, gps_time_t t,
                  double pos[3], double vel[3],
                  double *clock_err, double *clock_rate_err)
{
  assert(pos != NULL);
  assert(vel != NULL);
  assert(clock_err != NULL);
  assert(clock_rate_err != NULL);
  assert(ephemeris != NULL);

  /* Calculate satellite clock terms */

  /* Seconds from clock data reference time (toc) */
  double dt = gpsdifftime(t, ephemeris->toc);
  *clock_err = ephemeris->af0 + dt * (ephemeris->af1 + dt * ephemeris->af2)
               - ephemeris->tgd;
  *clock_rate_err = ephemeris->af1 + 2.0 * dt * ephemeris->af2;

  /* Seconds from the time from ephemeris reference epoch (toe) */
  dt = gpsdifftime(t, ephemeris->toe);

  /* If dt is greater than 4 hours our ephemeris isn't valid. */
  if (fabs(dt) > 4*3600) {
    log_warn("Using ephemeris older (or newer!) than 4 hours.\n");
    return -1;
  }

  /* Calculate position per IS-GPS-200D p 97 Table 20-IV */

  /* Semi-major axis in meters. */
  double a = ephemeris->sqrta * ephemeris->sqrta;
  /* Corrected mean motion in radians/sec. */
  double ma_dot = sqrt(GPS_GM / (a * a * a)) + ephemeris->dn;
  /* Corrected mean anomaly in radians. */
  double ma = ephemeris->m0 + ma_dot * dt;

  /* Iteratively solve for the Eccentric Anomaly
   * (from Keith Alter and David Johnston) */
  double ea = ma; /* Starting value for E. */
  double ea_old;
  double temp;
  double ecc = ephemeris->ecc;
  u32 count = 0;

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
  double einstein = GPS_F * ecc * ephemeris->sqrta * sin(ea);
  *clock_err += einstein;

  /* Begin calc for True Anomaly and Argument of Latitude */
  double temp2 = sqrt(1.0 - ecc * ecc);
  /* Argument of Latitude = True Anomaly + Argument of Perigee. */
  double al = atan2(temp2 * sin(ea), cos(ea) - ecc) + ephemeris->w;
  double al_dot = temp2 * ea_dot / temp;

  /* Calculate corrected argument of latitude based on position. */
  double cal = al + ephemeris->cus * sin(2.0 * al) + ephemeris->cuc * cos(2.0 * al);
  double cal_dot = al_dot * (1.0 + 2.0 * (ephemeris->cus * cos(2.0 * al)
                                          - ephemeris->cuc * sin(2.0 * al)));

  /* Calculate corrected radius based on argument of latitude. */
  double r = a * temp + ephemeris->crc * cos(2.0 * al)
             + ephemeris->crs * sin(2.0 * al);
  double r_dot = a * ecc * sin(ea) * ea_dot
                 + 2.0 * al_dot * (ephemeris->crs * cos(2.0 * al)
                                   - ephemeris->crc * sin(2.0 * al));

  /* Calculate inclination based on argument of latitude. */
  double inc = ephemeris->inc + ephemeris->inc_dot * dt
               + ephemeris->cic * cos(2.0 * al)
               + ephemeris->cis * sin(2.0 * al);
  double inc_dot = ephemeris->inc_dot
                   + 2.0 * al_dot * (ephemeris->cis * cos(2.0 * al)
                                     - ephemeris->cic * sin(2.0 * al));

  /* Calculate position and velocity in orbital plane. */
  double x = r * cos(cal);
  double y = r * sin(cal);
  double x_dot = r_dot * cos(cal) - y * cal_dot;
  double y_dot = r_dot * sin(cal) + x * cal_dot;

  /* Corrected longitude of ascenting node. */
  double om_dot = ephemeris->omegadot - GPS_OMEGAE_DOT;
  double om = ephemeris->omega0 + dt * om_dot
              - GPS_OMEGAE_DOT * ephemeris->toe.tow;

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

double predict_range(double rx_pos[3],
                     gps_time_t t,
                     ephemeris_t *ephemeris)
{
  double sat_pos[3];
  double sat_vel[3];
  double temp[3];
  double clock_err, clock_rate_err;

  calc_sat_state(ephemeris, t, sat_pos, sat_vel, &clock_err, &clock_rate_err);

  vector_subtract(3, sat_pos, rx_pos, temp); /* temp = sat_pos - rx_pos */
  return vector_norm(3, temp);
}

/** Is this ephemeris usable?
 *
 * \todo This should actually be more than just the "valid" flag.
 *       When we write an is_usable() function, lets use that instead
 *       of just es[prn].valid.
 *
 * \param eph Ephemeris struct
 * \return 1 if the ephemeris is valid and not too old.
 *         0 otherwise.
 */
u8 ephemeris_good(ephemeris_t *eph, gps_time_t t)
{
  /* Seconds from the time from ephemeris reference epoch (toe) */
  double dt = gpsdifftime(t, eph->toe);

  /* TODO: this doesn't exclude ephemerides older than a week so could be made
   * better. */
  return (eph->valid && eph->healthy && fabs(dt) < 4*3600);
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

  /* Subframe 1: SV health, T_GD, t_oc, a_f2, a_f1, a_f0 */

  /* GPS week number (mod 1024): Word 3, bits 20-30 */
  e->toe.wn = (frame_words[0][3-3] >> (30-10) & 0x3FF);
  e->toe.wn += GPS_WEEK_CYCLE * 1024;
  e->toc.wn = e->toe.wn;

  /* Health flag: Word 3, bit 17 */
  e->healthy = !(frame_words[0][3-3] >> (30-17) & 1);

  /* t_gd: Word 7, bits 17-24 */
  onebyte.u8 = frame_words[0][7-3] >> (30-24) & 0xFF;
  e->tgd = onebyte.s8 * pow(2,-31);

  /* t_oc: Word 8, bits 8-24 */
  e->toc.tow = (frame_words[0][8-3] >> (30-24) & 0xFFFF) * 16;

  /* a_f2: Word 9, bits 1-8 */
  onebyte.u8 = frame_words[0][9-3] >> (30-8) & 0xFF;
  e->af2 = onebyte.s8 * pow(2,-55);

  /* a_f1: Word 9, bits 9-24 */
  twobyte.u16 = frame_words[0][9-3] >> (30-24) & 0xFFFF;
  e->af1 = twobyte.s16 * pow(2,-43);

  /* a_f0: Word 10, bits 1-22 */
  fourbyte.u32 = frame_words[0][10-3] >> (30-22) & 0x3FFFFF;
  /* Shift to the left for sign extension */
  fourbyte.u32 <<= 10;
  /* Carry the sign bit back down and reduce to signed 22 bit value */
  fourbyte.s32 >>= 10;
  e->af0 = fourbyte.s32 * pow(2,-31);

  /* Subframe 2: crs, dn, m0, cuc, ecc, cus, sqrta, toe */

  /* crs: Word 3, bits 9-24 */
  twobyte.u16 = frame_words[1][3-3] >> (30-24) & 0xFFFF;
  e->crs = twobyte.s16 * pow(2,-5);

  /* dn: Word 4, bits 1-16 */
  twobyte.u16 = frame_words[1][4-3] >> (30-16) & 0xFFFF;
  e->dn = twobyte.s16 * pow(2,-43) * GPS_PI;

  /* m0: Word 4, bits 17-24 and word 5, bits 1-24 */
  fourbyte.u32 = ((frame_words[1][4-3] >> (30-24) & 0xFF) << 24)
               | (frame_words[1][5-3] >> (30-24) & 0xFFFFFF);
  e->m0 = fourbyte.s32 * pow(2,-31) * GPS_PI;

  /* cuc: Word 6, bits 1-16 */
  twobyte.u16 = frame_words[1][6-3] >> (30-16) & 0xFFFF;
  e->cuc = twobyte.s16 * pow(2,-29);

  /* ecc: Word 6, bits 17-24 and word 7, bits 1-24 */
  fourbyte.u32 = ((frame_words[1][6-3] >> (30-24) & 0xFF) << 24)
               | (frame_words[1][7-3] >> (30-24) & 0xFFFFFF);
  e->ecc = fourbyte.u32 * pow(2,-33);

  /* cus: Word 8, bits 1-16 */
  twobyte.u16 = frame_words[1][8-3] >> (30-16) & 0xFFFF;
  e->cus = twobyte.s16 * pow(2,-29);

  /* sqrta: Word 8, bits 17-24 and word 9, bits 1-24 */
  fourbyte.u32 = ((frame_words[1][8-3] >> (30-24) & 0xFF) << 24)
               | (frame_words[1][9-3] >> (30-24) & 0xFFFFFF);
  e->sqrta = fourbyte.u32 * pow(2,-19);

  /* t_oe: Word 10, bits 1-16 */
  e->toe.tow = (frame_words[1][10-3] >> (30-16) & 0xFFFF) * 16;


  /* Subframe 3: cic, omega0, cis, inc, crc, w, omegadot, inc_dot */

  /* cic: Word 3, bits 1-16 */
  twobyte.u16 = frame_words[2][3-3] >> (30-16) & 0xFFFF;
  e->cic = twobyte.s16 * pow(2,-29);

  /* omega0: Word 3, bits 17-24 and word 4, bits 1-24 */
  fourbyte.u32 = ((frame_words[2][3-3] >> (30-24) & 0xFF) << 24)
               | (frame_words[2][4-3] >> (30-24) & 0xFFFFFF);
  e->omega0 = fourbyte.s32 * pow(2,-31) * GPS_PI;

  /* cis: Word 5, bits 1-16 */
  twobyte.u16 = frame_words[2][5-3] >> (30-16) & 0xFFFF;
  e->cis = twobyte.s16 * pow(2,-29);

  /* inc (i0): Word 5, bits 17-24 and word 6, bits 1-24 */
  fourbyte.u32 = ((frame_words[2][5-3] >> (30-24) & 0xFF) << 24)
               | (frame_words[2][6-3] >> (30-24) & 0xFFFFFF);
  e->inc = fourbyte.s32 * pow(2,-31) * GPS_PI;

  /* crc: Word 7, bits 1-16 */
  twobyte.u16 = frame_words[2][7-3] >> (30-16) & 0xFFFF;
  e->crc = twobyte.s16 * pow(2,-5);

  /* w (omega): Word 7, bits 17-24 and word 8, bits 1-24 */
  fourbyte.u32 = ((frame_words[2][7-3] >> (30-24) & 0xFF) << 24)
               | (frame_words[2][8-3] >> (30-24) & 0xFFFFFF);
  e->w = fourbyte.s32 * pow(2,-31) * GPS_PI;

  /* Omega_dot: Word 9, bits 1-24 */
  fourbyte.u32 = frame_words[2][9-3] >> (30-24) & 0xFFFFFF;
  /* Shift left for sign extension */
  fourbyte.u32 <<= 8;
  /* Sign extend it */
  fourbyte.s32 >>= 8;
  e->omegadot = fourbyte.s32 * pow(2,-43) * GPS_PI;

  /* iode: Word 10, bits 1-8 */
  e->iode = frame_words[2][10-3] >> (30-8) & 0xFF;

  /* inc_dot (IDOT): Word 10, bits 9-22 */
  twobyte.u16 = frame_words[2][10-3] >> (30-22) & 0x3FFF;
  twobyte.u16 <<= 2;
  /* Sign extend */
  twobyte.s16 >>= 2;
  e->inc_dot = twobyte.s16 * pow(2,-43) * GPS_PI;

  e->valid = 1;
}

