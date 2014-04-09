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

#include "bits.h"
#include "edc.h"
#include "rtcm3.h"

#define RTCM3_PREAMBLE 0xD3 /**< RTCM v3 Frame sync / preamble byte. */
#define PRUNIT_GPS 299792.458 /**< RTCM v3 Unit of GPS Pseudorange (m) */

#define CLIGHT  299792458.0         /* speed of light (m/s) */
#define FREQ1   1.57542e9           /* L1/E1  frequency (Hz) */
#define LAMBDA1 (CLIGHT / FREQ1)

/** \addtogroup io Input / Output
 * \{ */

/** \defgroup rtcm3 RTCM v3
 * RTCM v3.1 Format message encoding and decoding.
 *
 * Implements RTCM Standard 10403.1
 *
 * DFxxx codes indicate a corresponding Data Field as described in RTCM 10403.1
 * Table 3.4-1.
 * \{ */

/** Check RTCM frame header and CRC valid.
 *
 * \param buff A pointer to the RTCM message buffer to check.
 * \return If valid then return length of the data message in range 0-1023.
 *         Returns a negative number if the frame is invalid:
 *          - `-1` : Preamble not valid
 *          - `-2` : CRC mismatch
 */
s16 rtcm3_check_frame(u8 *buff)
{
  /* Check preamble byte. */
  u8 preamble = getbitu(buff, 0, 8);
  if (preamble != RTCM3_PREAMBLE)
    return -1;

  s16 len = getbitu(buff, 14, 10);

  /* Calculate CRC of frame header and data message. */
  u32 crc_calc = crc24q(buff, len + 3, 0);

  /* Check CRC in frame. */
  u32 crc_frame = getbitu(buff, (len+3)*8, 24);

  if (crc_calc != crc_frame)
    return -2;

  return len;
}

/** Write RTCM frame header and CRC into a buffer.
 *
 * The buffer should already contain the data message starting at the
 * 4th byte of the buffer, i.e.
 *
 *     data_message = &buff[3]
 *
 * The buffer must have 3 bytes free at the start to contain the header and
 * must leave 3 bytes past the end free to contain the CRC. The total length of
 * the buffer should be `len+6` bytes.
 *
 * `len` should be in the range 0-1023 as per the RTCM v3 standard.
 *
 * \param len The length of the data message contained in the buffer.
 * \param buff A pointer to the RTCM message buffer.
 * \return Zero on success, -1 if `len` is too large.
 */
s8 rtcm3_write_frame(u16 len, u8 *buff)
{
  if (len > 1023)
    return -1;

  /* Set preamble, reserved and length bits. */
  setbitu(buff, 0, 8, RTCM3_PREAMBLE);
  setbitu(buff, 8, 6, 0);
  setbitu(buff, 14, 10, len);

  /* Calculate CRC of frame header and data message. */
  u32 crc = crc24q(buff, len + 3, 0);

  /* Write CRC to end of frame. */
  setbitu(buff, (len+3)*8, 24, crc);

  return 0;
}

/** Write RTCM header for observation message types 1001..1004.
 *
 * The data message header will be written starting from byte zero of the
 * buffer. If the buffer also contains a frame header then be sure to pass a
 * pointer to the start of the data message rather than a pointer to the start
 * of the frame buffer. The RTCM observation header is 8 bytes (64 bits) long.
 *
 * If the Synchronous GNSS Message Flag is set to `0`, it means that no further
 * GNSS observables referenced to the same Epoch Time will be transmitted. This
 * enables the receiver to begin processing the data immediately after decoding
 * the message. If it is set to `1`, it means that the next message will
 * contain observables of another GNSS source referenced to the same Epoch
 * Time.
 *
 * Divergence-free Smoothing Indicator values:
 *
 * Indicator | Meaning
 * --------- | ----------------------------------
 *     0     | Divergence-free smoothing not used
 *     1     | Divergence-free smoothing used
 *
 * GPS Smoothing Interval indicator values are listed in RTCM 10403.1 Table
 * 3.4-4, reproduced here:
 *
 * Indicator | Smoothing Interval
 * --------- | ------------------
 *  000 (0)  |   No smoothing
 *  001 (1)  |   < 30 s
 *  010 (2)  |   30-60 s
 *  011 (3)  |   1-2 min
 *  100 (4)  |   2-4 min
 *  101 (5)  |   4-8 min
 *  110 (6)  |   >8 min
 *  111 (7)  |   Unlimited
 *
 * \param buff A pointer to the RTCM data message buffer.
 * \param type Message type number, i.e. 1001..1004 (DF002).
 * \param id Reference station ID (DF003).
 * \param t GPS time of epoch (DF004).
 * \param sync Synchronous GNSS Flag (DF005).
 * \param n_sat Number of GPS satellites included in the message (DF006).
 * \param div_free GPS Divergence-free Smoothing Indicator (DF007).
 * \param smooth GPS Smoothing Interval indicator (DF008).
 */
void rtcm3_write_header(u8 *buff, u16 type, u16 id, gps_time_t t,
                        u8 sync, u8 n_sat, u8 div_free, u8 smooth)
{
  setbitu(buff, 0, 12, type);
  setbitu(buff, 12, 12, id);
  setbitu(buff, 24, 30, round(t.tow*1e3));
  setbitu(buff, 54, 1, sync);
  setbitu(buff, 55, 5, n_sat);
  setbitu(buff, 60, 1, div_free);
  setbitu(buff, 61, 3, smooth);
}

/** Read RTCM header for observation message types 1001..1004.
 *
 * The data message header will be read starting from byte zero of the
 * buffer. If the buffer also contains a frame header then be sure to pass a
 * pointer to the start of the data message rather than a pointer to the start
 * of the frame buffer. The RTCM observation header is 8 bytes (64 bits) long.
 *
 * All return values are written into the parameters passed by reference.
 *
 * \param buff A pointer to the RTCM data message buffer.
 * \param type Message type number, i.e. 1001..1004 (DF002).
 * \param id Reference station ID (DF003).
 * \param tow GPS time of week of the epoch (DF004).
 * \param sync Synchronous GNSS Flag (DF005).
 * \param n_sat Number of GPS satellites included in the message (DF006).
 * \param div_free GPS Divergence-free Smoothing Indicator (DF007).
 * \param smooth GPS Smoothing Interval indicator (DF008).
 */
void rtcm3_read_header(u8 *buff, u16 *type, u16 *id, double *tow,
                       u8 *sync, u8 *n_sat, u8 *div_free, u8 *smooth)
{
  *type = getbitu(buff, 0, 12);
  *id = getbitu(buff, 12, 12);
  *tow = getbitu(buff, 24, 30) / 1e3;
  *sync = getbitu(buff, 54, 1);
  *n_sat = getbitu(buff, 55, 5);
  *div_free = getbitu(buff, 60, 1);
  *smooth = getbitu(buff, 61, 3);
}

/** Convert a lock time in seconds into a RTCMv3 Lock Time Indicator value.
 * See RTCM 10403.1, Table 3.4-2.
 *
 * \param time Lock time in seconds.
 * \return Lock Time Indicator value.
 */
static u8 to_lock_ind(u32 time)
{
  if (time < 24)
    return time;
  if (time < 72)
    return (time + 24) / 2;
  if (time < 168)
    return (time + 120) / 4;
  if (time < 360)
    return (time + 408) / 8;
  if (time < 744)
    return (time + 1176) / 16;
  if (time < 937)
    return (time + 3096) / 32;
  return 127;
}

/** Convert a RTCMv3 Lock Time Indicator value into a minimum lock time in seconds.
 * See RTCM 10403.1, Table 3.4-2.
 *
 * \param lock Lock Time Indicator value.
 * \return Minimum lock time in seconds.
 */
static u32 from_lock_ind(u8 lock)
{
  if (lock < 24)
    return lock;
  if (lock < 48)
    return 2*lock - 24;
  if (lock < 72)
    return 4*lock - 120;
  if (lock < 96)
    return 8*lock - 408;
  if (lock < 120)
    return 16*lock - 1176;
  if (lock < 127)
    return 32*lock - 3096;
  return 937;
}

/** Generate RTCMv3 formatted GPS observation fields.
 * Currently L1 only.
 *
 * \param nm Struct containing the observation.
 * \param amb The GPS Integer L1 Pseudorange Modulus Ambiguity (DF014).
 * \param pr The GPS L1 Pseudorange (DF011).
 * \param ppr The GPS L1 PhaseRange – L1 Pseudorange (DF012).
 * \param lock The GPS L1 Lock Time Indicator (DF013).
 * \param cnr The GPS L1 CNR (DF015).
 */
static void gen_obs_gps(navigation_measurement_t *nm,
                        u8 *amb, u32 *pr, s32 *ppr, u8 *lock, u8 *cnr)
{
  /* Calculate GPS Integer L1 Pseudorange Modulus Ambiguity (DF014). */
  *amb = (u8)(nm->raw_pseudorange / PRUNIT_GPS);

  /* Calculate GPS L1 Pseudorange (DF011). */
  *pr = (u32)lround((nm->raw_pseudorange - *amb * PRUNIT_GPS) / 0.02);

  /* Construct pseudorange value as transmitted. */
  double prc = *pr * 0.02 + *amb * PRUNIT_GPS;

  /* L1 phaserange - L1 pseudorange */
  double cp_pr = nm->carrier_phase - prc / (CLIGHT / FREQ1);

  /* If the phaserange and pseudorange have diverged close to the limits of the
   * data field (20 bits) then we modify the carrier phase by an integer amount
   * to bring it back into range an reset the phase lock time to zero to reset
   * the integer ambiguity.
   * The spec suggests adjusting by 1500 cycles but I calculate the range to be
   * +/- 1379 cycles. Limit to just 1000 as that should still be plenty. */
  if (fabs(cp_pr) > 1000) {
    nm->lock_time = 0;
    nm->carrier_phase -= (s32)cp_pr;
    cp_pr -= (s32)cp_pr;
  }

  /* Calculate GPS L1 PhaseRange – L1 Pseudorange (DF012). */
  *ppr = lround(cp_pr * (CLIGHT / FREQ1) / 0.0005);

  if (lock)
    *lock = to_lock_ind(nm->lock_time);

  if (cnr)
    /* TODO: When we start actually using cnr values then we can directly use
     * the value here. */
    *cnr = (u8)((10.0*log10(nm->snr) + 40.0) * 4.0);
}

/** Encode an RTCMv3 message type 1002 (Extended L1-Only GPS RTK Observables)
 * Message type 1002 has length `64 + n_sat*74` bits. Returned message length
 * is rounded up to the nearest whole byte.
 *
 * \param buff A pointer to the RTCM data message buffer.
 * \param id Reference station ID (DF003).
 * \param t GPS time of epoch (DF004).
 * \param n_sat Number of GPS satellites included in the message (DF006).
 * \param nm Struct containing the observation.
 * \param sync Synchronous GNSS Flag (DF005).
 * \return The message length in bytes.
 */
u16 rtcm3_encode_1002(u8 *buff, u16 id, gps_time_t t, u8 n_sat,
                      navigation_measurement_t *nm, u8 sync)
{
  rtcm3_write_header(buff, 1002, id, t, sync, n_sat, 0, 0);

  u16 bit = 64; /* Start at end of header. */

  u32 pr;
  s32 ppr;
  u8 amb, lock, cnr;

  for (u8 i=0; i<n_sat; i++) {
    gen_obs_gps(&nm[i], &amb, &pr, &ppr, &lock, &cnr);

    setbitu(buff, bit, 6,  nm[i].prn + 1); bit += 6;
    /* TODO: set GPS code indicator if we ever support P(Y) code measurements. */
    setbitu(buff, bit, 1,  0);    bit += 1;
    setbitu(buff, bit, 24, pr);   bit += 24;
    setbits(buff, bit, 20, ppr);  bit += 20;
    setbitu(buff, bit, 7,  lock); bit += 7;
    setbitu(buff, bit, 8,  amb);  bit += 8;
    setbitu(buff, bit, 8,  cnr);  bit += 8;
  }

  /* Round number of bits up to nearest whole byte. */
  return (bit + 7) / 8;
}

/** Decode an RTCMv3 message type 1002 (Extended L1-Only GPS RTK Observables)
 *
 * \param buff A pointer to the RTCM data message buffer.
 * \param id Reference station ID (DF003).
 * \param tow GPS time of week of epoch (DF004).
 * \param n_sat Number of GPS satellites included in the message (DF006).
 * \param nm Struct containing the observation.
 * \param sync Synchronous GNSS Flag (DF005).
 * \return If valid then return 0.
 *         Returns a negative number if the message is invalid:
 *          - `-1` : Message type mismatch
 *          - `-2` : Message uses unsupported P(Y) code
 */
s8 rtcm3_decode_1002(u8 *buff, u16 *id, double *tow, u8 *n_sat,
                     navigation_measurement_t *nm, u8 *sync)
{
  u16 type;
  u8 div_free, smooth;

  rtcm3_read_header(buff, &type, id, tow, sync, n_sat, &div_free, &smooth);

  if (type != 1002)
    /* Unexpected message type. */
    return -1;

  /* TODO: Fill in t->wn. */

  if (!nm)
    /* No nav meas pointer, probably just interested in
     * n_sat so we are all done. */
    return 0;

  u16 bit = 64;
  for (u8 i=0; i<*n_sat; i++) {
    /* TODO: Handle SBAS prns properly, numbered differently in RTCM? */
    nm[i].prn = getbitu(buff, bit, 6) - 1; bit += 6;

    u8 code = getbitu(buff, bit, 1); bit += 1;
    /* TODO: When we start storing the signal/system etc. properly we can
     * store the code flag in the nav meas struct. */
    if (code == 1)
      /* P(Y) code not currently supported. */
      return -2;

    u32 pr = getbitu(buff, bit,24); bit += 24;
    s32 ppr = getbits(buff, bit,20); bit += 20;
    u8 lock = getbitu(buff, bit, 7); bit += 7;
    u8 amb = getbitu(buff, bit, 8); bit += 8;
    u8 cnr = getbitu(buff, bit, 8); bit += 8;

    nm[i].raw_pseudorange = 0.02*pr + PRUNIT_GPS*amb;
    nm[i].carrier_phase = (nm[i].raw_pseudorange + 0.0005*ppr) / (CLIGHT / FREQ1);
    nm[i].lock_time = from_lock_ind(lock);
    nm[i].snr = pow(10.0, ((cnr / 4.0) - 40.0) / 10.0);
  }

  return 0;
}


/** \} */
/** \} */


