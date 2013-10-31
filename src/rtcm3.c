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

#include "bits.h"
#include "edc.h"
#include "rtcm3.h"

#define RTCM3_PREAMBLE 0xD3

/** \defgroup rtcm3 RTCM v3
 * RTCM v3.1 Format message encoding and decoding.
 *
 * Implements RTCM Standard 10403.1
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

/** \} */

