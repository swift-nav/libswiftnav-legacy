/*
 * Copyright (C) 2016 Swift Navigation Inc.
 * Contact: Dmitry Tatarinov <dmitry.tatarinov@exafore.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */
#include <stdio.h>
#include <assert.h>

#include <libswiftnav/common.h>
#include <libswiftnav/nav_msg.h>
#include <libswiftnav/l2c_capability.h>

/** \defgroup l2c_capability L2C Capability
 * Functions used in l2c capability
 * \{ */

/** Decode L2C capability data from Nav message.
 * \param subframe4_words Array containing words 3 through 10
 *                        of subframe 4.
 *                        SV configuration data is in words 3 through 8
 * \param l2c_cpbl output Pointer to 32-bit L2C capability container
 *
 */
void decode_l2c_capability(const u32 *subframe4_words, u32 *l2c_cpbl)
{
  assert(subframe4_words != NULL);
  assert(l2c_cpbl != NULL);

  *l2c_cpbl = 0;

  struct sv_conf_location {
    u8 n_sv; /* for how many SVs config data in the Word */
    u8 end_bit_n; /* position of MSB for the 1st SV in the Word */
  };

  const struct sv_conf_location sv_conf_loc[6] = {
    {4,12},{6,4},{6,4},{6,4},{6,4},{4,4},
  };

  u8 sv_id = 0;

  /* go through words 3-8 */
  for (u8 i = 3; i <= 8; i++ ) {

    /* go through all sv inside the word */
    for (u8 j = 0; j < sv_conf_loc[i-3].n_sv; j++ ) {
      /* get SV config bits take into account only 3 LSB */
      u8 sv_conf =
      subframe4_words[i-3] >> (30 - (sv_conf_loc[i-3].end_bit_n + j*4)) & 7;

      /* set or clear appropriate capability bit,
       * refer pg. 115-116 of IS-200H for the criteria
       * uses an open upper bound to ensure we track L2C on future satellite
       * generations launched after ICD was updated */
      if (sv_conf >= 2)
        *l2c_cpbl |= (1 << sv_id);

      sv_id++;
    }
  }
}

/** \} */
