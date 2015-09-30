/*
 * Copyright (C) 2015 Swift Navigation Inc.
 * Contact: Vlad Ungureanu <vvu@vdev.ro>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */


#include <stdlib.h>
#include <string.h>

#include "signal.h"

signal_t sbas_index_to_sid(u8 index)
{
  signal_t sid;

  memset(&sid, 0, sizeof(signal_t));

  sid.band = L1_BAND;
  sid.constellation = SBAS_CONSTELLATION;

  switch (index) {
    case 0:
      sid.prn = 132;
      break;
    case 1:
      sid.prn = 134;
      break;
    case 2:
      sid.prn = 137;
      break;
  }

  return sid;
}

u8 sbas_sid_to_index(signal_t sid)
{
  switch(sid.prn) {
    case 132:
      return 0;
    case 134:
      return 1;
    case 137:
      return 2;
  }

  //TODO FIX THIS ASAP
  return 100;
}

