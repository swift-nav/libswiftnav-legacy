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

/*****************************************************************************
 * Automatically generated from sbp.yaml with generate.py, do not hand edit! *
 *****************************************************************************/

#ifndef LIBSWIFTNAV_SBP_MESSAGES_H
#define LIBSWIFTNAV_SBP_MESSAGES_H

#include "common.h"

((* for m in msgs *))
/** (((m.short_desc)))
(((m.desc|commentify)))
 */
#define SBP_(((m.name))) ((('0x%04X'|format(m.id))))

typedef struct __attribute__((packed)) {
((*- for f in m.fields *))
  (((f.type.ljust(m.max_type_len)))) ((((f.name+';').ljust(m.max_name_len+1)))) /**< (((f.desc))) ((* if f.units *))[(((f.units)))] ((* endif *))*/
((*- endfor *))
} sbp_(((m.name|lower)))_t;

((* endfor *))
#endif /* LIBSWIFTNAV_SBP_MESSAGES_H */


