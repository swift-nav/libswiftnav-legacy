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

#ifndef LIBSWIFTNAV_L2C_CAPABILITY_H
#define LIBSWIFTNAV_L2C_CAPABILITY_H

#include <libswiftnav/common.h>
#include <libswiftnav/nav_msg.h>

void decode_l2c_capability(const u32 *subframe4_words, u32 *l2c_cpbl);

#endif /* LIBSWIFTNAV_L2C_CAPABILITY_H */
