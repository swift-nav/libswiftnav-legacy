/*
 * Copyright (C) 2014-2015 Swift Navigation Inc.
 * Contact: Ian Horn <ian@swift-nav.com>
 *          Fergus Noble <fergus@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */


#ifndef LIBSWIFTNAV_FILTER_UTILS_H
#define LIBSWIFTNAV_FILTER_UTILS_H

#include <common.h>

s8 assign_de_mtx(u8 num_sats, const sdiff_t *sats_with_ref_first,
                 const double ref_ecef[3], double *DE);

#endif /* LIBSWIFTNAV_FILTER_UTILS_H */

