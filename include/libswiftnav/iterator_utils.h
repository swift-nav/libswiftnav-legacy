/*
 * Copyright (C) 2014 Swift Navigation Inc.
 * Contact: Ian Horn <ian@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */


#ifndef LIBSWIFTNAV_ITERATOR_UTILS_H
#define LIBSWIFTNAV_ITERATOR_UTILS_H
 
#include "single_diff.h"
#include "common.h"

void map_sdiff_prn(const void *arg, const void *it_current, void *state_current);

#endif /* LIBSWIFTNAV_ITERATOR_UTILS */
