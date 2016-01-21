/*
* Copyright (C) 2016 Swift Navigation Inc.
* Contact: Pasi Miettinen <pasi.miettinen@exafore.com>
*
* This source is subject to the license found in the file 'LICENSE' which must
* be be distributed together with this source. All other rights reserved.
*
* THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
* EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
*/

#include <libswiftnav/common.h>

#include <stdbool.h>
#include <stdint.h>

#define CNAVMSG_LEN_BYTES 38

u8 get_l2c_message_length(void);
bool get_l2c_message(u8 *au_message);
s32 write_l2c_to_file(int i_fileno, u64 t_msg_amount);
