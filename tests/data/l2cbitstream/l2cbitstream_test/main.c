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


#include <l2cbitstream.h>
#include <libswiftnav/edc.h>

#include <inttypes.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>

#define SLEEP_US 1000000

static volatile s32 i_keep_running = 1;

void int_handler(int i_dummy)
{
  i_keep_running = 0;
}

int main(int argc, char **argv)
{
  u8 au_message[CNAVMSG_LEN_BYTES] = {0};
  u8 u_idx = 0;
  printf("l2cbitstream unit tests starting\n");

  FILE *pz_file = fopen("test.bin", "wb");
  write_l2c_to_file(fileno(pz_file), 10);
  fclose(pz_file);

  signal(SIGINT, int_handler);
  while (i_keep_running) {
    if (get_l2c_message(au_message)) {
      for (u_idx = 0; u_idx < CNAVMSG_LEN_BYTES; u_idx++) {
        printf("%02" PRIx8, au_message[u_idx]);
      }

      printf("\n\n");
      usleep(SLEEP_US);
    }
  }

  printf("\nl2cbitstream unit tests done\n");
  return 0;
}
