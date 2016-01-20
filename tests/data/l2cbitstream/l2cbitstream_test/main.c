/*
* Copyright (C) 2016 Swift Navigation Inc.
* Contact: Pasi Miettinen <pmiettinen@exafore.com>
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

#include <stdint.h>
#include <stdio.h>
#include <unistd.h>

#define SLEEP_US 1000000

int main(int argc, char **argv)
{
    uint8_t au_Message[GPS_CNAV_MESSAGE_LEN_BYTES] = {0};
    uint8_t u_Idx = 0;
    printf("l2cbitstream unit tests starting\n");

    FILE *pz_File = fopen("test.bin", "wb");
    writeL2CToFile(fileno(pz_File), 10);
    fclose(pz_File);

    while (true)
    {
        if (getL2CMessage(au_Message))
        {
            for (u_Idx = 0; u_Idx < GPS_CNAV_MESSAGE_LEN_BYTES; u_Idx++)
            {
                printf("%02x", au_Message[u_Idx]);
            }
            
            printf("\n\n");
            usleep(SLEEP_US);
        }
    }

    printf("l2cbitstream unit tests done\n");
    return 0;
}
