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


#include "l2cbitstream.h"
#include <libswiftnav/edc.h>

#include <time.h>
#include <stdlib.h>
#include <string.h>

#define GPS_CNAV_MSG_PREAMBLE 0x8B
#define GPS_CNAV_MSG_OFFSET_IN_BITS 4
#define GPS_CNAV_MSG_CRC_LEN_BYTES 3
#define GPS_CNAV_MSG_TOW_COUNT_INC 8 // equals 12 seconds
#define GPS_CNAV_MSG_TOW_COUNT_MAX 403199
#define GPS_CNAV_MSG_TOW_COUNT_BITS 17
#define GPS_CNAV_MSG_TOW_COUNT_IGNORE_LSB 2

#define BITS_IN_BYTE 8


/* http://stackoverflow.com/a/23019348 */
static uint64_t extractBits(uint64_t t_In, uint8_t u_Len, uint8_t u_Idx)
{
    return t_In << (64 - u_Idx - u_Len) >> (64 - u_Len);
}


uint8_t getL2CMessageLength()
{
    return GPS_CNAV_MESSAGE_LEN_BYTES;
}


bool getL2CMessage(uint8_t *au_Message)
{
    uint8_t u_Idx = 0;
    uint32_t q_Crc = 0;
    struct timespec z_Spec;
    static uint32_t q_TowCount = 0;

    q_TowCount += GPS_CNAV_MSG_TOW_COUNT_INC;

    if (GPS_CNAV_MSG_TOW_COUNT_MAX < q_TowCount)
    {
        q_TowCount = 0;
    }

    clock_gettime(CLOCK_REALTIME, &z_Spec);
    srand(time(NULL) * z_Spec.tv_sec * z_Spec.tv_nsec);

    au_Message[u_Idx++] = GPS_CNAV_MSG_PREAMBLE >>\
        BITS_IN_BYTE - GPS_CNAV_MSG_OFFSET_IN_BITS;

    au_Message[u_Idx] = 0xF0 &\
        (GPS_CNAV_MSG_PREAMBLE << GPS_CNAV_MSG_OFFSET_IN_BITS);
    au_Message[u_Idx++] |= rand() & 0x0F;
    au_Message[u_Idx++] = rand();

    /* 17 MSBs of the TOW count -> ignore 2 LSBs */
    au_Message[u_Idx++] = extractBits(q_TowCount, BITS_IN_BYTE,
        GPS_CNAV_MSG_TOW_COUNT_BITS\
        - 1 * BITS_IN_BYTE + GPS_CNAV_MSG_TOW_COUNT_IGNORE_LSB);
    au_Message[u_Idx++] = extractBits(q_TowCount, BITS_IN_BYTE,
        GPS_CNAV_MSG_TOW_COUNT_BITS\
        - 2 * BITS_IN_BYTE + GPS_CNAV_MSG_TOW_COUNT_IGNORE_LSB);
    // last bit
    au_Message[u_Idx] = extractBits(q_TowCount, 1,
        GPS_CNAV_MSG_TOW_COUNT_IGNORE_LSB) << 7;

    au_Message[u_Idx++] |= (rand() & 0x7F);

    for (u_Idx; u_Idx < GPS_CNAV_MESSAGE_LEN_BYTES - GPS_CNAV_MSG_CRC_LEN_BYTES;
        u_Idx++)
    {
        au_Message[u_Idx] = rand();
    }

    q_Crc = crc24q(au_Message,
        GPS_CNAV_MESSAGE_LEN_BYTES - GPS_CNAV_MSG_CRC_LEN_BYTES, q_Crc);

    au_Message[u_Idx++] = extractBits(q_Crc, BITS_IN_BYTE, 2 * BITS_IN_BYTE);
    au_Message[u_Idx++] = extractBits(q_Crc, BITS_IN_BYTE, 1 * BITS_IN_BYTE);
    au_Message[u_Idx++] = extractBits(q_Crc, BITS_IN_BYTE, 0 * BITS_IN_BYTE);

    return true;
}


int32_t writeL2CToFile(int i_Fileno, uint64_t t_WantedMessageAmount)
{
    bool b_Offset = true;
    uint32_t q_Ret = 0;
    uint64_t t_MessageAmount = 0;
    uint8_t u_Idx = 0;
    uint8_t au_MessageToWrite[GPS_CNAV_MESSAGE_LEN_BYTES] = {0};
    uint8_t au_NextMessage[GPS_CNAV_MESSAGE_LEN_BYTES] = {0};
    FILE *pz_File = fdopen(i_Fileno, "wb");

    if (NULL == pz_File ||
        !getL2CMessage(au_NextMessage))
    {
        printf("ERROR\n");
        fflush(pz_File);
        return -1;
    }

    t_WantedMessageAmount++;
    while (t_MessageAmount < t_WantedMessageAmount)
    {
        uint8_t u_BytesReady = 0;
        while (u_BytesReady < GPS_CNAV_MESSAGE_LEN_BYTES)
        {
            if (b_Offset)
            {
                au_MessageToWrite[u_BytesReady] =\
                    au_NextMessage[u_Idx++] << GPS_CNAV_MSG_OFFSET_IN_BITS;
            }
            else
            {
                au_MessageToWrite[u_BytesReady++] = au_NextMessage[u_Idx++];
            }

            if (u_Idx >= GPS_CNAV_MESSAGE_LEN_BYTES)
            {
                if (!getL2CMessage(au_NextMessage))
                {
                    printf("ERROR\n");
                    fflush(pz_File);
                    return t_MessageAmount;
                }
                else
                {
                    u_Idx = 0;
                    t_MessageAmount++;
                    if (b_Offset)
                    {
                        au_MessageToWrite[u_BytesReady++] |=\
                            au_NextMessage[u_Idx++] & 0x0F;
                        b_Offset = false;
                    }
                    else
                    {
                        b_Offset = true;
                        continue;
                    }
                }
            }

            if (b_Offset)
            {
                au_MessageToWrite[u_BytesReady++] |=\
                    au_NextMessage[u_Idx] >> GPS_CNAV_MSG_OFFSET_IN_BITS;
            }
        }

        if (fwrite(au_MessageToWrite, sizeof(au_MessageToWrite[0]),
            sizeof(au_MessageToWrite), pz_File) != sizeof(au_MessageToWrite))
        {
            printf("ERROR writing file\n");
            fflush(pz_File);
            return t_MessageAmount;
        }
        memset(au_MessageToWrite, 0, sizeof(au_MessageToWrite)); 
    }

    fflush(pz_File);

    return t_MessageAmount;
}
