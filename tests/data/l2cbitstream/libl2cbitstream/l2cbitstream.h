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


#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

#define GPS_CNAV_MESSAGE_LEN_BYTES 38


uint8_t getL2CMessageLength();


/*------------------------------------------------- getL2CMessage -----
 *  Function getL2CMessage
 *
 *  Purpose:  Fills buffer with L2C message frame including preamble
 *            and CRC-24q. Bits are right-aligned so four leftmost bits
 *            of the 38 byte buffer are always zero.
 *
 *  Parameters:
 *      au_Message (OUT) -- buffer containing the L2C message.
 * 
 *  Dependencies: crc24q() from libswiftnav
 *
 *  Returns:  True if returned message is valid. Otherwise returns
 *            False.
 *-------------------------------------------------------------------*/
bool getL2CMessage(uint8_t *au_Message);


/*------------------------------------------------- writeL2CToFile -----
 *  Function writeL2CToFile
 *
 *  Purpose:  Writes L2C messages into given file. Utilizes
 *            getL2CMessage function. Flushes, but doesn't close
 *            the file.
 *
 *  Parameters:
 *      i_Fileno (IN) -- File number (given by fileno()) which is
 *                        used to save the data.
 *      t_WantedMessageAmount (IN) -- Indicates how many messages
 *                                    should be written. Depending
 *                                    on the offset the final result
 *                                    can be 1 message off.
 *
 *  Returns:  How many messages was requested from getL2CMessage.
 *-------------------------------------------------------------------*/
int32_t writeL2CToFile(int i_Fileno, uint64_t t_WantedMessageAmount);