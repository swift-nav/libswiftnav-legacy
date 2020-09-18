'''
Copyright (C) 2016 Swift Navigation Inc.
Contact: Dmitry Tatarinov <dmitry.tatarinov@exafore.com>

This source is subject to the license found in the file 'LICENSE' which must
be be distributed together with this source. All other rights reserved.

THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

'''

'''
The Python script performs convolutional encoding of a file for Viterbi
decoder testing purposes.

USAGE: 
  python v27_enc.py <input_file_name> <output_file_name>

OUTPUT:
  Symbol file <output_file_name> which has encoded data, each byte is
  represented as 0xff for '1' and 0x00 for '0' parity bits.
'''

import struct
import os
import sys

G1 = 0171 #generator polinomial for p1
G2 = 0133 #generator polinomial for p2

def parity(st):
  return (0x6996 >> ((st ^ (st >> 4)) & 15)) & 1

state = 0

if len(sys.argv) < 3:
  print 'Not enough input data. Usage: v27_enc.py <input_file> <output_file>'
  exit()
else:
  f_in_name = str(sys.argv[1])
  f_out_name = str(sys.argv[2])

f_out_hard = open(f_out_name,'wb')
f_in = open(f_in_name,'rb')

f_in_size = os.path.getsize(f_in_name)

byte = 0

pre =  [0,0,0,0]
post = [0,0,0,0]

# 2 bits delay to have data be byte aligned
j = 3
while j >= 0:
  f_out_hard.write(chr(0))
  j -= 1

#encode prefix first
for b in pre:
  j = 7
  while j >= 0:
    state = (((b >> j) & 1) << 6) | (state >> 1)
    p1 = parity(state & G1)
    p2 = parity(state & G2)
    #here is output SN-like for Viterbi decoder testing 
    if p1 == 1:
      f_out_hard.write(chr(0xff))
    else:
      f_out_hard.write(chr(0))

    if p2 == 1:
      f_out_hard.write(chr(0xff))
    else:
      f_out_hard.write(chr(0))

    j -= 1

#now encode data from file
with open(f_in_name,'rb') as f_in:
  while f_in_size:
    byte = struct.unpack('B',f_in.read(1))[0] #read byte from input file

    j = 7
    while j >= 0: #for each bit in the byte
      state = (((byte >> j) & 1) << 6) | (state >> 1)
      p1 = parity(state & G1) #calculate parity
      p2 = parity(state & G2)

      #here is output SN-like for Viterbi decoder testing 
      if p1 == 1:
        f_out_hard.write(chr(0xff))
      else:
        f_out_hard.write(chr(0))

      if p2 == 1:
        f_out_hard.write(chr(0xff))
      else:
        f_out_hard.write(chr(0))

      j -= 1
        
    f_in_size -= 1;

#encode postfix
for b in post:
  j = 7
  while j >= 0:
    state = (((b >> j) & 1) << 6) | (state >> 1)
    p1 = parity(state & G1)
    p2 = parity(state & G2)
    #here is output SN-like for Viterbi decoder testing 
    if p1 == 1:
      f_out_hard.write(chr(0xff))
    else:
      f_out_hard.write(chr(0))

    if p2 == 1:
      f_out_hard.write(chr(0xff))
    else:
      f_out_hard.write(chr(0))

    j -= 1

f_out_hard.close()
f_in.close()
