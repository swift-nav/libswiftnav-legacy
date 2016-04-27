# Copyright (C) 2015 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

"""Bit Utils

Bit field packing, unpacking and utility functions.

"""
from swiftnav.bits cimport parity as c_parity

def parity(x):
  '''
  Cython wrapper for parity function

  Parameters
  ----------
  x : int
    Value for parity computation

  Returns
  -------
  int
    Parity value: 1 or 0.
  '''
  return c_parity(x)

def getbitu(buff,  pos,  length):
  raise NotImplementedError

def getbits(buff,  pos,  length):
  raise NotImplementedError

def setbitu(buff,  pos,  length,  data):
  raise NotImplementedError

def setbits(buff,  pos, length, data):
  raise NotImplementedError
