# Copyright (C) 2016 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

from fmt_utils import fmt_repr
from libc.string cimport memcpy, memset
cimport numpy as np
import numpy as np

cdef class CNavMsgDecoder:

  def __init__(self, **kwargs):
    memset(&self._thisptr, 0, sizeof(cnav_msg_decoder_t))
    cnav_msg_decoder_init(&self._thisptr)

  def __repr__(self):
    return fmt_repr(self)

  def decode(self, u8 symbol, CNavMsg msg):
    cdef u32 delay = 0
    res = cnav_msg_decoder_add_symbol(&self._thisptr,
                                      symbol,
                                      &msg._thisptr,
                                      &delay)
    return res, delay

cdef class CNavMsg:
  def __init__(self, **kwargs):
    memset(&self._thisptr, 0, sizeof(cnav_msg_t))
    if 'prn' in kwargs:
      self._thisptr.prn = kwargs.pop('prn')
    if 'msg_id' in kwargs:
      self._thisptr.msg_id = kwargs.pop('msg_id')
    if 'tow' in kwargs:
      self._thisptr.tow = kwargs.pop('tow')
    if 'alert' in kwargs:
      self._thisptr.alert = kwargs.pop('alert')

  def __getattr__(self, k):
    return self._thisptr.get(k)

  def __repr__(self):
    return fmt_repr(self)

  def to_dict(self):
    return self._thisptr

  def from_dict(self, d):
    self._thisptr = d

  def __reduce__(self):
    return (rebuild_CNavMsg, tuple([tuple(self.to_dict().items())]))
  
  def getPrn(self):
    return self._thisptr.prn

  def getTow(self):
    return self._thisptr.tow

  def getMsgId(self):
    return self._thisptr.msg_id

  def getAlert(self):
    return self._thisptr.alert

def rebuild_CNavMsg(reduced):
  """
  Rebuild CNavMsg for unpickling.

  Parameters
  ----------
  reduced: tuple
    Tuple of dict of NavMsg cnav_msg_t struct fields

  Returns
  -------
  out: :class:`CNavMsg` instance
    Rebuilt :class:`CNavMsg` instance
  """
  nm = CNavMsg()
  nm.from_dict(dict(reduced))
  return nm


cdef class CNavRawMsg:
  @staticmethod
  def generate(prn, msg_id, tow):
    cdef np.ndarray[np.uint8_t, ndim=1, mode="c"] tmp_ = np.ndarray(38, dtype=np.uint8)
    res = get_l2c_message(&tmp_[0], prn, msg_id, tow)
    return np.unpackbits(tmp_)[4:]
