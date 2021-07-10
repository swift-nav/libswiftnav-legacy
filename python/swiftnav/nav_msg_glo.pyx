# Copyright (C) 2016 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

# cython: embedsignature=True

from ephemeris cimport ephemeris_t
from fmt_utils import fmt_repr
from time cimport gps_time_t
from signal cimport gnss_signal_t
from swiftnav.ephemeris import Ephemeris


cdef class NavMsgGlo:

  def __cinit__(self):
    nav_msg_init_glo(&self._thisptr)

  def __getattr__(self, k):
    return self._thisptr.get(k)

  def __repr__(self):
    return fmt_repr(self)

  def to_dict(self):
    return self._thisptr

  def from_dict(self, d):
    self._thisptr = d

  def __reduce__(self):
    return (rebuild_NavMsgGlo, tuple(self.to_dict().items()))

  def updateEphemeris(self, Ephemeris e):
    return process_string_glo(&self._thisptr, &e._thisptr)

  def update(self, bit_val):
    return nav_msg_update_glo(&self._thisptr, bit_val)

  def getTow(self):
    return nav_msg_get_tow_glo(&self._thisptr)

  def detectError(self):
    return error_detection_glo(&self._thisptr)

  def isDecodeDone(self):
    return self._thisptr.decode_done != 0

  def __richcmp__(self, other, op):
    """
    Weird Cython comparison method. See
    http://docs.cython.org/src/userguide/special_methods.html.
    """
    if op == 2:
      return self._equal(other)
    elif op == 3:
      return not self._equal(other)
    else:
      raise NotImplementedError

  def _equal(self, other):
    """
    Compare equality between self and another :class:`NavMsg` object.

    Parameters
    ----------
    other : :class:`NavMsgGlo` object
      The :class:`NavMsgGlo` to test equality against.

    Return
    ------
    out : bool
      True if the passed :class:`NavMsg` object is identical.

    """
    if self.to_dict().keys() != other.to_dict().keys():
      return False

    for k in self.to_dict().keys():
      if self.to_dict()[k] != other.to_dict()[k]:
        return False

    return True


def rebuild_NavMsgGlo(*reduced):
  """
  Rebuild NavMsg for unpickling.

  Parameters
  ----------
  reduced: tuple
    Tuple of dict of NavMsg nav_msg_t struct fields

  Returns
  -------
  out: :class:`NavMsg` instance
    Rebuilt :class:`NavMsg` instance
  """
  nm = NavMsgGlo()
  nm.from_dict(dict(reduced))
  return nm

