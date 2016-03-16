# Copyright (C) 2012 Swift Navigation Inc.
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

cdef class NavMsg:

  def __cinit__(self):
    nav_msg_init(&self._thisptr)

  def __getattr__(self, k):
    return self._thisptr.get(k)

  def __repr__(self):
    return fmt_repr(self)

  def to_dict(self):
    return self._thisptr

  def from_dict(self, d):
    self._thisptr = d

  def __reduce__(self):
    return (rebuild_NavMsg, tuple(self.to_dict().items()))

  def update(self, bit_val):
    return nav_msg_update(&self._thisptr, bit_val)

  def subframe_ready(self):
    return subframe_ready(&self._thisptr)

  def process_subframe(self, e):
    cdef ephemeris_t tmp = e._thisptr
    return process_subframe(&self._thisptr, &tmp)

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
    other : :class:`NavMsg` object
      The :class:`NavMsg` to test equality against.

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


def rebuild_NavMsg(*reduced):
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
  nm = NavMsg()
  nm.from_dict(dict(reduced))
  return nm

