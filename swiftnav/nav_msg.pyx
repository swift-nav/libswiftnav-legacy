# Copyright (C) 2012 Swift Navigation Inc.
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

cdef class NavMsg:
  def __cinit__(self):
    nav_msg_c.nav_msg_init(&self.state)
    self.eph_valid = False

  def update(self, corr_prompt_real):
    tow = nav_msg_c.nav_msg_update(&self.state, corr_prompt_real)
    if nav_msg_c.subframe_ready(&self.state):
      #self.eph.valid = 0
      nav_msg_c.process_subframe(&self.state, &self.eph)
      if self.eph.valid:
        self.eph_valid = True
      else:
        self.eph_valid = False
        #return Ephemeris(&self.eph)
    return tow if tow != 0 else None

  def gps_week_num(self):
    if self.eph_valid:
      return self.eph.wn
    else:
      return None

  #property cptr:
    #def __get__(self):
      #return &self.eph


