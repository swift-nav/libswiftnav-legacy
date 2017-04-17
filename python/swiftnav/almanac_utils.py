#!/usr/bin/env python
# Copyright (C) 2011-2014 Swift Navigation Inc.
# Contact: Fergus Noble <fergus@swift-nav.cogps>
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.

import math as m
from .time import datetime2gpst, gpst_components2datetime
from .almanac import Almanac
import datetime
import urllib2
import os

WPR = (-2712219.0, -4316338.0, 3820996.0)
SECONDS_IN_WEEK = (7*24*60*60)
NUMBER_ROLLOVERS =  1

MAX_DOWNLOAD_ATTEMPTS = 50
# Not sure the best way to specify max TOA diff
# but I think getting within a few days is good enough
MAX_TOA_DIFF = SECONDS_IN_WEEK/2.0

"""
Tools for downloading and operating on GPS almanac information.
"""

class AlmanacError(Exception):
  """Base class for exceptions in this module."""
  pass

class UsageError(AlmanacError):
  pass

class ConnectivityError(AlmanacError):
  pass


class AlmanacData(object):
  """
  Class to represent GPS alamanac information
  """

  def __init__(self, week=None, tow=None, cachedir=None, cacheext='alm'):
    """
    Construct an almanac.  If week or TOW are not supplied,
    they default to the current system time
    """
    self.sats={}
    self.text=""
    if week:
      self.week_number = week
      if tow:
        self.tow = tow
      else:
        self.tow = 0
    # if no weeknumber is specified, we set weeknumber and tow to current time
    else:
      gpstime = datetime2gpst(datetime.datetime.today())
      self.week_number = gpstime.wn%1024
      self.tow = gpstime.tow
      print ("No weeknumber specified.  Defaulting to current time.")
      # warn user if something was done illogically
      if tow:
        print ("Ignoring time of week without week number passed to "
                 "Almanc constructor.  Defaulting to current time.")
    self._get_almanac(cachedir, cacheext)

  def prn_vis_dop(self, prn, tow=None, pos=WPR, el_mask=0.0):
    """
    Get the doppler and elevation for one satellite only in the almanac.
    If TOW is not specified the tow from almanac construction is used
    Returns:
    --------
    A Tuple of floats representing
    doppler (hz)
    elevation (radians)
    """
    if self.sats:
      if not tow:
        tow = self.tow
      if self.sats.get(prn, None):
        doppler = self.sats[prn].calc_doppler(tow, pos, week=self.week_number)
        azimuth, elevation = self.sats[prn].calc_az_el(tow, pos, week=self.week_number)
        if not ((elevation > el_mask*(m.pi/180)) and self.sats[prn].healthy):
          doppler, elevation = None, None
        return doppler, elevation
      else:
        return (None, None)
    else:
      raise UsageError("No satellite dictionary initialized.")

  def _alm_week_number(self):
    """
    Return the week number of the almanac stored in our sats dictionary.
    Uses smallest PRN
    """
    if self.sats:
      return self.sats[min(self.sats)].gps['week']
    else:
      raise UsageError("No satellite dictionary initialized.")


  def _alm_toa(self):
    """
    Return the toa stored in our sats dictionary.
    Uses smallest PRN.
    """
    if self.sats:
      return self.sats[min(self.sats)].gps['toa']
    else:
      raise UsageError("No satellite dictionary initialized.")


  def get_distance(self):
    """
    Return the difference between the almanac's toa and that
    requested by the user in float number of seconds.
    """
    if self.sats:
      return (self._alm_toa() - self.tow
      + (self._alm_week_number()  - self.week_number) * SECONDS_IN_WEEK)
    else:
      raise UsageError("No satellite dictionary initialized.")


  def valid(self, maxd=MAX_TOA_DIFF):
    """
    Checks if the almanac is valid by comparing the TOA difference to a max value
    to the first satellite's week number.  We always want the TOA to be in the future
    where the difference would be postive
    """

    try:
      current_difference = self.get_distance()
      return current_difference > 0 and current_difference < maxd
    except UsageError:
      return false

  def clear(self):
    """
    remove sats dictionary and text from object
    """
    self.sats.clear()
    self.text = ""

  def _get_almanac(self, cachedir=None, ext='alm'):
    """
    get the almanac specified by the constructor
    """
    # First check the cache
    # cached files are named according wo WN.TOA.extension
    if ext==None:
      ext=''
    if cachedir:
      print ("Attempting to retreive cached almanac"
             " from directory {0} with extension \'.{1}\'.").format(cachedir, ext)
      files = [f for f in os.listdir(cachedir) if (f.endswith(ext) and
                                                  f.startswith(str(self.week_number)))]
      for f in files:
        self.load_file(os.path.join(cachedir, f))
        if self.valid():
          print "Valid almanac found with filename {0}".format(f)
          return True
      print "No cached Almanac found."
     # if the cache failed, download and maybe cache the results each time
    print "Downloading Almanac."
    if self._download(cachedir, ext):
       print ("Downloaded alamanac with week number {0} "
              "and toa {1}.").format(self._alm_week_number(), self._alm_toa())
       return True
    return False


  def _download(self, cachedir, ext):
    """
    Get the alamanac from the web
    Cache anything downloaded if cachedir is specified
    --------
    True upon successful download
    False upon failure
    """
    year = gpst_components2datetime(self.week_number + 1024 * NUMBER_ROLLOVERS,
                                    self.tow).year
    # almanac is update by navcen at least 2 times per week, so it is reasonable
    # to start at the year's week number * 2 when searching for the almanac
    i = (self.week_number % 52) * 2
    counter = 0
    try:
      handle = urllib2.urlopen('http://www.navcen.uscg.gov/?pageName=currentAlmanac&format=yuma')
      self.process_handle(handle)
      while not self.valid():
        if counter > MAX_DOWNLOAD_ATTEMPTS:
          print ("Error downloading almanac for year {0} and week {1}:"
               "Reached max interations.").format(year, self.week_number)
          raise(ConnectivityError("Reached maximum iterations in attempt to download Almanac."))
        url = 'http://www.navcen.uscg.gov/?Do=getAlmanac&almanac=' + str(i) + '&year=' +str(year)
        handle = urllib2.urlopen(url, timeout=60)
        self.process_handle(handle)
        # if specified cache our results
        if cachedir:
          self.save_file(os.path.join(cachedir, ".".join([str(self._alm_week_number()),
                                                          str(int(self._alm_toa())),
                                                          ext])))
        # increment our url iteration number
        # if we are more than one week before, we can increment by 2 * the week difference without going too far
        if -self.get_distance() > SECONDS_IN_WEEK:
          i += 2 *(self.week_number - self._alm_week_number())
        # as we get closer, we should only increment by 1 issue of almanac
        else:
          i += 1
        # always increment our iteration counter to exit loop
        counter += 1
    except urllib2.URLError as e:
        raise(ConnectivityError(e.reason))
    # we will only hit here with a a valid almanac
    return True


  def load_file(self, filename):
    """
    Load the almanac from a file
    """
    with open(filename, 'r') as fp:
      try:
        self.process_handle(fp)
      except IOError:
        errorstr = "IOError in load_almanac_file for file %s".format(filename)
        raise(ConnectivityError(errorstr))

  def save_file(self, filename):
    """
    Save the alamanac to a file
    """
    with open(filename, 'w+') as fp:
      if self.text:
        try:
          fp.write(self.text)
        except IOError:
          errorstr = "IOError in save_almanac_file for file %s".format(filename)
          raise(ConnectivityError(errorstr))

  def process_handle(self, handle):
    """
    Process a handle to the yuma block
    into a prn keyed dict of satellite objects
    Store the text as self.text for possible file saving
    """
    self.clear()
    yuma = handle.readlines()

    if yuma:
      blocks = []
      self.text = ''.join(yuma)
      for (n, line) in enumerate(yuma):
        if line[:3] == "ID:":
          bl = yuma[n:n+13]
          fields = map(lambda x: x[25:], bl)
          gnss_signal = dict(
            sat = int(fields[0]),
            band = 0,
            constellation = 0
          )
          gps_alm = dict(
            ecc     = float(fields[2]),
            toa     = float(fields[3]),
            inc     = float(fields[4]),
            rora    = float(fields[5]),
            a       = float(fields[6])**2,
            raaw    = float(fields[7]),
            argp    = float(fields[8]),
            ma      = float(fields[9]),
            af0     = float(fields[10]),
            af1     = float(fields[11]),
            week    = int(fields[12]),
          )
          alm = dict(
            sid     = gnss_signal,
            healthy = (int(fields[1]) == 0),
            valid   = True,
            gps     = gps_alm,
          )
          satAlmanac = Almanac(**alm)
          self.sats[satAlmanac.sid['sat']] = satAlmanac

  def get_dopps(self, tow=None, location=WPR, el_mask=0.0):
    """
    Give a location and TOW, return a list of the sats and their orbit info
    """
    if not tow:
      tow = self.tow
    if self.sats:
      dopps = map(lambda k: (k,)+self.prn_vis_dop(k, tow, location, el_mask),
                                                  self.sats.keys())
      dopps = filter(lambda (prn, dopp, el): (dopp != None), dopps)
      return dopps
    else:
      return None

if __name__ == "__main__":
  """
  Main method for using almanac info
  """
  import argparse
  parser = argparse.ArgumentParser(description='Swift Nav SBP log converter.')
  parser.add_argument('-w', '--week', nargs=1, type=int, default=[None],
                      help='Specify the week number expected.')
  parser.add_argument('-t', '--tow', nargs=1, type=int, default=[None],
                      help='Specify the time of week.')
  parser.add_argument("-p", "--position", nargs=3, type=float, default=WPR,
                      help="Specify ECE3 coordinates")
  args = parser.parse_args()
  alm = AlmanacData(args.week[0], args.tow[0])
  print "Dopplers: "
  print alm.get_dopps(args.tow[0], args.position[0:3])

