
import struct

# Automatically generated from sbp.yaml with generate.py, do not hand edit!



SBP_GPS_TIME = 0x0100
class GPSTime:
  """
  SBP class for message GPS_TIME (0x0100)

  GPS Time.

  """

  def __init__(self, d):
    self.from_binary(d)

  def from_binary(self, d):
    (
      self.wn,
      self.tow,
      self.ns,
      self.flags,
    ) = struct.unpack('<HIiB', d)

  def to_binary(self, d):
    return struct.pack('<HIiB', (
      self.wn,
      self.tow,
      self.ns,
      self.flags,
    ))


SBP_DOPS = 0x0206
class Dops:
  """
  SBP class for message DOPS (0x0206)

  Dilution of Precision.

  """

  def __init__(self, d):
    self.from_binary(d)

  def from_binary(self, d):
    (
      self.tow,
      self.gdop,
      self.pdop,
      self.tdop,
      self.hdop,
      self.vdop,
    ) = struct.unpack('<IHHHHH', d)

  def to_binary(self, d):
    return struct.pack('<IHHHHH', (
      self.tow,
      self.gdop,
      self.pdop,
      self.tdop,
      self.hdop,
      self.vdop,
    ))


SBP_POS_ECEF = 0x0200
class PosECEF:
  """
  SBP class for message POS_ECEF (0x0200)

  Position solution in absolute Earth Centered Earth Fixed (ECEF) coordinates.

  """

  def __init__(self, d):
    self.from_binary(d)

  def from_binary(self, d):
    (
      self.tow,
      self.x,
      self.y,
      self.z,
      self.accuracy,
      self.n_sats,
      self.flags,
    ) = struct.unpack('<IdddHBB', d)

  def to_binary(self, d):
    return struct.pack('<IdddHBB', (
      self.tow,
      self.x,
      self.y,
      self.z,
      self.accuracy,
      self.n_sats,
      self.flags,
    ))


SBP_POS_LLH = 0x0201
class PosLLH:
  """
  SBP class for message POS_LLH (0x0201)

  Geodetic position solution.

  """

  def __init__(self, d):
    self.from_binary(d)

  def from_binary(self, d):
    (
      self.tow,
      self.lat,
      self.lon,
      self.height,
      self.h_accuracy,
      self.v_accuracy,
      self.n_sats,
      self.flags,
    ) = struct.unpack('<IdddHHBB', d)

  def to_binary(self, d):
    return struct.pack('<IdddHHBB', (
      self.tow,
      self.lat,
      self.lon,
      self.height,
      self.h_accuracy,
      self.v_accuracy,
      self.n_sats,
      self.flags,
    ))


SBP_BASELINE_ECEF = 0x0202
class BaselineECEF:
  """
  SBP class for message BASELINE_ECEF (0x0202)

  Baseline in Earth Centered Earth Fixed (ECEF) coordinates.

  """

  def __init__(self, d):
    self.from_binary(d)

  def from_binary(self, d):
    (
      self.tow,
      self.x,
      self.y,
      self.z,
      self.accuracy,
      self.n_sats,
      self.flags,
    ) = struct.unpack('<IiiiHBB', d)

  def to_binary(self, d):
    return struct.pack('<IiiiHBB', (
      self.tow,
      self.x,
      self.y,
      self.z,
      self.accuracy,
      self.n_sats,
      self.flags,
    ))


SBP_BASELINE_NED = 0x0203
class BaselineNED:
  """
  SBP class for message BASELINE_NED (0x0203)

  Baseline in local North East Down (NED) coordinates.

  """

  def __init__(self, d):
    self.from_binary(d)

  def from_binary(self, d):
    (
      self.tow,
      self.n,
      self.e,
      self.d,
      self.h_accuracy,
      self.v_accuracy,
      self.n_sats,
      self.flags,
    ) = struct.unpack('<IiiiHHBB', d)

  def to_binary(self, d):
    return struct.pack('<IiiiHHBB', (
      self.tow,
      self.n,
      self.e,
      self.d,
      self.h_accuracy,
      self.v_accuracy,
      self.n_sats,
      self.flags,
    ))


SBP_VEL_ECEF = 0x0204
class VelECEF:
  """
  SBP class for message VEL_ECEF (0x0204)

  Velocity in Earth Centered Earth Fixed (ECEF) coordinates.

  """

  def __init__(self, d):
    self.from_binary(d)

  def from_binary(self, d):
    (
      self.tow,
      self.x,
      self.y,
      self.z,
      self.accuracy,
      self.n_sats,
      self.flags,
    ) = struct.unpack('<IiiiHBB', d)

  def to_binary(self, d):
    return struct.pack('<IiiiHBB', (
      self.tow,
      self.x,
      self.y,
      self.z,
      self.accuracy,
      self.n_sats,
      self.flags,
    ))


SBP_VEL_NED = 0x0205
class VelNED:
  """
  SBP class for message VEL_NED (0x0205)

  Velocity in local North East Down (NED) coordinates.

  """

  def __init__(self, d):
    self.from_binary(d)

  def from_binary(self, d):
    (
      self.tow,
      self.n,
      self.e,
      self.d,
      self.h_accuracy,
      self.v_accuracy,
      self.n_sats,
      self.flags,
    ) = struct.unpack('<IiiiHHBB', d)

  def to_binary(self, d):
    return struct.pack('<IiiiHHBB', (
      self.tow,
      self.n,
      self.e,
      self.d,
      self.h_accuracy,
      self.v_accuracy,
      self.n_sats,
      self.flags,
    ))



msg_classses = {
  0x0100: GPSTime,
  0x0206: Dops,
  0x0200: PosECEF,
  0x0201: PosLLH,
  0x0202: BaselineECEF,
  0x0203: BaselineNED,
  0x0204: VelECEF,
  0x0205: VelNED,
}

def sbp_decode(t, d):
  return msg_classses[t](d)
