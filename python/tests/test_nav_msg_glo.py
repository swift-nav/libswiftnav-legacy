

from swiftnav.nav_msg_glo import NavMsgGlo
from swiftnav.ephemeris import Ephemeris
import numpy


def test_imports():
  """Verify that distributed packages survive setuptools installation.

  """
  pass

GLO_STRINGS = numpy.asanyarray([[1, 1, 1],  # dummy words
                                # 01074396999b05c3a850b5
                                [0x010743, 0x96999b05, 0xc3a850b5],
                                # 021760a5256204d9c15f66
                                [0x021760, 0xa5256204, 0xd9c15f66],
                                # 0380269d60899a6d0e3123
                                [0x038026, 0x6d0e3123, 0x9d60899a],
                                # 04865d1cc0000000344918
                                [0x04865d, 0x00344918, 0x1cc00000],
                                # 050d100000000340000895
                                [0x050d10, 0x3, 0x40000895]
                                ],
                               dtype=numpy.dtype('>u4'))
GLO_STRING_BITS = [numpy.unpackbits(s.view(numpy.uint8))[-85:]
                   for s in GLO_STRINGS]

GLO_TM = numpy.asanyarray([0x3E375096], dtype=numpy.dtype('>u4'))
GLO_TM_BITS = numpy.unpackbits(GLO_TM.view(numpy.uint8))[-30:]


def test_init():
  '''
  Object construction test
  '''
  msg = NavMsgGlo()
  assert msg.isDecodeDone() == False
  assert msg.getTow() == -1.


def test_update():
  '''
  Object update test
  '''
  msg = NavMsgGlo()
  assert msg.update(0) == -1
  assert msg.update(0) == -1
  assert msg.update(0) == -1
  assert msg.update(1) == -1

  assert msg.isDecodeDone() == False
  assert msg.getTow() == -1.


def test_decode():
  '''
  Ephemeris and ToW decoding test.
  '''
  msg = NavMsgGlo()
  n_strings = 0
  n_tow = 0
  tow = -1
  e = Ephemeris()
  print GLO_STRING_BITS
  for bit_string in GLO_STRING_BITS:
    print bit_string
    for bit in bit_string:
      assert msg.update(bit ^ 1) == -1
      if msg.update(bit) == 1:
        n_strings += 1
        if msg.updateEphemeris(e) == 1:
          n_tow += 1
          tow = msg.getTow()

    # now pass time mark bit by bit to receiver (MSB first), no line code
    # needed
    for bit in GLO_TM_BITS:
      assert msg.update(bit) == -1

  assert n_strings == 5
  assert n_tow == 1
  assert isinstance(tow, float)
  assert tow == 262707.0
  assert e.toe['tow'] == 262707.0 - 10.
