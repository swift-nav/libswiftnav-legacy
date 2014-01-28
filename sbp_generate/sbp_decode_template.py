
import struct

# Automatically generated from sbp.yaml with generate.py, do not hand edit!

msg_classses = {
((*- for m in msgs *))
  ((('0x%04X'|format(m.id)))): ((( m.name | classnameify ))),
((*- endfor *))
}

def sbp_decode(t, d):
  return msg_classses[t](d)

((* for m in msgs *))
class ((( m.name | classnameify ))):
  """
  SBP class for message (((m.name))) ((('(0x%04X)'|format(m.id))))

  (((m.desc)))
  """

  def __init__(self, d):
    self.from_binary(d)

  def from_binary(self, d):
    (
    ((*- for f in m.fields *))
      self.(((f.name))),
    ((*- endfor *))
    ) = struct.unpack('((( m.fields | pystruct )))', d)

  def to_binary(self, d):
    return struct.pack('((( m.fields | pystruct )))', (
    ((*- for f in m.fields *))
      self.(((f.name))),
    ((*- endfor *))
    ))

((* endfor *))
