# Copyright (C) 2016 Swift Navigation Inc.
# Contact: Pasi Miettinen <pmiettinen@exafore.com>
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.


from ctypes import *

try:
    lib = cdll.LoadLibrary("./libl2cbitstream.so")
    r"""library: tries to load default library name from current working dir
    """
except:
    lib = None


def loadlibl2c(libfullpath):
    r"""Loads the given library for the module functions to use

    Parameters
    ----------
    libfullpath: string

    Returns
    -------
    None

    Raises
    ------
    Exception
    """
    global lib
    lib = cdll.LoadLibrary(libfullpath)


def getl2cmessage():
    r"""Returns a ctypes string buffer containing L2C message frame.

    Buffer is 38 bytes containing 300 bit message which is rightaligned ie.
    the leftmost four bits are always zero

    Parameters
    ----------
    None

    Returns
    -------
    buffer
        ctypes 38 byte size string buffer. In case of error 0 length buffer.

    Raises
    ------
    Exception
    """

    buf = create_string_buffer(lib.getL2CMessageLength())
    if lib.getL2CMessage(buf):
        return buf
    else:
        return create_string_buffer(0)


def writel2c_to_file(filename, messageamount):
    r"""Writes given amount of messages to a file.

    Parameters
    ----------
    filename: string

    messageamount: int
        How many messages should be written into the file

    Returns
    -------
    int
        amount of fetched messages

    Raises
    ------
    Exception
    """

    l2cfile = open(filename, "wb")
    ret = lib.writeL2CToFile(l2cfile.fileno(), messageamount)
    l2cfile.flush()
    l2cfile.close()
    return ret