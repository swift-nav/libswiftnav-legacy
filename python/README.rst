==================
libswiftnav-python
==================

Python bindings to the `libswiftnav
<http://github.com/swift-nav/libswiftnav>`_ GNSS library. Full
documentation is available online at
http://docs.swift-nav.com/libswiftnav-python.

Installation
============

Obtaining the source
--------------------

The libswiftnav-python source and release tarballs are available from
GitHub, https://github.com/swift-nav/libswiftnav-python. The latest
development version of libswiftnav-python can be cloned from github
using this command::

   $ git clone git://github.com/swift-nav/libswiftnav-python.git

Requirements
--------------------

libswiftnav-python requires libswiftnav C libary available at
`libswiftnav <https://github.com/swift-nav/libswiftnav>`_. If you've
checkout this repository, get this as a submodule::

    $ git submodule init
    $ git submodule update

The Python dependencies are included in ``requirements.txt``, which
you can install via::

    $ (sudo) pip install -r requirements

Building and Installing
-----------------------

This library uses `Swig <https://github.com/swig/swig>`_ to generate
Python bindings. You can install Swig from source or using your
platform's package manager of choice. To install libswiftnav-python
(from the root of the source tree)::

    $ python setup.py build
    $ (sudo) python setup.py install
