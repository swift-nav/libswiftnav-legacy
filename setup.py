#!/usr/bin/env python

"""
============
Installation
============

Requirements
============

libswiftnav-python has the following requirements, they must be installed
before building and installing libswiftnav-python:

- `Python <http://www.python.org/>`_ 2.6 or 2.7

- `Numpy <http://www.numpy.org/>`_ 1.6 or later

- `Cython <http://www.cython.org>`_ version 0.17 or later

- `libswiftnav <https://github.com/swift-nav/libswiftnav>`_

Installing libswiftnav-python
=============================

Obtaining the source
--------------------

The libswiftnav-python source and release tarballs are available from GitHub,
https://github.com/swift-nav/libswiftnav-python.

The latest development version of libswiftnav-python can be cloned from github
using this command::

   git clone git://github.com/swift-nav/libswiftnav-python.git


Building and Installing
-----------------------

.. note::
  Installation requires `Distribute <http://pypi.python.org/pypi/distribute>`_,
  if your python installation doesn't provide this it will automatically be
  installed.

To install libswiftnav-python (from the root of the source tree)::

    $ python setup.py install

Building documentation
----------------------

.. note::
    The latest version of libswiftnav-python's documentation should be
    available online at http://docs.swift-nav.com/libswiftnav-python.

Building the documentation requires the libswiftnav-python source code and some
additional packages:

    - `Sphinx <http://sphinx.pocoo.org>`_ (and its dependencies) 1.0 or later
    - `numpydoc <http://pypi.python.org/pypi/numpydoc>`_ 0.4 or later
    - `matplotlib <http://matplotlib.org/>`_ 1.1 or later
    - `IPython <http://ipython.org/>`_ 0.13.1 or later
    - `sphinx-doxylink <http://pypi.python.org/pypi/sphinxcontrib-doxylink>`_
      1.3 or later

To build the libswiftnav-python documentation, execute the following commands::

    $ cd docs
    $ make html

The documentation will be built in the ``docs/_build/html`` directory, and can
be read by pointing a web browser to ``docs/_build/html/index.html``.

"""

from version import get_git_version

setup_args = dict(
  name = 'swiftnav',
  version = get_git_version(),
  description = 'Python bindings to the libswiftnav library.',
  license = 'LGPLv3',
  url = 'http://www.swift-nav.com',

  author = 'Swift Navigation Inc.',
  author_email = 'info@swift-nav.com',
  maintainer = 'Fergus Noble',
  maintainer_email = 'fergus@swift-nav.com',

  packages = ['swiftnav'],
)

if __name__ == "__main__":
  # Bootstrap Distribute if the user doesn't have it
  from distribute_setup import use_setuptools
  use_setuptools()

  import numpy as np
  import os, sys

  from setuptools import setup, Extension

  # We'd better have Cython installed, or it's a no-go
  try:
    from Cython.Distutils import build_ext
  except:
    print "You don't seem to have Cython installed."
    sys.exit(1)

  # Override ARCHFLAGS to build native, not universal on OS X.
  os.environ['ARCHFLAGS'] = ""

  def make_extension(ext_name):
    ext_path = ext_name.replace('.', os.path.sep) + '.pyx'
    return Extension(
      ext_name, [ext_path],
      include_dirs = [np.get_include(), '.'],
      extra_compile_args = ['-O3', '-Wall', '-Wno-unused-function'],
      extra_link_args = ['-g'],
      libraries = ['m', 'swiftnav'],
    )

  ext_names = [
    'swiftnav.coord_system',
    'swiftnav.nav_msg',
    'swiftnav.pvt',
    'swiftnav.correlate',
    'swiftnav.track',
    'swiftnav.almanac',
    'swiftnav.lam',
    'swiftnav.float_kf'
  ]

  extensions = [make_extension(name) for name in ext_names]

  setup_args['ext_modules'] = extensions
  setup_args['cmdclass'] = {'build_ext': build_ext}

  setup(**setup_args)

