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
GitHub, https://github.com/swift-nav/libswiftnav. The latest
development version of libswiftnav-python can be cloned from github
using this command::

   $ git clone git://github.com/swift-nav/libswiftnav.git
   $ cd python

Requirements
--------------------

libswiftnav-python requires libswiftnav C libary available at
`libswiftnav <https://github.com/swift-nav/libswiftnav>`_.

The Python dependencies are included in ``requirements.txt``, which
you can install via::

    $ (sudo) pip install -r requirements

Building and Installing
-----------------------

To install libswiftnav-python (from the root of the source tree)::

    $ python setup.py build
    $ (sudo) python setup.py install

Building Documentation
----------------------

Building the documentation requires the libswiftnav-python source code
and some additional packages:

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
