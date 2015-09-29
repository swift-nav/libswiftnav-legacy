libswiftnav development
=======================

Follow these instructions if you intend to make changes to libswiftnav.

Tools needed:
 - doxygen
 - convert from ImageMagick
 - pip install gcovr diff-cover

To get started, run::

  ./checks/setup-hooks.sh [DIR]

from within the libswiftnav root directory. The default build-dir will be `build`.

This makes the directory if it doesn't already exist, runs cmake with `Coverage` mode enabled
inside it, and installs a git pre-commit hook that runs style and code-coverage checkers
in this build directory. These checks require doxygen, gcovr, and diff-cover

To manually run the coverage task, use `make check-coverage`. For syntax, use `make check-style`.

Issues
======

**"No lines with coverage information in this diff."**

You may see this in the output from diff-cover.

Did you set up your build with::

  cmake -DCMAKE_BUILD_TYPE=Coverage ..

or with the :code:`setup-hooks.sh` script? The :code:`Coverage` flag is
required.  If you are manually running the coverage script, make sure you run
it from within your build directory.
