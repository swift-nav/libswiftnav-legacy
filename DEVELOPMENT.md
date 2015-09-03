libswiftnav development
=======================

Follow these instructions if you intend to make changes to libswiftnav.

Tools needed:
 * doxygen
 * pip install gcovr diff-cover

To get started, run::
  ./checks/setup-hooks.sh [DIR]

from within the libswiftnav root directory. The default build-dir will be `build`.

This makes the directory if it doesn't already exist, runs cmake with `Coverage` mode enabled
inside it, and installs a git pre-commit hook that runs style and code-coverage checkers
in this build directory. These checks require doxygen, gcovr, and diff-cover

To manually run the coverage task, use `make coverage`. For syntax, use `make style`.
