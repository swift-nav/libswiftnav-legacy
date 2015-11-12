#!/usr/bin/env python

from version import get_git_version
import os
import sys

try:
  from Cython.Build import cythonize
  from Cython.Distutils import build_ext
  from setuptools import Command, Extension, setup
  from setuptools.command.install import install
  from subprocess import call
  from distutils.command.build import build
  import numpy as np
except ImportError as exc_info:
  print "Oops! %s." % exc_info.message
  sys.exit(1)

# Get the README and requirements
CWD = os.path.abspath(os.path.dirname(__file__))
with open(CWD + '/README.rst') as f:
  readme = f.read()

with open(CWD + '/requirements.txt') as f:
  INSTALL_REQUIRES = [i.strip() for i in f.readlines()]

PLATFORMS = ['linux', 'osx', 'win32']

class SwiftCInstall(build):
  """Compile and install libswiftnav C bindings. Why does someone have
  to write something like this from scratch?

  """

  def run(self):
    """Before building extensions, build and copy libswiftnav.

    """
    print "Building libswiftnav C library!"
    path = os.path.join(CWD, '')
    cmd = ['cd ../',
           'mkdir -v -p build/ && cd build/ && cmake ../ && make',
           'cd %s' % path]
    lib = "libswiftnav-static.a"
    target_files = []
    print "Produced...%s\n" % target_files
    def compile():
      print '*' * 80
      os.system(";\n".join(cmd))
      print '*' * 80
    self.execute(compile, [], '\nCompiling libswiftnav C!\n')
    for root, dirs, files in os.walk(r'../build/'):
      for name in files:
        if name == lib:
          target_files = [os.path.abspath(os.path.join(root, name))]
    # copy resulting tool to library build folder
    self.build_lib = 'build/lib'
    self.mkpath(self.build_lib)
    if not self.dry_run:
      for target in target_files:
        print "\nCopying %s to %s.\n" % (target, self.build_lib)
        self.copy_file(target, self.build_lib)
    build.run(self)
    print("\n\n\n\nSuccessfully built libswiftnav C libraries!\n\n\n\n")

SETUP_ARGS = dict(name='swiftnav',
                  version=get_git_version(),
                  description='Python bindings to the libswiftnav library',
                  long_description=readme,
                  license='LGPLv3',
                  url='https://github.com/swift-nav/libswiftnav',
                  author='Swift Navigation',
                  author_email='dev@swiftnav.com',
                  maintainer='Swift Navigation',
                  maintainer_email='dev@swiftnav.com',
                  install_requires=INSTALL_REQUIRES,
                  platforms=PLATFORMS,
                  use_2to3=False,
                  zip_safe=False,
                  packages=['swiftnav'])

def setup_package():
  """Setup the package and required Cython extensions.

  """
  # Override ARCHFLAGS to build native, not universal on OS X.
  os.environ['ARCHFLAGS'] = ""
  extensions = Extension("*",
                         ["swiftnav/*.pyx"],
                         include_dirs=[np.get_include(), '.'],
                         extra_compile_args=['-O3', '-g', '-Wno-unused-function'],
                         extra_link_args=['-g'],
                         libraries=['m', 'swiftnav'])
  SETUP_ARGS['ext_modules'] = cythonize(extensions)
  SETUP_ARGS['cmdclass'] = {'build': SwiftCInstall,
                            'build_ext': build_ext}
  setup(**SETUP_ARGS)

if __name__ == "__main__":
  setup_package()
