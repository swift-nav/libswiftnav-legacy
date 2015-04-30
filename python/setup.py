#!/usr/bin/env python

import os
import sys

try:
  from setuptools import Command, Extension, setup
  from setuptools.command.install import install
  from distutils.command.build import build
  from distutils.command.build_py import build_py
  from distutils.command.clean import clean
except ImportError as exc_info:
  print "Fucking shit! %s." % exc_info.message
  sys.exit(1)

# Package version and platform info for Pyi
VERSION = "0.13"
PLATFORMS = ['linux', 'osx']

# Get the README and requirements
CWD = os.path.abspath(os.path.dirname(__file__))
with open(CWD + '/README.rst') as f:
  readme = f.read()

with open(CWD + '/requirements.txt') as f:
  INSTALL_REQUIRES = [i.strip() for i in f.readlines()]


# Swig and build flags.
SWIG_INCLUDE_PATH = '../include/libswiftnav/'
SWIG_OPTS = ['-modern', '-Wall', '-O', '-I%s' % SWIG_INCLUDE_PATH, '-castmode']
EXTRA_COMPILE_ARGS = ['-O3', '-g', '-Wno-unused-function']
SWIG_PATH = '../swig/'
LIBRARY_DIRS = ['build/lib/']


def get_extensions(path=SWIG_PATH):
  """Create a C extension for Python, given an extension name.

  """
  exts = []
  modules = []
  fs = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
  print "Copying c and py files from %s." % path
  for file_name in fs:
    name, ftype = file_name.split('.')
    if ftype != 'i' or name == "numpy" or name == "swiftnav_ext":
      print "Skipping file %s." % file_name
      continue
    modules.append(name)
    args = EXTRA_COMPILE_ARGS
    exts.append(Extension('_%s' % name,
                          [path + file_name],
                          extra_compile_args=args,
                          swig_opts=SWIG_OPTS,
                          include_dirs=[SWIG_INCLUDE_PATH, '.'],
                          library_dirs=LIBRARY_DIRS,
                          extra_link_args=['-g'],
                          libraries=['m', 'swiftnav']))
  return (modules, exts)


class SwiftPyBuild(build_py):
  """Copies over swig-generated files.

  """

  def run(self):
    if not self.dry_run:
      path = SWIG_PATH
      fs = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
      for target in fs:
        name, ftype = target.split(".")
        self.copy_file(SWIG_PATH + target, 'swiftnav/')
      build_py.run(self)


class SwiftClean(clean):
  """Cleans swig-generated files.

  """

  def run(self):
    if not self.dry_run:
      paths = ['swiftnav/', SWIG_PATH]
      for path in paths:
        fs = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
        for target in fs:
          name, ftype = target.split(".")
          if name == '__init__' or ftype == 'i':
            continue
          else:
            os.remove(path + target)
      clean.run(self)


class SwiftCInstall(build):
  """Compile and install libswiftnav C bindings. Why does someone have
  to write something like this from scratch?

  """

  def run(self):
    """Before building extensions, build and copy libswiftnav.

    """
    print "Building libswiftnav C library!"
    path = os.path.join(CWD, '../')
    cmd = ['mkdir -v -p %s/build' % path,
           'cd %s/build/' % path,
           'cmake ../',
           'make',
           'cd %s' % path]
    target_files = [os.path.join(path, 'build/src/libswiftnav-static.a')]
    print "Produced...%s\n" % target_files
    def compile():
      print '*' * 80
      os.system(";\n".join(cmd))
      print '*' * 80
    self.execute(compile, [], '\nCompiling libswiftnav C!\n')
    # copy resulting tool to library build folder
    self.build_lib = 'build/lib'
    self.mkpath(self.build_lib)
    if not self.dry_run:
      for target in target_files:
        print "\nCopying %s to %s.\n" % (target, self.build_lib)
        if os.path.isfile(target):
          self.copy_file(target, self.build_lib)
        else:
          assert False, "Expected library at %s" % target
    build.run(self)
    print("\n\n\n\nSuccessfully built libswiftnav C libraries!\n\n\n\n")

SETUP_ARGS = dict(name='swiftnav',
                  version=VERSION,
                  description='Python bindings to the libswiftnav library',
                  long_description=readme,
                  license='LGPLv3',
                  url='https://github.com/swift-nav/libswiftnav-python',
                  author='Swift Navigation',
                  author_email='dev@swiftnav.com',
                  maintainer='Swift Navigation',
                  maintainer_email='dev@swiftnav.com',
                  install_requires=INSTALL_REQUIRES,
                  platforms=PLATFORMS,
                  use_2to3=False,
                  zip_safe=False,
                  packages=['swiftnav'],
                  py_modules=['swiftnav'])


def setup_package():
  """Setup the package and required Cython extensions.

  """
  # Override ARCHFLAGS to build native, not universal on OS X.
  os.environ['ARCHFLAGS'] = ""
  modules, extensions = get_extensions()
  SETUP_ARGS['ext_modules'] = extensions
  SETUP_ARGS['py_modules'] = modules
  SETUP_ARGS['cmdclass'] = {'build': SwiftCInstall,
                            'build_py': SwiftPyBuild,
                            'clean': SwiftClean}
  setup(**SETUP_ARGS)

if __name__ == "__main__":
  setup_package()
