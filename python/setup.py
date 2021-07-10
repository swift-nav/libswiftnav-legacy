#!/usr/bin/python
from version import get_git_version

setup_args = dict(
  name='swiftnav',
  version=get_git_version(),
  description='Python bindings to the libswiftnav library.',
  license='LGPLv3',
  url='http://www.swiftnav.com',
  author='Swift Navigation Inc.',
  author_email='dev@swiftnav.com',
  maintainer='Swift Navigation',
  maintainer_email='dev@swiftnav.com',
  packages=['swiftnav'],
)

if __name__ == "__main__":
  import numpy as np
  import os
  import sys
  from setuptools import setup, Extension
  try:
    from Cython.Distutils import build_ext
  except:
    print "You don't seem to have Cython installed."
    sys.exit(1)
  os.environ['ARCHFLAGS'] = ""

  # Additional library search directories:
  # - LD_LIBRARY_PATH (if present)
  # - User local enviroment
  library_dirs = []
  # If LD_LIBRARY_PATH has been manually specified, add it to the
  # library search path
  if 'LD_LIBRARY_PATH' in os.environ:
    library_dirs.append(os.environ['LD_LIBRARY_PATH'])
  library_dirs.append(os.path.expanduser('~/.local/lib'))

  # Additional include directories:
  # - Numpy includes
  # - User local enviroment
  # - Current directory
  include_dirs = []
  include_dirs.append(np.get_include())
  include_dirs.append(os.path.expanduser('~/.local/include'))
  include_dirs.append('.')
  # three more includes for travis builds as it does not install libraries
  include_dirs.append('../include/')
  include_dirs.append('../libfec/include/')
  include_dirs.append('../tests/data/l2cbitstream/')

  def make_extension(ext_name):
    ext_path = ext_name.replace('.', os.path.sep) + '.pyx'
    return Extension(
      ext_name, [ext_path],
      include_dirs=include_dirs,
      extra_compile_args=['-O0', '-g'],
      extra_link_args=['-g'],
      libraries=['m', 'swiftnav', 'l2cbitstream'],
      library_dirs=library_dirs,
    )
  ext_names = [
    'swiftnav.edc',
    'swiftnav.signal',
    'swiftnav.coord_system',
    'swiftnav.constants',
    'swiftnav.cnav_msg',
    'swiftnav.nav_msg',
    'swiftnav.nav_msg_glo',
    'swiftnav.pvt',
    'swiftnav.correlate',
    'swiftnav.track',
    'swiftnav.almanac',
    'swiftnav.lambda_',
    'swiftnav.ephemeris',
    'swiftnav.linear_algebra',
    'swiftnav.amb_kf',
    'swiftnav.time',
    'swiftnav.observation',
    'swiftnav.dgnss_management',
    'swiftnav.ambiguity_test',
    'swiftnav.baseline',
    'swiftnav.bits',
    'swiftnav.filter_utils',
    'swiftnav.memory_pool',
    'swiftnav.prns',
    'swiftnav.sats_management',
    'swiftnav.iono',
    'swiftnav.tropo',
    'swiftnav.set',
    'swiftnav.bit_sync',
  ]
  extensions = [make_extension(name) for name in ext_names]
  setup_args['ext_modules'] = extensions
  setup_args['cmdclass'] = {'build_ext': build_ext}
  setup(**setup_args)
