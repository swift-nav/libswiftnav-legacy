#!/usr/bin/python
from version import get_git_version

setup_args = dict(
  name = 'swiftnav',
  version = get_git_version(),
  description = 'Python bindings to the libswiftnav library.',
  license = 'LGPLv3',
  url = 'http://www.swiftnav.com',
  author = 'Swift Navigation Inc.',
  author_email = 'dev@swiftnav.com',
  maintainer = 'Swift Navigation',
  maintainer_email = 'dev@swiftnav.com',
  packages = ['swiftnav'],
)

if __name__ == "__main__":
  import numpy as np
  import os, sys
  from setuptools import setup, Extension
  try:
    from Cython.Distutils import build_ext
  except:
    print "You don't seem to have Cython installed."
    sys.exit(1)
  os.environ['ARCHFLAGS'] = ""
  def make_extension(ext_name):
    ext_path = ext_name.replace('.', os.path.sep) + '.pyx'
    return Extension(
      ext_name, [ext_path],
      include_dirs = [np.get_include(), '.', '../include/'],
      extra_compile_args = ['-O0', '-g'],
      extra_link_args = ['-g'],
      libraries = ['m', 'swiftnav'],
    )
  ext_names = [
    'swiftnav.edc',
    'swiftnav.signal',
    'swiftnav.coord_system',
    'swiftnav.constants',
    'swiftnav.nav_msg',
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
    'swiftnav.tropo',
    'swiftnav.set',
    'swiftnav.bit_sync',
  ]
  extensions = [make_extension(name) for name in ext_names]
  setup_args['ext_modules'] = extensions
  setup_args['cmdclass'] = {'build_ext': build_ext}
  setup(**setup_args)
