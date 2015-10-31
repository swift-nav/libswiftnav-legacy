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
      include_dirs = [np.get_include(), '.'],
      extra_compile_args = ['-O0', '-g'],
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
    'swiftnav.ephemeris',
    'swiftnav.linear_algebra',
    'swiftnav.amb_kf',
    'swiftnav.gpstime',
    'swiftnav.observation',
    'swiftnav.dgnss_management',
    'swiftnav.ambiguity_test'
  ]
  extensions = [make_extension(name) for name in ext_names]
  setup_args['ext_modules'] = extensions
  setup_args['cmdclass'] = {'build_ext': build_ext}
  setup(**setup_args)
