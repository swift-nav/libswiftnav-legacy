import sys, os

# Append the project root to the path so our packages can be found.
sys.path.insert(0, os.path.abspath('..'))
# Append the Sphinx extensions directory to the path.
sys.path.insert(0, os.path.abspath('extensions'))

from setup import setup_args

# -- General configuration -----------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#needs_sphinx = '1.0'

# Sphinx extension modules.
extensions = [
  'sphinx.ext.autosummary',
  'sphinx.ext.autodoc',
  'sphinx.ext.intersphinx',
  'sphinx.ext.todo',
  'sphinx.ext.coverage',
  'sphinx.ext.pngmath',
  'sphinx.ext.viewcode',
  'matplotlib.sphinxext.only_directives',
  'matplotlib.sphinxext.plot_directive',
  'ipython_directive',
  'matplotlib.sphinxext.ipython_console_highlighting',
  'numpydoc',
  'sphinxcontrib.doxylink.doxylink'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The encoding of source files.
#source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'libswiftnav-python'
copyright = u'2012, Swift Navigation Inc.'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
# Extract the version from our setup.py file
version = setup_args['version']
# The full version, including alpha/beta/rc tags.
release = version

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build', '_templates', 'extensions']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# -- Options for HTML output ---------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'default'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#html_theme_options = {}

# Add any paths that contain custom themes here, relative to this directory.
#html_theme_path = []

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
#html_title = None

# A shorter title for the navigation bar.  Default is the same as html_title.
#html_short_title = None

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
#html_logo = None

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#html_favicon = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
#html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
#html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
#html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
#html_additional_pages = {}

# If false, no module index is generated.
#html_domain_indices = True

# If false, no index is generated.
#html_use_index = True

# If true, the index is split into individual pages for each letter.
#html_split_index = False

# If true, links to the reST sources are added to the pages.
#html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
#html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
#html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
#html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
#html_file_suffix = None

# Output file base name for HTML help builder.
htmlhelp_basename = 'libswiftnav-pythondoc'

# -- Misc. Options -------------------------------------------------------------

# Configure intersphinx.
intersphinx_mapping = {
  'python': ('http://docs.python.org/', None),
  'scipy': ('http://docs.scipy.org/doc/scipy/reference/', None),
  'matplotlib': ('http://matplotlib.sourceforge.net/', None),
  'numpy': ('http://docs.scipy.org/doc/numpy/', None)
}

# Grab the latest tagfile from the libswiftnav Doxygen docs.
print "Downloading libswiftnav.tag"
import urllib
urllib.urlretrieve('http://docs.swift-nav.com/libswiftnav/libswiftnav.tag',
                   'libswiftnav.tag')

doxylink = {
  'libswiftnav' : ('libswiftnav.tag', 'http://docs.swift-nav.com/libswiftnav/')
}

# Generate stub files from autosummary.
autosummary_generate = True

# Run our patched sphinx-apidoc script that will document .pyx files.
import apidoc
apidoc.main(['fakin it', '--force', '--no-toc', '-o', '.', '../swiftnav'])

