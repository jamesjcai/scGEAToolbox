# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
# sys.path.insert(0, os.path.abspath('.'))


# https://sphinx-rtd-theme.readthedocs.io/en/stable/
import sphinx_rtd_theme


# -- Project information -----------------------------------------------------

project = 'scGEAToolbox'
copyright = '2026, James Cai'
author = 'James Cai'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx_rtd_theme',
]

# https://github.com/sphinx-contrib/matlabdomain  'sphinxcontrib.matlab', 'sphinx.ext.autodoc'


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'


# No custom static files or local themes: both '_static' and '_themes' were
# declared here without existing, which Sphinx warns about on every build.
# Re-add html_static_path when there is actually something to put in it.
html_static_path = []
