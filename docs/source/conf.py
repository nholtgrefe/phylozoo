# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
from pathlib import Path

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'phylozoo'
copyright = '2025, N. Holtgrefe'
author = 'N. Holtgrefe'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Add the src directory to the Python path for autodoc
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',  # For NumPy-style docstrings
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autosummary',
    'sphinx.ext.mathjax',  # For mathematical notation
    'sphinx.ext.githubpages',  # Adds .nojekyll for GitHub Pages compatibility
    'sphinxcontrib.bibtex',  # BibTeX citations
]

# Optional extensions
try:
    import sphinx_copybutton
    extensions.append('sphinx_copybutton')  # Copy button for code blocks
except ImportError:
    pass  # Optional extension, skip if not available

# BibTeX configuration
bibtex_bibfiles = ['bibliography.bib']  # Relative to source directory
bibtex_default_style = 'unsrt'
bibtex_reference_style = 'author_year'

# Napoleon settings for NumPy-style docstrings
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
    'inherited-members': False,
}

autodoc_mock_imports = []

templates_path = ['_templates']
exclude_patterns = []

language = 'en'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# Read the Docs theme (works on both RTD and GitHub Pages)
try:
    import sphinx_rtd_theme
    html_theme = 'sphinx_rtd_theme'
    # Note: get_html_theme_path() is deprecated in newer versions, but still works
    # The theme is automatically found if installed
    try:
        html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
    except AttributeError:
        # Newer versions don't need this
        html_theme_path = []
except ImportError:
    html_theme = 'default'
    html_theme_path = []

html_static_path = ['_static']

# Theme options
if html_theme == 'sphinx_rtd_theme':
    html_theme_options = {
        'collapse_navigation': False,
        'sticky_navigation': True,
        'navigation_depth': 4,
        'includehidden': True,
        'titles_only': False,
        'logo_only': False,
        'display_version': True,
        'prev_next_buttons_location': 'bottom',
        'style_external_links': True,
        'style_nav_header_background': '#2980B9',
    }

# GitHub Pages compatibility (optional, for future migration)
# Uncomment and update if migrating to GitHub Pages:
# html_baseurl = 'https://yourusername.github.io/phylozoo/'

# Read the Docs automatically handles this, but good to have for local builds
html_context = {
    'display_github': True,  # Enable "Edit on GitHub" link
    'github_user': 'yourusername',  # Update with your GitHub username
    'github_repo': 'phylozoo',
    'github_version': 'main',
    'conf_py_path': '/docs/source/',
}

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'networkx': ('https://networkx.org/documentation/stable/', None),
}
