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

try:
    import sphinx_togglebutton
    extensions.append('sphinx_togglebutton')  # Toggle button for collapsible content
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
    'undoc-members': True,  # Show all members, even if undocumented
    'show-inheritance': True,
    'inherited-members': False,
}

# Suppress common autodoc warnings
# This suppresses warnings about undocumented members, missing imports, etc.
suppress_warnings = [
    'autodoc.import_object',  # Suppress import warnings
    'autodoc',  # Suppress other autodoc warnings
]

autodoc_mock_imports = []

templates_path = ['_templates']
exclude_patterns = []

language = 'en'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# PyData Sphinx Theme (similar to Aesara's documentation style)
try:
    import pydata_sphinx_theme
    html_theme = 'pydata_sphinx_theme'
    # Modern versions don't need html_theme_path - Sphinx finds it automatically
    html_theme_path = []
except ImportError:
    # Fallback to Read the Docs theme if PyData theme not available
    try:
        import sphinx_rtd_theme
        html_theme = 'sphinx_rtd_theme'
        try:
            html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
        except AttributeError:
            html_theme_path = []
    except ImportError:
        html_theme = 'default'
        html_theme_path = []

html_static_path = ['_static']

# PyData Sphinx Theme options (similar to Aesara's style)
if html_theme == 'pydata_sphinx_theme':
    html_theme_options = {
        'github_url': 'https://github.com/nholtgrefe/phylozoo',
        'logo': {
            'text': 'PhyloZoo',
        },
        'use_edit_page_button': True,
        'show_toc_level': 2,
        'navbar_align': 'left',
        'navbar_end': ['theme-switcher', 'navbar-icon-links'],
        'icon_links': [
            {
                'name': 'GitHub',
                'url': 'https://github.com/nholtgrefe/phylozoo',
                'icon': 'fa-brands fa-github',
            },
        ],
        'show_nav_level': 2,
        'navigation_depth': 4,
        'collapse_navigation': False,
    }
    html_context = {
        'github_user': 'nholtgrefe',
        'github_repo': 'phylozoo',
        'github_version': 'main',
        'doc_path': 'docs/source',
    }
elif html_theme == 'sphinx_rtd_theme':
    # Fallback theme options
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

# GitHub context (used by both themes)
# Note: PyData theme uses its own html_context above, but we keep this for RTD fallback
if html_theme != 'pydata_sphinx_theme':
    html_context = {
        'display_github': True,  # Enable "Edit on GitHub" link
        'github_user': 'nholtgrefe',  # Update with your GitHub username
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
