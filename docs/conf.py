"""Sphinx configuration for Cellucid documentation."""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

# -- Path setup --------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
sys.path.insert(0, str(SRC_DIR))

# -- Project information -----------------------------------------------------
project = "Cellucid"
author = "Kemal Inecik"
copyright = f"{datetime.now().year}, {author}"

# Version
try:
    from cellucid import __version__

    release = version = __version__
except Exception:
    release = version = "0.0.1a2"

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",
    "sphinx_design",
    "myst_nb",
    "sphinx_copybutton",
]

# -- MyST configuration ------------------------------------------------------
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
    "html_admonition",
    "html_image",
    "tasklist",
]
myst_heading_anchors = 3

# -- MyST-NB configuration ---------------------------------------------------
nb_execution_mode = "off"
nb_execution_timeout = 300
nb_output_stderr = "remove"

# -- Autodoc configuration ---------------------------------------------------
autodoc_default_options = {
    "members": True,
    "member-order": "bysource",
    "special-members": "__init__",
    "undoc-members": True,
    "exclude-members": "__weakref__",
    "show-inheritance": True,
}
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_typehints_format = "short"

autodoc_type_aliases = {
    "anndata.AnnData": "anndata.AnnData",
    "AnnData": "anndata.AnnData",
}

# -- Autosummary configuration -----------------------------------------------
autosummary_generate = True

# -- Napoleon configuration --------------------------------------------------
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True

# -- Intersphinx configuration -----------------------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "anndata": ("https://anndata.readthedocs.io/en/latest/", None),
    "scanpy": ("https://scanpy.readthedocs.io/en/stable/", None),
}

# -- HTML output configuration -----------------------------------------------
html_theme = "pydata_sphinx_theme"
html_title = "Cellucid"
html_logo = "_static/cellucid-logo.svg"
html_favicon = "_static/cellucid-logo.svg"

html_theme_options = {
    # Note: analytics config removed due to pydata-sphinx-theme bug causing duplicate scripts
    # Using custom template in _templates/layout.html instead
    "logo": {
        "text": "Cellucid",
    },
    "header_links_before_dropdown": 5,
    "navbar_align": "left",
    "navbar_end": ["theme-switcher", "navbar-icon-links", "version-switcher"],
    "icon_links": [
        {
            "name": "Cellucid App",
            "url": "https://cellucid.com",
            "icon": "fa-solid fa-globe",
        },
        {
            "name": "GitHub: cellucid",
            "url": "https://github.com/theislab/cellucid",
            "icon": "fa-brands fa-github",
        },
        {
            "name": "GitHub: cellucid-annotation",
            "url": "https://github.com/theislab/cellucid-annotation",
            "icon": "fa-solid fa-pen-to-square",
        },
        {
            "name": "GitHub: cellucid-python",
            "url": "https://github.com/theislab/cellucid-python",
            "icon": "fa-brands fa-python",
        },
        {
            "name": "GitHub: cellucid-r",
            "url": "https://github.com/theislab/cellucid-r",
            "icon": "fa-brands fa-r-project",
        },
    ],
    "show_version_warning_banner": True,
    "footer_start": ["copyright"],
    "footer_end": ["sphinx-version"],
    "secondary_sidebar_items": ["page-toc", "edit-this-page"],
    "use_edit_page_button": True,
    # Version switcher
    "switcher": {
        "json_url": "https://cellucid.readthedocs.io/en/latest/_static/switcher.json",
        "version_match": version,
    },
    # Navigation
    "navigation_with_keys": True,
    "show_nav_level": 2,
    "show_toc_level": 2,
}

html_context = {
    "github_user": "theislab",
    "github_repo": "cellucid-python",
    "github_version": "main",
    "doc_path": "docs",
}

html_static_path = ["_static"]
html_css_files = ["custom.css"]

# Sidebar - show navigation on all pages
html_sidebars = {
    "**": ["sidebar-nav-bs"],
}

# -- Source configuration ----------------------------------------------------
templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "**.ipynb_checkpoints",
    "Thumbs.db",
    ".DS_Store",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "myst-nb",
    ".ipynb": "myst-nb",
}

master_doc = "index"

# -- Copybutton configuration ------------------------------------------------
copybutton_prompt_text = r">>> |\.\.\.|\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True
copybutton_remove_prompts = True

# -- Suppress warnings -------------------------------------------------------
suppress_warnings = ["myst.header", "mystnb.unknown_mime_type"]

# -- Nitpicky mode -----------------------------------------------------------
nitpicky = False
nitpick_ignore = [
    ("py:class", "anndata.AnnData"),
    ("py:class", "AnnData"),
]
