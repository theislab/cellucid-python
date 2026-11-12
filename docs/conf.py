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
    release = version = "0.9.8"  # CELLUCID_VERSION

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
    "sphinx_sitemap",
    "sphinxext.opengraph",
    "sphinx_llms_txt",
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

# mudata is an optional extra (``cellucid[mudata]``) the docs build does not
# install; the signatures name ``mudata.MuData`` only for type checkers.
autodoc_mock_imports = ["mudata"]

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
# The app's favicon (cellucid/assets/img/favicon.svg): the mark, not the wordmark.
html_favicon = "_static/favicon.svg"

html_theme_options = {
    # Note: analytics config removed due to pydata-sphinx-theme bug causing duplicate scripts
    # Using custom template in _templates/layout.html instead
    # The wordmark already says "Cellucid"; no text beside it.
    "logo": {
        "alt_text": "Cellucid",
    },
    "header_links_before_dropdown": 5,
    "navbar_align": "left",
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "icon_links": [
        {
            "name": "Cellucid App",
            # Navigational link: the canonical browser host is the `www` form,
            # which is what `cellucid/CNAME` publishes and what the app's own
            # <link rel="canonical"> declares. The bare apex is reserved for
            # frozen identifiers (schema $id, figure provenance) and is not
            # used for links. See docs/contributing/web_app.md.
            "url": "https://www.cellucid.com",
            "icon": "fa-solid fa-globe",
        },
    ],
    "footer_start": ["copyright"],
    "footer_end": ["sphinx-version"],
    # The home page (index.md) is a landing page laid out at full width: no
    # "On this page" column, and no section navigation (see html_sidebars).
    "secondary_sidebar_items": {"**": ["page-toc"], "index": []},
    # Navigation
    "navigation_with_keys": True,
    "show_nav_level": 1,
    "show_toc_level": 2,
}

html_context = {
    "default_mode": "light",
}

html_static_path = ["_static"]
# Copied verbatim to the root of the HTML output. Read the Docs serves
# `/robots.txt` from the root of the default version (`latest`), so the
# `_extra/robots.txt` of the `latest` build is the one crawlers see.
html_extra_path = ["_extra"]
html_css_files = ["custom.css", "home.css"]

# Sidebar: section navigation on every page but the home page.
html_sidebars = {
    "**": ["sidebar-nav-bs"],
    "index": [],
}

# -- Canonical URL and sitemap ---------------------------------------------
# The documentation is served from the custom domain docs.cellucid.com (Read
# the Docs canonical domain). `html_baseurl` makes every page declare
# <link rel="canonical"> on that host, and sphinx-sitemap writes the page list
# to https://docs.cellucid.com/en/latest/sitemap.xml, which `_extra/robots.txt`
# advertises.
html_baseurl = "https://docs.cellucid.com/en/latest/"
sitemap_url_scheme = "{link}"
sitemap_excludes = ["search.html", "genindex.html", "py-modindex.html"]

# -- Search results and link previews -----------------------------------------
# sphinxext-opengraph writes each page's <meta name="description"> (from its
# opening text, or its `myst: html_meta: description` front matter) and its
# Open Graph tags,
# so a search result has a snippet and a link to any page shared on X,
# LinkedIn, Slack, Bluesky or Mastodon shows a card. The card image is the web
# app's own (cellucid/assets/img/og-preview.jpg, 1200 x 630, pinned at that URL
# by cellucid's tests/seo-head-contract.test.mjs). og:url is the page's
# canonical address on html_baseurl.
ogp_site_url = html_baseurl
ogp_site_name = "Cellucid"
ogp_image = "https://www.cellucid.com/assets/img/og-preview.jpg"
ogp_image_alt = (
    "A 3D view of a single-cell atlas in Cellucid: thousands of cells drawn as "
    "spheres, colored by cell type, receding into depth."
)
ogp_social_cards = {"enable": False}
ogp_custom_meta_tags = [
    '<meta property="og:image:width" content="1200" />',
    '<meta property="og:image:height" content="630" />',
    '<meta name="twitter:card" content="summary_large_image" />',
]

# Every page's <title> names the guide it belongs to: the Python, R and web
# app guides share page names ("Installation", "Testing and CI"), and a title
# that differs only by its URL is a duplicate to a search engine. The home page
# is "Cellucid documentation", not "Cellucid — Cellucid".
TITLE_SUFFIXES = (
    ("user_guide/web_app/", " — Cellucid web app"),
    ("user_guide/python_package/", " — Cellucid Python package"),
    ("user_guide/r_package/", " — Cellucid R package"),
)
HOME_TITLE_SUFFIX = " documentation: interactive single-cell data visualization"


def _title_suffix(app, pagename, templatename, context, doctree):
    # Read by the htmltitle block of _templates/layout.html; the theme's own
    # suffix is a template variable (" — " + html_title), out of a hook's reach.
    context["cellucid_title_suffix"] = f" — {app.config.html_title}"
    if pagename == app.config.root_doc:
        context["cellucid_title_suffix"] = HOME_TITLE_SUFFIX
        return
    for prefix, suffix in TITLE_SUFFIXES:
        if pagename.startswith(prefix):
            context["cellucid_title_suffix"] = suffix
            return


# -- llms.txt -----------------------------------------------------------------
# sphinx-llms-txt writes /en/latest/llms.txt (llmstxt.org: this summary, then
# every page with a link to its Markdown source) and llms-full.txt (every page's
# source in one file) for AI assistants. The summary is the one the web app's
# llms.txt and structured data carry.
llms_txt_title = "Cellucid documentation"
llms_txt_summary = (
    "Cellucid is a free, open-source web app for interactive single-cell data "
    "visualization. It renders millions of cells on the GPU in the browser, in "
    "1D, 2D and 3D embeddings such as UMAP, colored by gene expression or cell "
    "metadata, and reads AnnData, MuData, Loom, 10x Cell Ranger, MatrixMarket "
    "and Seurat or SingleCellExperiment files in place, without uploading them. "
    "The cellucid Python package (Scanpy/AnnData, Jupyter) and R package "
    "(Seurat, SingleCellExperiment) prepare data for it and show data in it. "
    "The app runs at https://www.cellucid.com/."
)
llms_txt_exclude = ["search", "genindex", "py-modindex"]


# -- Source configuration ----------------------------------------------------
templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "**.ipynb_checkpoints",
    "Thumbs.db",
    ".DS_Store",
    # Authoring material for the documentation itself, not documentation. Both
    # live inside the source tree so they sit beside what they describe, and
    # both are excluded so Sphinx neither reads them nor reports them as pages
    # missing from a toctree — which, under `-W`, would fail the build.
    #
    # Exclusion here does not exempt them from `tests/`: the Python contract
    # suite walks `docs/**/*.md` directly. See
    # `_tooling/screenshots/capture.mjs` and the checklist's conventions page.
    "_tooling",
    "user_guide/web_app/r_screenshot_checklist",
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


# -- Page descriptions without the reading card ------------------------------
# Most guide pages open with a reading card: **Audience:**, **Time:**, then
# **Goal:**, **Prerequisites:** or **What you'll learn/get/do:** and a list.
# Read as a description it is the same boilerplate on 172 pages. Before
# sphinxext-opengraph reads the page (priority 500), the card leaves the
# doctree it is given; the page itself is already written, so only the
# description changes. A Goal becomes the description's first sentence.
PREAMBLE_LABELS = frozenset(
    label.replace("'", "’") if curly else label
    for label in (
        "Audience:",
        "Time:",
        "Goal:",
        "Prerequisites:",
        "What you'll learn:",
        "What you'll get:",
        "What you'll do:",
    )
    for curly in (False, True)
)


def _drop_reading_card(app, pagename, templatename, context, doctree):
    from docutils import nodes

    if doctree is None:
        return
    card = next(iter(doctree.findall(nodes.paragraph)), None)
    if card is None:
        return
    labels = [child.astext().strip() for child in card.children if isinstance(child, nodes.strong)]
    if not labels or labels[0] != "Audience:" or not set(labels) <= PREAMBLE_LABELS:
        return
    goal, label = [], None
    for child in card.children:
        if isinstance(child, nodes.strong):
            label = child.astext().strip()
        elif label == "Goal:" and not isinstance(child, nodes.raw):
            goal.append(child.astext())
    following = card.next_node(descend=False, siblings=True)
    if labels[-1] != "Goal:" and isinstance(following, nodes.bullet_list):
        following.parent.remove(following)
    text = " ".join("".join(goal).split())
    if text:
        card.replace_self(nodes.paragraph(text=text[0].upper() + text[1:].rstrip(".") + "."))
    else:
        card.parent.remove(card)


def setup(app):
    app.connect("html-page-context", _drop_reading_card, priority=400)
    # After sphinxext-opengraph (priority 500), whose og:title is the page title.
    app.connect("html-page-context", _title_suffix, priority=600)
