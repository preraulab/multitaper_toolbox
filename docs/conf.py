"""Sphinx configuration for multitaper."""
import os
import sys
from pathlib import Path

project = "multitaper"
author = "Michael J. Prerau Laboratory"
copyright = "2011-present, Michael J. Prerau Laboratory"

# Path setup for matlabdomain -- points to repo root so it finds .m files.
REPO_ROOT = Path(__file__).resolve().parent.parent
matlab_src_dir = str(REPO_ROOT)
primary_domain = "mat"

extensions = [
    "sphinxcontrib.matlab",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "myst_parser",
    "sphinx_copybutton",
    "sphinx_design",
]

autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
}

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
    "tasklist",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

html_theme = "furo"
html_static_path = ["_static"]
html_title = f"{project}"

# Intersphinx -- no peers for a standalone submodule; peers added in labcode_main.
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}
