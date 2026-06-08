from pathlib import Path

project = "GridKit"
author = "GridKit Developers"

docs_dir = Path(__file__).parent.resolve()

extensions = ["breathe", "myst_parser", "sphinx.ext.graphviz"]

breathe_projects = {"GridKit": str(docs_dir / "xml")}
breathe_default_project = "GridKit"

html_theme = "sphinx_rtd_theme"
html_extra_path = ["Figures"]
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 6,
    "titles_only": True,
}

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

myst_enable_extensions = [
    "amsmath",
    "dollarmath",
    "html_image",
]
myst_heading_anchors = 5

exclude_patterns = ["_build", "xml", "README.md", "superpowers/**"]
