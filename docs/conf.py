from pathlib import Path

from exhale import utils as exhale_utils

project = "GridKit"
author = "GridKit Developers"

docs_dir = Path(__file__).parent.resolve()

extensions = ["breathe", "exhale", "myst_parser", "sphinx.ext.graphviz"]

breathe_projects = {"GridKit": str(docs_dir / "xml")}
breathe_default_project = "GridKit"


def exhale_specs(kind):
    if kind in {"class", "struct"}:
        return [":members:"]
    return []


exhale_args = {
    "containmentFolder": "./api/reference",
    "rootFileName": "EXCLUDE",
    "doxygenStripFromPath": str(docs_dir.parent),
    "createTreeView": False,
    "customSpecificationsMapping": exhale_utils.makeCustomSpecificationsMapping(
        exhale_specs
    ),
    "contentsDirectives": False,
    "fullToctreeMaxDepth": 1,
    "pageLevelConfigMeta": ":orphan:",
}

primary_domain = "cpp"

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

exclude_patterns = [
    "_build",
    "api/generated/**",
    "xml",
    "README.md",
    "superpowers/**",
]


def rewrite_included_root_readme_links(app, relative_path, parent_docname, source):
    if Path(str(relative_path)).name == "README.md":
        source[0] = source[0].replace("](INSTALL.md", "](generated/install.md")


def setup(app):
    app.connect("include-read", rewrite_included_root_readme_links)
