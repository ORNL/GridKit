from pathlib import Path
import sys

from exhale import utils as exhale_utils

project = "GridKit"
author = "GridKit Developers"

docs_dir = Path(__file__).parent.resolve()
sys.path.insert(0, str(docs_dir / "_ext"))

extensions = ["breathe", "exhale", "myst_parser", "sphinx_design", "gridkit_links"]

breathe_projects = {"GridKit": str(docs_dir / "xml")}
breathe_default_project = "GridKit"


# Exhale's default class/struct pages also include protected and undocumented members.
def public_member_specs(kind):
    if kind in {"class", "struct"}:
        return [":members:"]
    return []


exhale_args = {
    "containmentFolder": "./reference/api/generated",
    "rootFileName": "EXCLUDE",
    "doxygenStripFromPath": str(docs_dir.parent),
    "customSpecificationsMapping": exhale_utils.makeCustomSpecificationsMapping(
        public_member_specs
    ),
    "contentsDirectives": False,
    "pageLevelConfigMeta": ":orphan:",
}

primary_domain = "cpp"

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_css_files = ["css/gridkit.css"]
html_theme_options = {
    "collapse_navigation": False,
    "includehidden": True,
    "navigation_depth": 4,
    "titles_only": True,
}

myst_enable_extensions = [
    "alert",
    "amsmath",
    "colon_fence",
    "dollarmath",
    "html_image",
]
myst_fence_as_directive = ["math"]
myst_heading_anchors = 5


exclude_patterns = [
    "_build",
    "README.md",
]

# Breathe renders public nested types with their parent and Exhale also gives
# those types standalone pages. Sphinx otherwise reports the intentional
# duplicate declarations when both pages are read.
suppress_warnings = ["duplicate_declaration.cpp"]
