from pathlib import Path

from exhale import utils as exhale_utils

project = "GridKit"
author = "GridKit Developers"

docs_dir = Path(__file__).parent.resolve()

extensions = ["breathe", "exhale", "myst_parser"]

breathe_projects = {"GridKit": str(docs_dir / "xml")}
breathe_default_project = "GridKit"


# Exhale's default class/struct pages also include protected and undocumented members.
def public_member_specs(kind):
    if kind in {"class", "struct"}:
        return [":members:"]
    return []


exhale_args = {
    "containmentFolder": "./api/reference",
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
html_theme_options = {
    "collapse_navigation": True,
    "navigation_depth": 6,
    "titles_only": True,
}

myst_enable_extensions = [
    "amsmath",
    "dollarmath",
    "html_image",
]
myst_fence_as_directive = ["math"]
myst_heading_anchors = 5


exclude_patterns = [
    "_build",
    "README.md",
]
