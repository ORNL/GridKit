import os
import sys
from pathlib import Path

from exhale import utils as exhale_utils

project = "GridKit"
author = "GridKit Developers"

docs_dir = Path(__file__).parent.resolve()
sys.path.insert(0, str(docs_dir / "_ext"))

extensions = [
    "breathe",
    "exhale",
    "gridkit",
    "myst_parser",
    "sphinx-jsonschema",
    "sphinx_design",
]

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
html_extra_path = ["../GridKit/Model/case.schema.json"]
html_static_path = ["_static"]
html_css_files = ["css/gridkit.css"]
html_theme_options = {
    "collapse_navigation": True,
    "navigation_depth": 2,
    "titles_only": True,
}

myst_enable_extensions = [
    "alert",
    "amsmath",
    "colon_fence",
    "dollarmath",
    "html_image",
    "substitution",
]

# The published schema URL follows the version RTD is building. The fallback
# matches the schema version slug configured on RTD for local builds.
schema_url = (
    os.environ.get(
        "READTHEDOCS_CANONICAL_URL", "https://gridkit.readthedocs.io/en/schema/"
    ).rstrip("/")
    + "/case.schema.json"
)
myst_substitutions = {
    "schema_url": f"<{schema_url}>",
    "schema_example": f'```json\n{{\n  "$schema": "{schema_url}"\n}}\n```',
}
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
