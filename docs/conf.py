from pathlib import Path
from xml.etree import ElementTree

from exhale import utils as exhale_utils

project = "GridKit"
author = "GridKit Developers"

docs_dir = Path(__file__).parent.resolve()

for xml_path in (docs_dir / "xml").glob("*.xml"):
    tree = ElementTree.parse(xml_path)
    occurrences = {}
    changed = False
    for member in tree.iter("memberdef"):
        member_id = member.get("id")
        if member_id is None:
            continue
        occurrence = occurrences.get(member_id, 0)
        occurrences[member_id] = occurrence + 1
        if occurrence:
            member.set("id", f"{member_id}_overload_{occurrence}")
            changed = True
    if changed:
        tree.write(xml_path, encoding="utf-8", xml_declaration=True)

extensions = ["breathe", "exhale", "myst_parser", "sphinx_design"]

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
html_static_path = ["_static"]
html_css_files = ["css/gridkit.css"]
html_theme_options = {
    "collapse_navigation": False,
    "includehidden": True,
    "navigation_depth": 6,
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
