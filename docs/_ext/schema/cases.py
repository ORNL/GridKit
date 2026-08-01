from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path

from docutils import nodes
from sphinx.util.docutils import SphinxDirective


class CaseModelsDirective(SphinxDirective):
    """Render a model-count table for a case file.

    The single argument is a repository-root-relative path to the case JSON,
    for MyST documents only:

        ```{case-models} examples/PhasorDynamics/Small/TwoArea/twoarea.json
        ```
    """

    required_arguments = 1

    def run(self) -> list[nodes.Node]:
        path = Path(self.env.srcdir).parent / self.arguments[0]
        self.env.note_dependency(path)
        case = json.loads(path.read_text(encoding="utf-8"))
        counts = Counter(
            entry["class"]
            for section in ("buses", "devices")
            for entry in case.get(section, ())
        )
        rows = "\n".join(
            f"| {{ref}}`{name} <model-phasor-dynamics-{_slug(name)}>` | {count} |"
            for name, count in counts.items()
        )
        return self.parse_text_to_nodes(f"| Model | Count |\n| --- | --- |\n{rows}")


def _slug(class_name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", class_name.lower())


SIZE_TIERS = {"Tiny", "Small", "Medium", "Large", "Huge"}


class CaseExamplesDirective(SphinxDirective):
    """Render a table of examples that run a case file.

    The single argument is a repository-root-relative path to the case JSON.
    An example uses the case when the case file sits in the example directory
    or when the example's CMake or solver configuration names it.
    """

    required_arguments = 1

    def run(self) -> list[nodes.Node]:
        repository = Path(self.env.srcdir).parent
        case = repository / self.arguments[0]
        rows = []
        for cmake in sorted((repository / "examples").rglob("CMakeLists.txt")):
            directory = cmake.parent
            configs = [cmake, *directory.glob("*.solver.json")]
            if case.parent == directory or any(
                case.name in config.read_text(encoding="utf-8") for config in configs
            ):
                parts = directory.relative_to(repository / "examples").parts
                domain = _kebab(parts[0])
                name = [part for part in parts[1:] if part not in SIZE_TIERS]
                slug = "-".join(_kebab(part) for part in name)
                title = " ".join(name)
                rows.append(f"| {{ref}}`{title} <example-{domain}-{slug}>` |")
        if not rows:
            return self.parse_text_to_nodes("No example currently runs this case.")
        return self.parse_text_to_nodes("| Example |\n| --- |\n" + "\n".join(rows))


def _kebab(name: str) -> str:
    return re.sub(r"(?<=[a-z0-9])(?=[A-Z])", "-", name).lower()
