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
            f"| [{name}](#model-phasor-dynamics-{_slug(name)}) | {count} |"
            for name, count in counts.items()
        )
        return self.parse_text_to_nodes(f"| Model | Count |\n| --- | --- |\n{rows}")


def _slug(class_name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", class_name.lower())
