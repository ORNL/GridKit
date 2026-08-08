"""Resolve links between source READMEs to the pages that publish them.

A README in the source tree links to its neighbours with ordinary relative
paths, because it has to stay readable on GitHub. The documentation tree is
organised by topic rather than by source layout, so those paths point outside
the Sphinx source directory and resolve to nothing. Every included file knows
which page includes it, which is enough to redirect the link to that page.

The work happens ahead of MyST's own resolver, so anything left unresolved
still produces its usual warning.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

from docutils import nodes
from sphinx import addnodes
from sphinx.transforms.post_transforms import SphinxPostTransform
from sphinx.util import logging
from sphinx.util.nodes import make_refnode

if TYPE_CHECKING:
    from sphinx.application import Sphinx

logger = logging.getLogger(__name__)

_INCLUDE = re.compile(r"^```\{include\}\s+(\S+)", re.M)


@dataclass(frozen=True, slots=True)
class Pages:
    """The page that publishes each included source file."""

    publisher: dict[Path, str]

    def of(self, path: Path) -> str | None:
        return self.publisher.get(path)


def _scan(source: Path) -> Pages:
    publisher: dict[Path, str] = {}
    for page in sorted(source.rglob("*.md")):
        docname = page.relative_to(source).with_suffix("").as_posix()
        for reference in _INCLUDE.findall(page.read_text(encoding="utf-8")):
            path = (page.parent / reference).resolve()
            if path in publisher:
                logger.warning(
                    "%s is published by both %s and %s",
                    path, publisher[path], docname, type="gridkit",
                )
                continue
            publisher[path] = docname
    return Pages(publisher)


def pages(app: Sphinx) -> Pages:
    found = getattr(app, "_gridkit_pages", None)
    if found is None:
        found = _scan(Path(app.srcdir))
        app._gridkit_pages = found
    return found


class SourceLinks(SphinxPostTransform):
    """Point links between source READMEs at the pages that publish them."""

    default_priority = 5  # ahead of MystReferenceResolver, which runs at 9

    def run(self, **kwargs: Any) -> None:
        found = pages(self.app)
        for node in list(self.document.findall(addnodes.pending_xref)):
            if node.get("reftype") != "myst":
                continue

            target, anchor = self._target(node)
            if target is None:
                continue
            docname = found.of(target)
            if docname is None:
                continue

            content = node[0].deepcopy()
            node.replace_self(
                make_refnode(
                    self.app.builder,
                    node.get("refdoc") or self.env.docname,
                    docname,
                    anchor,
                    content,
                )
            )

    def _target(self, node: nodes.Element) -> tuple[Path | None, str | None]:
        reference = node.get("reftarget", "")
        if not reference:
            return None, None

        # A link MyST already read as a document carries an absolute stem.
        if node.get("refdomain") == "doc":
            return Path(reference + ".md"), node.get("reftargetid")

        path, _, anchor = reference.partition("#")
        if not path.endswith(".md"):
            return None, None
        # The node remembers the file its text came from, which is the file the
        # link is relative to even when several sources share one page.
        source = node.source or self.env.doc2path(self.env.docname)
        return (Path(source).parent / path).resolve(), anchor or None


def setup(app: Sphinx) -> None:
    app.add_post_transform(SourceLinks)
