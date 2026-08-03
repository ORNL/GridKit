"""Directives that render views of the repository.

Each one is a query plus a shape. Tables and galleries are emitted as MyST and
handed back to the parser, so cross references resolve the ordinary way and an
unresolved label is a build warning rather than a silent gap.
"""

from __future__ import annotations

import os
from collections.abc import Iterable, Sequence
from pathlib import Path

import yaml
from docutils import nodes
from sphinx.util.docutils import SphinxDirective

from .doxygen import Model
from .repository import Case, Repository, RepositoryError, repository

GRID = "::::{{grid}} {columns}\n:gutter: 3\n\n{cards}\n::::"

CARD = """\
:::{{grid-item-card}} {title}
:img-top: {image}
:link: {label}
:link-type: ref
:class-card: gk-card

{description}

+++
{footer}
:::
"""


def _ref(label: str) -> str:
    """A link that takes its text from the page it points at.

    Sphinx supplies the title, so a table cannot disagree with the page it
    links to, and an undefined label becomes a build warning.
    """
    return f"{{ref}}`{label}`"


def _class_ref(model: Model) -> str:
    """A link labeled with the class string used in a case file.

    Model pages are not one per model. Bus and BusInfinite share a page, so a
    link that borrows the page title cannot distinguish them.
    """
    return f"[`{model.name}`](#{model.label})"


def _example_refs(repository: Repository, case: Case) -> str:
    return ", ".join(
        _ref(example.label) for example in repository.examples_of(case)
    )


class GridKitDirective(SphinxDirective):
    """Base class for directives that read the repository."""

    empty = "Nothing to show."

    @property
    def repository(self) -> Repository:
        found = repository(self.env.app)
        for source in found.sources:
            self.env.note_dependency(str(source))
        return found

    def model(self, name: str):
        try:
            return self.repository.model(name)
        except RepositoryError as error:
            raise self.error(str(error)) from error

    def case(self, relative: str) -> Case:
        try:
            return self.repository.case(relative)
        except (RepositoryError, OSError) as error:
            raise self.error(str(error)) from error

    def table(
        self,
        headers: Sequence[str],
        rows: Iterable[Sequence[object]],
        align: Sequence[str] | None = None,
    ) -> list[nodes.Node]:
        """Render rows as a MyST table, omitting columns with no content.

        `align` gives one entry per header: "right" or "center" sets that
        column's alignment, anything else keeps the default.
        """
        cells = [[str(value).replace("|", r"\|") for value in row] for row in rows]
        if not cells:
            return self.parse_text_to_nodes(self.empty)

        keep = [index for index in range(len(headers)) if any(row[index] for row in cells)]
        head = [headers[index] for index in keep]
        marks = {"right": "---:", "center": ":---:"}
        rules = [marks.get(align[index] if align else "", "---") for index in keep]
        lines = [
            f"| {' | '.join(head)} |",
            f"| {' | '.join(rules)} |",
            *(f"| {' | '.join(row[index] for index in keep)} |" for row in cells),
        ]
        return self.parse_text_to_nodes(
            ":::{div} gk-table\n" + "\n".join(lines) + "\n:::"
        )

    def gallery(self, cards: Sequence[str]) -> list[nodes.Node]:
        if not cards:
            return self.parse_text_to_nodes(self.empty)
        columns = self.options.get("columns", "1 2 3 3")
        return self.parse_text_to_nodes(
            GRID.format(columns=columns, cards="\n".join(cards))
        )

    def relative(self, target: Path) -> str:
        """The path to `target` relative to the document being rendered."""
        source = Path(self.env.doc2path(self.env.docname)).parent
        return os.path.relpath(target, source).replace(os.sep, "/")


class ModelParameters(GridKitDirective):
    """Parameter table for one model, including its default value."""

    required_arguments = 1
    empty = "This model has no parameters."

    def run(self) -> list[nodes.Node]:
        model = self.model(self.arguments[0])
        return self.table(
            ("Parameter", "Symbol", "Units", "Description", "Default"),
            (
                (
                    f"`{item.name}`",
                    item.symbol,
                    item.unit,
                    item.description,
                    item.default,
                )
                for item in model.parameters
            ),
        )


class ModelPorts(GridKitDirective):
    """Buses and signals one model connects to."""

    required_arguments = 1
    empty = "This model has no ports."

    def run(self) -> list[nodes.Node]:
        model = self.model(self.arguments[0])
        groups = (
            ("Bus", model.buses),
            ("Signal in", model.signal_inputs),
            ("Signal out", model.signal_outputs),
        )
        return self.table(
            ("Port", "Kind", "Description"),
            (
                (f"`{item.name}`", kind, item.description)
                for kind, items in groups
                for item in items
            ),
        )


class ModelMonitors(GridKitDirective):
    """Monitorable variables of one model."""

    required_arguments = 1
    empty = "This model has no monitorable variables."

    def run(self) -> list[nodes.Node]:
        model = self.model(self.arguments[0])
        return self.table(
            ("Variable", "Symbol", "Description"),
            (
                (f"`{item.name}`", item.symbol, item.description)
                for item in model.monitors
            ),
        )


class ModelCases(GridKitDirective):
    """Cases containing one model, and the examples that run them."""

    required_arguments = 1
    empty = "No case uses this model."

    def run(self) -> list[nodes.Node]:
        found = self.repository
        model = self.model(self.arguments[0])
        return self.table(
            ("Case", "Instances", "Examples"),
            (
                (
                    _ref(case.label),
                    case.counts[model.name],
                    _example_refs(found, case),
                )
                for case in found.uses(model)
            ),
        )


class ModelCatalog(GridKitDirective):
    """Every model class a case file may declare."""

    empty = "No models found in the Doxygen output."

    def run(self) -> list[nodes.Node]:
        found = self.repository
        return self.table(
            ("Model", "Family", "Parameters", "Ports", "Monitors", "Cases", "Instances"),
            (
                (
                    _class_ref(model),
                    model.family,
                    len(model.parameters),
                    len(model.ports),
                    len(model.monitors),
                    len(found.uses(model)) or "",
                    found.instances(model) or "",
                )
                for model in found.models.values()
            ),
        )


class CaseModels(GridKitDirective):
    """Model counts for one case.

    The argument is the case file, relative to the repository root.
    """

    required_arguments = 1
    empty = "This case declares no buses or devices."

    def run(self) -> list[nodes.Node]:
        found = self.repository
        case = self.case(self.arguments[0])
        rows = []
        for name, count in case.counts.most_common():
            model = found.models.get(name)
            rows.append(
                (
                    _class_ref(model) if model else f"`{name}`",
                    model.family if model else "",
                    count,
                )
            )
        return self.table(("Model", "Family", "Count"), rows)


class CaseExamples(GridKitDirective):
    """Examples that run one case."""

    required_arguments = 1
    empty = "No example currently runs this case."

    def run(self) -> list[nodes.Node]:
        found = self.repository
        case = self.case(self.arguments[0])
        return self.table(
            ("Example", "Simulated time", "Events"),
            (
                (
                    _ref(example.label),
                    f"{example.duration:g} s" if example.duration else "",
                    example.events or "",
                )
                for example in found.examples_of(case)
            ),
        )


class CaseCatalog(GridKitDirective):
    """Every case with a documentation page."""

    empty = "No documented cases."

    def run(self) -> list[nodes.Node]:
        found = self.repository
        return self.table(
            ("Case", "Buses", "Devices", "Model types", "Examples"),
            (
                (
                    _ref(case.label),
                    case.buses,
                    case.devices,
                    len(case.counts),
                    _example_refs(found, case),
                )
                for case in found.documented
            ),
        )


class CatalogDirective(GridKitDirective):
    """Base class for directives that read an overhead-line type catalog.

    The first argument is the catalog file, relative to the repository root.
    """

    required_arguments = 1

    def catalog(self) -> dict:
        path = Path(self.env.app.srcdir).parent / self.arguments[0]
        self.env.note_dependency(str(path))
        try:
            found = yaml.safe_load(path.read_text(encoding="utf-8"))
        except (OSError, yaml.YAMLError) as error:
            raise self.error(
                f"cannot read catalog {self.arguments[0]}: {error}"
            ) from error
        if not isinstance(found, dict):
            raise self.error(f"{self.arguments[0]} is not a catalog document")
        return found

    @staticmethod
    def number(value: object) -> str:
        """Format a quantity, in math notation once `:g` turns exponential."""
        if not isinstance(value, (int, float)):
            return ""
        mantissa, _, exponent = f"{value:g}".partition("e")
        if not exponent:
            return mantissa
        return rf"${mantissa} \times 10^{{{int(exponent)}}}$"


class CatalogConductors(CatalogDirective):
    """Conductor types of one catalog file."""

    empty = "This catalog defines no conductor types."

    def run(self) -> list[nodes.Node]:
        conductors = self.catalog().get("conductors") or {}
        return self.table(
            (
                "Type", "Description", "Radius (outer / inner) [m]",
                "Conductivity [S/m]", "Permeability [-]", "Weight [N/m]",
            ),
            (
                (
                    f"`{name}`",
                    item.get("description", ""),
                    self._radius(item),
                    self.number(item.get("conductivity")),
                    self._permeability(item),
                    self.number(item.get("weight")),
                )
                for name, item in conductors.items()
            ),
            align=("", "", "right", "right", "right", "right"),
        )

    @classmethod
    def _radius(cls, item: dict) -> str:
        """One cell for both radii; solid types show the outer alone."""
        radius = item.get("radius") or {}
        outer = cls.number(radius.get("outer"))
        inner = radius.get("inner")
        return f"{outer} / {cls.number(inner)}" if inner else outer

    @classmethod
    def _permeability(cls, item: dict) -> str:
        """Blank at the default, so an all-default column drops out."""
        value = item.get("permeability", 1)
        return "" if value == 1 else cls.number(value)


class CatalogTowers(CatalogDirective):
    """Tower types of one catalog file."""

    empty = "This catalog defines no tower types."

    def run(self) -> list[nodes.Node]:
        towers = self.catalog().get("towers") or {}
        return self.table(
            ("Type", "Description", "Attachment points"),
            (
                (
                    f"`{name}`",
                    item.get("description", ""),
                    len(item.get("attachments") or {}),
                )
                for name, item in towers.items()
            ),
            align=("", "", "right"),
        )


class CatalogTower(CatalogDirective):
    """Attachment points of one tower type from one catalog file."""

    required_arguments = 2
    empty = "This tower type defines no attachment points."

    def run(self) -> list[nodes.Node]:
        towers = self.catalog().get("towers") or {}
        tower = towers.get(self.arguments[1])
        if tower is None:
            raise self.error(
                f"no tower type {self.arguments[1]} in {self.arguments[0]}"
            )
        return self.table(
            ("Attachment", "x [m]", "h [m]"),
            (
                (f"`{name}`", self.number(point.get("x")),
                 self.number(point.get("h")))
                for name, point in (tower.get("attachments") or {}).items()
            ),
            align=("", "right", "right"),
        )


class CaseGallery(GridKitDirective):
    """The same cases as cards, each led by its one-line diagram."""

    empty = "No documented cases."
    option_spec = {"columns": lambda argument: argument or "1 2 3 3"}

    def run(self) -> list[nodes.Node]:
        return self.gallery(
            [
                CARD.format(
                    title=case.title,
                    image=self.relative(case.oneline),
                    label=case.label,
                    description=case.description,
                    footer=f"{case.buses} buses, {case.devices} devices",
                )
                for case in self.repository.documented
                if case.oneline
            ]
        )


DIRECTIVES = {
    "model-parameters": ModelParameters,
    "model-ports": ModelPorts,
    "model-monitors": ModelMonitors,
    "model-cases": ModelCases,
    "model-catalog": ModelCatalog,
    "case-models": CaseModels,
    "case-examples": CaseExamples,
    "case-catalog": CaseCatalog,
    "case-gallery": CaseGallery,
    "catalog-conductors": CatalogConductors,
    "catalog-towers": CatalogTowers,
    "catalog-tower": CatalogTower,
}
