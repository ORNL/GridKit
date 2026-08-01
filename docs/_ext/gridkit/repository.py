"""Models, cases, and examples, as the documentation sees them.

Every directive renders a view over one `Repository`, built once per Sphinx
process, so each fact is parsed once and derived in one place. Nothing here
touches the build environment. Cases and examples are discovered from the tree
rather than declared by the pages that present them.
"""

from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from collections.abc import Iterator
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING

from .doxygen import Model, read_models

if TYPE_CHECKING:
    from sphinx.application import Sphinx

# Where things live. These are the only paths to change if the tree moves.
EXAMPLES = "examples"
CASE_FILES = f"{EXAMPLES}/**/*.case.json"
EXAMPLE_FILES = f"{EXAMPLES}/**/CMakeLists.txt"
SOLVER_FILES = f"{EXAMPLES}/**/*.solver.json"
CASE_PAGES = "cases"

SIZE_TIERS = {"Tiny", "Small", "Medium", "Large", "Huge"}

_KEBAB = re.compile(r"(?<=[a-z0-9])(?=[A-Z])")
_TEST = re.compile(r"add_test\s*\(\s*NAME\s+(\S+)")


class RepositoryError(RuntimeError):
    pass


def kebab(name: str) -> str:
    return _KEBAB.sub("-", name).replace("_", "-").lower()


@dataclass(frozen=True, slots=True)
class Range:
    """The interval a parameter spans across the case corpus."""

    low: float
    high: float
    count: int

    def __or__(self, other: Range) -> Range:
        return Range(
            min(self.low, other.low),
            max(self.high, other.high),
            self.count + other.count,
        )

    def __str__(self) -> str:
        span = f"{self.low:g}" if self.low == self.high else f"{self.low:g} to {self.high:g}"
        return f"{span} (n={self.count})"


@dataclass(frozen=True, slots=True)
class Case:
    """A case file, its summary attributes, and its contents."""

    path: Path
    domain: str
    slug: str
    title: str
    description: str
    buses: int
    devices: int
    counts: Counter[str]
    ranges: dict[str, dict[str, Range]]
    page: Path | None

    @property
    def label(self) -> str:
        return f"case-{self.domain}-{self.slug}"

    @property
    def oneline(self) -> Path | None:
        return next(iter(sorted(self.path.parent.glob("oneline.*"))), None)


@dataclass(frozen=True, slots=True)
class Example:
    """A registered test and the study it executes."""

    directory: Path
    domain: str
    slug: str
    title: str
    tests: tuple[str, ...]
    case: Path | None
    duration: float | None
    events: int

    @property
    def label(self) -> str:
        return f"example-{self.domain}-{self.slug}"


def _heading(readme: Path) -> str:
    """The title of a README, which is also the title of the page including it."""
    if not readme.is_file():
        return ""
    for line in readme.read_text(encoding="utf-8").splitlines():
        if line.strip().startswith("# "):
            return line.strip().removeprefix("# ").strip().strip("*")
    return ""


def _read_case(root: Path, path: Path, pages: Path) -> Case:
    data = json.loads(path.read_text(encoding="utf-8"))
    header = data.get("header", {})
    counts: Counter[str] = Counter()
    ranges: dict[str, dict[str, Range]] = defaultdict(dict)

    for entry in (*data.get("buses", ()), *data.get("devices", ())):
        name = entry["class"]
        counts[name] += 1
        for parameter, value in entry.get("params", {}).items():
            # A parameter may also be `true`, meaning the model derives it.
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                continue
            span = Range(value, value, 1)
            current = ranges[name].get(parameter)
            ranges[name][parameter] = span if current is None else current | span

    relative = path.relative_to(root / EXAMPLES)
    domain = kebab(relative.parts[0])
    slug = kebab(path.parent.name)
    page = pages / domain / f"{slug}.md"
    return Case(
        path=path,
        domain=domain,
        slug=slug,
        title=_heading(path.parent / "README.md") or header.get("case_name") or path.parent.name,
        description=header.get("case_description", ""),
        buses=len(data.get("buses", ())),
        devices=len(data.get("devices", ())),
        counts=counts,
        ranges=dict(ranges),
        page=page if page.is_file() else None,
    )


def _read_examples(root: Path, cases: dict[str, Case]) -> Iterator[Example]:
    for cmake in sorted(root.glob(EXAMPLE_FILES)):
        directory = cmake.parent
        build = cmake.read_text(encoding="utf-8")
        tests = tuple(_TEST.findall(build))
        # A directory that only gathers subdirectories is not an example.
        if not tests:
            continue

        solvers = sorted(directory.glob("*.solver.json"))
        study = json.loads(solvers[0].read_text(encoding="utf-8")) if solvers else {}
        named = f"{build}\n{json.dumps(study)}"
        case = next(
            (
                candidate
                for name, candidate in cases.items()
                if candidate.path.parent == directory or name in named
            ),
            None,
        )

        relative = directory.relative_to(root / EXAMPLES)
        name = [part for part in relative.parts[1:] if part not in SIZE_TIERS]
        yield Example(
            directory=directory,
            domain=kebab(relative.parts[0]),
            slug="-".join(kebab(part) for part in name),
            title=" ".join(name),
            tests=tests,
            case=case.path if case else None,
            duration=study.get("tmax"),
            events=len(study.get("events", ())),
        )


class Repository:
    def __init__(self, root: Path, xml: Path, pages: Path):
        self.root = root
        self.models = read_models(xml)
        self.cases = tuple(
            _read_case(root, path, pages) for path in sorted(root.glob(CASE_FILES))
        )
        by_name = {case.path.name: case for case in self.cases}
        self.examples = tuple(_read_examples(root, by_name))

    @cached_property
    def sources(self) -> tuple[Path, ...]:
        """Files the tables are derived from, for incremental rebuilds.

        The Doxygen XML is deliberately excluded. Regenerating it also
        regenerates the API pages, which rebuilds the site regardless.
        """
        return (
            *(case.path for case in self.cases),
            *sorted(self.root.glob(EXAMPLE_FILES)),
            *sorted(self.root.glob(SOLVER_FILES)),
        )

    @cached_property
    def by_json_name(self) -> dict[str, Model]:
        return {model.json_name: model for model in self.models.values()}

    @cached_property
    def documented(self) -> tuple[Case, ...]:
        """Cases with a page under `docs/cases/`, largest first."""
        return tuple(
            sorted(
                (case for case in self.cases if case.page is not None),
                key=lambda case: -case.buses,
            )
        )

    @cached_property
    def _uses(self) -> dict[str, tuple[Case, ...]]:
        uses: dict[str, list[Case]] = defaultdict(list)
        for case in self.documented:
            for name in case.counts:
                uses[name].append(case)
        return {
            name: tuple(sorted(found, key=lambda case: -case.counts[name]))
            for name, found in uses.items()
        }

    @cached_property
    def _ranges(self) -> dict[str, dict[str, Range]]:
        merged: dict[str, dict[str, Range]] = defaultdict(dict)
        for case in self.documented:
            for name, parameters in case.ranges.items():
                for parameter, span in parameters.items():
                    current = merged[name].get(parameter)
                    merged[name][parameter] = span if current is None else current | span
        return merged

    def model(self, name: str) -> Model:
        try:
            return self.models[name]
        except KeyError:
            known = ", ".join(self.models)
            raise RepositoryError(f"unknown model {name!r}; known models: {known}") from None

    def case(self, relative: str) -> Case:
        path = self.root / relative
        for case in self.cases:
            if case.path == path:
                return case
        raise RepositoryError(f"{relative} is not a case file")

    def uses(self, model: Model) -> tuple[Case, ...]:
        """Documented cases containing this model.

        Case files also sit beside the examples that validate a single model.
        Those are fixtures, not systems anyone studies, so they are excluded
        here by the same rule the catalog and gallery use: a case counts once
        it has a page.
        """
        return self._uses.get(model.json_name, ())

    def instances(self, model: Model) -> int:
        return sum(case.counts[model.json_name] for case in self.uses(model))

    def ranges(self, model: Model) -> dict[str, Range]:
        return self._ranges.get(model.json_name, {})

    def examples_of(self, case: Case) -> tuple[Example, ...]:
        return tuple(example for example in self.examples if example.case == case.path)


def repository(app: Sphinx) -> Repository:
    """The shared repository, constructed on first use."""
    repo = getattr(app, "_gridkit_repository", None)
    if repo is None:
        source = Path(app.srcdir)
        repo = Repository(source.parent, source / "xml", source / CASE_PAGES)
        app._gridkit_repository = repo
    return repo
