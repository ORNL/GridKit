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

    domain: str
    slug: str
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

    for entry in (*data.get("buses", ()), *data.get("devices", ())):
        counts[entry["class"]] += 1

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
        page=page if page.is_file() else None,
    )


def _read_examples(root: Path, cases: dict[str, Case]) -> Iterator[Example]:
    for cmake in sorted(root.glob(EXAMPLE_FILES)):
        directory = cmake.parent
        build = cmake.read_text(encoding="utf-8")
        # A directory that only gathers subdirectories is not an example.
        if not _TEST.search(build):
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
            domain=kebab(relative.parts[0]),
            slug="-".join(kebab(part) for part in name),
            case=case.path if case else None,
            duration=study.get("tmax"),
            events=len(study.get("events", ())),
        )


class Repository:
    def __init__(self, root: Path, xml: Path, pages: Path):
        self.root = root
        self.xml = xml
        self.pages = pages

    @cached_property
    def models(self) -> dict[str, Model]:
        return read_models(self.xml)

    @cached_property
    def cases(self) -> tuple[Case, ...]:
        return tuple(
            _read_case(self.root, path, self.pages)
            for path in sorted(self.root.glob(CASE_FILES))
        )

    @cached_property
    def examples(self) -> tuple[Example, ...]:
        by_name = {case.path.name: case for case in self.cases}
        return tuple(_read_examples(self.root, by_name))

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
        return self._uses.get(model.name, ())

    def instances(self, model: Model) -> int:
        return sum(case.counts[model.name] for case in self.uses(model))

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
