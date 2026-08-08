"""GridKit model definitions, read from Doxygen XML.

A model declares its parameters, ports, and monitorable variables as documented
enums, reached through its `ModelDataT` alias. Reading them here keeps the JSON
schema and the documentation tables from diverging.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Literal


DOMAIN = "phasor-dynamics"
NAMESPACE = "GridKit::PhasorDynamics"
BUS_BASE = f"{NAMESPACE}::BusBase"
SOURCE_ROOT = Path("GridKit/Model/PhasorDynamics")

ALIASES = {
    "Parameters": "parameters",
    "Buses": "buses",
    "SignalInputs": "signal_inputs",
    "SignalOutputs": "signal_outputs",
    "MonitorableVariables": "monitors",
}

_UNIT = re.compile(r"^\[([^]]*)\]\s*")
_ALNUM = re.compile(r"[^a-z0-9]")


class InventoryError(RuntimeError):
    pass


@dataclass(frozen=True, slots=True)
class Item:
    """One documented enum value: a parameter, a port, or a monitor."""

    name: str
    symbol: str = ""
    unit: str = ""
    description: str = ""
    default: str = ""


@dataclass(frozen=True, slots=True)
class Model:
    name: str
    kind: Literal["bus", "device"]
    family: str = ""
    parameters: tuple[Item, ...] = ()
    buses: tuple[Item, ...] = ()
    signal_inputs: tuple[Item, ...] = ()
    signal_outputs: tuple[Item, ...] = ()
    monitors: tuple[Item, ...] = ()

    @property
    def slug(self) -> str:
        return _ALNUM.sub("", self.name.lower())

    @property
    def label(self) -> str:
        return f"model-{DOMAIN}-{self.slug}"

    @property
    def ports(self) -> tuple[Item, ...]:
        return (*self.buses, *self.signal_inputs, *self.signal_outputs)


def _text(node: ET.Element | None) -> str:
    if node is None:
        return ""
    return re.sub(r"\s+", " ", "".join(node.itertext())).strip()


def _number(text: str) -> str:
    try:
        return f"{float(text):g}"
    except ValueError:
        return text


def _item(node: ET.Element, defaults: dict[str, str]) -> Item:
    brief = node.find("briefdescription")
    description = _text(brief)

    # A leading formula is the quantity's symbol, and a bracketed word after it
    # is its unit: `///< \f$H\f$ [s] Rotor inertia`.
    symbol = _text(brief.find(".//formula")) if brief is not None else ""
    if symbol:
        description = description.removeprefix(symbol).strip(" -:.")
    unit = ""
    if match := _UNIT.match(description):
        unit, description = match.group(1), description[match.end() :]

    name = _text(node.find("name"))
    return Item(name, symbol, unit, description, defaults.get(f"{name}_", ""))


def _members(compound: ET.Element, kind: str):
    return compound.findall(f"./sectiondef/memberdef[@kind='{kind}']")


class _Index:
    def __init__(self, directory: Path):
        self.compounds: dict[str, ET.Element] = {}
        self.members: dict[str, ET.Element] = {}

        try:
            entries = ET.parse(directory / "index.xml").getroot()
            for entry in entries.findall("compound"):
                refid = entry.attrib["refid"]
                root = ET.parse(directory / f"{refid}.xml").getroot()
                compound = root.find("compounddef")
                if compound is None:
                    raise InventoryError(f"{refid}.xml has no compounddef")
                self.compounds[refid] = compound
                for member in compound.findall("./sectiondef/memberdef"):
                    if member_id := member.get("id"):
                        self.members[member_id] = member
        except (OSError, ET.ParseError, KeyError) as error:
            raise InventoryError(f"cannot read Doxygen XML: {error}") from error

    def compound(self, refid: str) -> ET.Element:
        try:
            return self.compounds[refid]
        except KeyError as error:
            raise InventoryError(f"unknown Doxygen compound {refid}") from error

    def member(self, refid: str) -> ET.Element:
        try:
            return self.members[refid]
        except KeyError as error:
            raise InventoryError(f"unknown Doxygen member {refid}") from error


def _alias(compound: ET.Element, name: str, kind: str) -> str | None:
    for member in _members(compound, "typedef"):
        if _text(member.find("name")) != name:
            continue
        for ref in member.findall("./type/ref"):
            if ref.get("kindref") == kind:
                return ref.get("refid")
        raise InventoryError(f"{name} alias has no {kind} reference")
    return None


def _inherits(index: _Index, refid: str, base_name: str) -> bool:
    pending = [refid]
    visited = set()
    while pending:
        current = pending.pop()
        if current in visited:
            continue
        visited.add(current)
        compound = index.compound(current)
        if _text(compound.find("compoundname")) == base_name:
            return True
        pending.extend(
            base.attrib["refid"]
            for base in compound.findall("basecompoundref")
            if "refid" in base.attrib
        )
    return False


def _items(
    index: _Index,
    owner: ET.Element,
    alias: str,
    defaults: dict[str, str],
) -> tuple[Item, ...]:
    enum_refid = _alias(owner, alias, "member")
    if enum_refid is None:
        return ()

    enum = index.member(enum_refid)
    if enum.get("kind") != "enum":
        raise InventoryError(f"{alias} does not name a documented enum")

    return tuple(
        _item(value, defaults)
        for value in enum.findall("enumvalue")
        if _text(value.find("name")) not in {"", "NONE", "SIZE"}
    )


def _defaults(compound: ET.Element) -> dict[str, str]:
    """Member initializers, keyed by member name."""
    values = {}
    for member in _members(compound, "variable"):
        initializer = _text(member.find("initializer")).strip("{}= ")
        if initializer:
            values[_text(member.find("name"))] = _number(initializer)
    return values


def _family(compound: ET.Element) -> str:
    location = compound.find("location")
    directory = Path(location.get("file", "")).parent if location is not None else Path()
    try:
        parts = directory.relative_to(SOURCE_ROOT).parts
    except ValueError:
        return ""
    # A model nested one level deeper than the domain root belongs to a family.
    return parts[0] if len(parts) > 1 else ""


def read_models(directory: Path) -> dict[str, Model]:
    index = _Index(directory)
    models: dict[str, Model] = {}

    for refid, compound in index.compounds.items():
        cpp_name = _text(compound.find("compoundname"))
        if not cpp_name.startswith(f"{NAMESPACE}::"):
            continue

        data_refid = _alias(compound, "ModelDataT", "compound")
        if data_refid is None:
            continue
        data = index.compound(data_refid)
        name = cpp_name.rsplit("::", 1)[-1]
        if name in models:
            raise InventoryError(f"duplicate model class {name}")

        defaults = _defaults(compound)
        sections = {
            field: _items(index, data, alias, defaults if field == "parameters" else {})
            for alias, field in ALIASES.items()
        }
        models[name] = Model(
            name=name,
            kind="bus" if _inherits(index, refid, BUS_BASE) else "device",
            family=_family(compound),
            **sections,
        )

    if not models:
        raise InventoryError(f"no models found in {directory}")
    return dict(sorted(models.items()))
