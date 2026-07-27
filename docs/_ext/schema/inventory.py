from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Literal


NAMESPACE = "GridKit::PhasorDynamics"
BUS_BASE = f"{NAMESPACE}::BusBase"
ALIASES = {
    "Parameters": "parameters",
    "Buses": "buses",
    "SignalInputs": "signal_inputs",
    "SignalOutputs": "signal_outputs",
    "MonitorableVariables": "monitors",
}


class InventoryError(RuntimeError):
    pass


@dataclass(frozen=True, slots=True)
class Item:
    name: str
    description: str = ""


@dataclass(frozen=True, slots=True)
class Model:
    name: str
    kind: Literal["bus", "device"]
    parameters: tuple[Item, ...] = ()
    buses: tuple[Item, ...] = ()
    signal_inputs: tuple[Item, ...] = ()
    signal_outputs: tuple[Item, ...] = ()
    monitors: tuple[Item, ...] = ()


def _text(node: ET.Element | None) -> str:
    if node is None:
        return ""
    return re.sub(r"\s+", " ", "".join(node.itertext())).strip()


def _item(node: ET.Element) -> Item:
    brief = node.find("briefdescription")
    description = _text(brief)
    formula = brief.find(".//formula") if brief is not None else None
    if formula is not None:
        description = description.removeprefix(_text(formula)).strip(" -:.")
    return Item(_text(node.find("name")), description)


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


def _items(index: _Index, owner: ET.Element, alias: str) -> tuple[Item, ...]:
    enum_refid = _alias(owner, alias, "member")
    if enum_refid is None:
        return ()

    enum = index.member(enum_refid)
    if enum.get("kind") != "enum":
        raise InventoryError(f"{alias} does not name a documented enum")

    return tuple(
        _item(value)
        for value in enum.findall("enumvalue")
        if _text(value.find("name")) not in {"", "NONE", "SIZE"}
    )


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

        sections = {
            field: _items(index, data, alias)
            for alias, field in ALIASES.items()
        }
        models[name] = Model(
            name=name,
            kind="bus" if _inherits(index, refid, BUS_BASE) else "device",
            **sections,
        )

    if not models:
        raise InventoryError(f"no models found in {directory}")
    return dict(sorted(models.items()))
