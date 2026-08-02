"""The JSON Schema for overhead-line inputs and their type catalogs.

The schema is authored here and generated at documentation build time; the
repository does not commit the JSON artifact. Read the Docs copies the
generated file to a stable URL for editors, validators, and CI, and the
validation tests import :func:`build_schema` directly so the two never
drift.

One schema describes both document kinds. The root is the line document.
``$defs/catalog`` is the catalog of named conductor and tower types as it
appears embedded in a line document, and ``$defs/catalog_document`` is the
standalone form written as a YAML catalog file.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

FALLBACK_ID = "https://gridkit.readthedocs.io/en/schema/line.schema.json"

#: Catalog names and attachment names are lower-case kebab-case.
NAME = "^[a-z0-9][a-z0-9-]*$"

_PROVENANCE = (
    "Provenance: drawing number, data sheet, standard, publication, or an "
    "explicit statement that the data was assumed."
)


def _named_section(ref: str, description: str) -> dict[str, Any]:
    return {
        "type": "object",
        "minProperties": 1,
        "propertyNames": {"pattern": NAME},
        "additionalProperties": {"$ref": ref},
        "description": description,
    }


def _catalog_sections() -> dict[str, Any]:
    return {
        "conductors": _named_section(
            "#/$defs/conductor_type",
            "Conductor types by name. Conductor entries reference these "
            "names.",
        ),
        "towers": _named_section(
            "#/$defs/tower_type",
            "Tower types by name. The line document's tower field "
            "references one of these names.",
        ),
    }


def _conductor_type() -> dict[str, Any]:
    return {
        "title": "Conductor type",
        "type": "object",
        "additionalProperties": False,
        "required": ["radius", "conductivity", "permeability", "weight"],
        "description": (
            "Conductor dimensions, material, and weight. Phase, circuit, "
            "tension, and the attachment point are line data and belong to "
            "the conductor entries of the line document."
        ),
        "properties": {
            "description": {"type": "string"},
            "source": {"type": "string", "description": _PROVENANCE},
            "radius": {
                "type": "number",
                "exclusiveMinimum": 0,
                "description": "Outer radius [m].",
            },
            "inner_radius": {
                "type": "number",
                "minimum": 0,
                "default": 0,
                "description": (
                    "Inner radius [m]. Zero for solid conductors; positive "
                    "for tubular conduction regions such as the aluminum "
                    "annulus of ACSR."
                ),
            },
            "conductivity": {
                "type": "number",
                "exclusiveMinimum": 0,
                "description": "Conductor conductivity [S/m].",
            },
            "permeability": {
                "type": "number",
                "exclusiveMinimum": 0,
                "description": "Conductor permeability [H/m].",
            },
            "weight": {
                "type": "number",
                "exclusiveMinimum": 0,
                "description": (
                    "Conductor weight per unit length [N/m]. Used with "
                    "tension to compute catenary sag."
                ),
            },
        },
    }


def _attachment() -> dict[str, Any]:
    return {
        "title": "Attachment point",
        "type": "object",
        "additionalProperties": False,
        "required": ["x", "height"],
        "properties": {
            "x": {
                "type": "number",
                "description": (
                    "Horizontal offset from the tower centreline [m]. "
                    "Positive to the right when looking along the route."
                ),
            },
            "height": {
                "type": "number",
                "exclusiveMinimum": 0,
                "description": (
                    "Attachment height above earth [m]. This is where the "
                    "conductor hangs, already net of insulator string "
                    "length; it is not the crossarm elevation."
                ),
            },
        },
    }


def _tower_type() -> dict[str, Any]:
    return {
        "title": "Tower type",
        "type": "object",
        "additionalProperties": False,
        "required": ["attachments"],
        "description": (
            "Cross-sectional attachment geometry of a representative "
            "structure: one named attachment point per physical conductor, "
            "so bundle sub-conductor positions are separate points. Span "
            "length is route data and lives in the line document's path."
        ),
        "properties": {
            "description": {"type": "string"},
            "source": {"type": "string", "description": _PROVENANCE},
            "attachments": _named_section(
                "#/$defs/attachment",
                "Attachment points by name. Line documents hang conductors "
                "on these names.",
            ),
        },
    }


def _catalog() -> dict[str, Any]:
    return {
        "title": "Catalog",
        "type": "object",
        "additionalProperties": False,
        "minProperties": 1,
        "description": (
            "Named conductor and tower types defined by this document. Same "
            "sections as a catalog file; names of each kind must not "
            "collide with included catalogs."
        ),
        "properties": _catalog_sections(),
    }


def _catalog_document() -> dict[str, Any]:
    return {
        "title": "Catalog document",
        "type": "object",
        "additionalProperties": False,
        "required": ["catalog"],
        "anyOf": [{"required": ["conductors"]}, {"required": ["towers"]}],
        "description": (
            "A standalone catalog file, written in YAML: a header plus the "
            "same sections a line document embeds in its catalog object. "
            "Line documents pull catalog files in through include."
        ),
        "properties": {
            "catalog": {
                "const": "1.0",
                "description": (
                    "Format identifier and version. Present on every "
                    "catalog file."
                ),
            },
            "name": {"type": "string"},
            "description": {"type": "string"},
            "source": {"type": "string", "description": _PROVENANCE},
            **_catalog_sections(),
        },
    }


def _conductor() -> dict[str, Any]:
    return {
        "title": "Conductor",
        "type": "object",
        "additionalProperties": False,
        "required": ["at", "phase", "type"],
        "description": (
            "One physical conductor: an attachment point, a phase label, "
            "and a conductor type."
        ),
        "properties": {
            "at": {
                "type": "string",
                "pattern": NAME,
                "description": (
                    "Name of the tower-type attachment point where this "
                    "conductor hangs. Each attachment point carries at most "
                    "one conductor."
                ),
            },
            "phase": {"$ref": "#/$defs/phase"},
            "circuit": {"$ref": "#/$defs/circuit"},
            "type": {
                "type": "string",
                "pattern": NAME,
                "description": (
                    "Name of a conductor type from this document's catalog "
                    "or an included catalog."
                ),
            },
            "tension": {
                "type": "number",
                "exclusiveMinimum": 0,
                "description": (
                    "Horizontal tension [N]. When present, sag and the "
                    "effective height are computed from the catenary; when "
                    "absent, the attachment height is used directly."
                ),
            },
        },
    }


def _phase() -> dict[str, Any]:
    return {
        "enum": ["a", "b", "c", "n", "g"],
        "description": (
            "a, b, c: phase conductors. n: neutral, a return conductor "
            "carrying load current. g: grounded conductor, such as a shield "
            "or static wire, not intended to carry load current. Conductors "
            "sharing a phase and circuit are at the same potential."
        ),
    }


def _circuit() -> dict[str, Any]:
    return {
        "type": "integer",
        "minimum": 1,
        "default": 1,
        "description": (
            "Circuit index on multi-circuit structures. Meaningful for a, "
            "b and c; readers ignore it for n and g."
        ),
    }


def _path() -> dict[str, Any]:
    return {
        "type": "object",
        "additionalProperties": False,
        "required": ["span"],
        "oneOf": [{"required": ["length"]}, {"required": ["points"]}],
        "description": (
            "Route of the line: the representative support-to-support "
            "span, and exactly one of length or points: an explicit "
            "length, or an ordered GIS centreline whose great-circle "
            "length is scaled by the sagged-to-span ratio."
        ),
        "properties": {
            "span": {
                "type": "number",
                "exclusiveMinimum": 0,
                "description": (
                    "Representative support-to-support span length [m]."
                ),
            },
            "length": {
                "type": "number",
                "exclusiveMinimum": 0,
                "description": "Explicit line length [m].",
            },
            "points": {
                "type": "array",
                "minItems": 2,
                "items": {
                    "type": "object",
                    "additionalProperties": False,
                    "required": ["latitude", "longitude"],
                    "properties": {
                        "latitude": {
                            "type": "number",
                            "minimum": -90,
                            "maximum": 90,
                        },
                        "longitude": {
                            "type": "number",
                            "minimum": -180,
                            "maximum": 180,
                        },
                    },
                },
                "description": (
                    "Ordered GIS centreline [deg]. Intervals are route "
                    "segments, not tower spans."
                ),
            },
        },
    }


def _earth() -> dict[str, Any]:
    return {
        "type": "object",
        "additionalProperties": False,
        "required": ["conductivity", "permittivity"],
        "properties": {
            "conductivity": {
                "type": "number",
                "minimum": 0,
                "description": (
                    "Homogeneous earth conductivity [S/m] for the Carson "
                    "return path."
                ),
            },
            "permittivity": {
                "type": "number",
                "minimum": 8.8541878128e-12,
                "description": (
                    "Earth permittivity [F/m]. At least the vacuum "
                    "permittivity."
                ),
            },
        },
    }


def build_schema() -> dict[str, Any]:
    """Build the complete line schema, ready to serialize."""
    schema: dict[str, Any] = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": FALLBACK_ID,
        "title": "Overhead line",
        "description": (
            "One homogeneous overhead-line section: physical conductors "
            "mapped onto the attachment points of a tower type, the route, "
            "and the earth return. Conductor and tower types come from the "
            "document's own catalog section or from included catalog "
            "files. Field names and units match the EMT parameter model "
            "documentation. All values are SI."
        ),
        "type": "object",
        "additionalProperties": False,
        "required": ["line", "tower", "conductors", "path", "earth"],
        "properties": {
            "line": {
                "const": "1.0",
                "description": (
                    "Format identifier and version. Present on every line "
                    "document."
                ),
            },
            "$schema": {
                "type": "string",
                "format": "uri",
                "description": (
                    "Optional. Points an editor at this schema for "
                    "completion and inline validation."
                ),
            },
            "name": {"type": "string"},
            "description": {"type": "string"},
            "source": {"type": "string", "description": _PROVENANCE},
            "include": {
                "type": "array",
                "minItems": 1,
                "items": {"type": "string"},
                "description": (
                    "Catalog files (YAML), relative to this document. "
                    "Names of each kind must be unique across this "
                    "document's catalog section and every included "
                    "catalog; the loader rejects collisions rather than "
                    "shadowing."
                ),
            },
            "catalog": {"$ref": "#/$defs/catalog"},
            "tower": {
                "type": "string",
                "pattern": NAME,
                "description": (
                    "Name of a tower type from this document's catalog or "
                    "an included catalog. Its attachment points are where "
                    "the conductor entries hang."
                ),
            },
            "conductors": {
                "type": "array",
                "minItems": 1,
                "items": {"$ref": "#/$defs/conductor"},
                "description": (
                    "Every physical conductor on the structure, one entry "
                    "each. Array order defines the row and column order of "
                    "every matrix built from this line. Bundles are not a "
                    "distinct construct: sub-conductors are listed "
                    "individually, each on its own attachment point, and "
                    "are related only by sharing a phase and circuit."
                ),
            },
            "path": {"$ref": "#/$defs/path"},
            "earth": {"$ref": "#/$defs/earth"},
        },
        "$defs": {
            "phase": _phase(),
            "circuit": _circuit(),
            "conductor": _conductor(),
            "conductor_type": _conductor_type(),
            "attachment": _attachment(),
            "tower_type": _tower_type(),
            "catalog": _catalog(),
            "catalog_document": _catalog_document(),
            "path": _path(),
            "earth": _earth(),
        },
    }
    if canonical := os.environ.get("READTHEDOCS_CANONICAL_URL"):
        schema["$id"] = canonical.rstrip("/") + "/line.schema.json"
    schema["$comment"] = (
        "Generated from docs/_ext/gridkit/line_schema.py by the GridKit "
        "Sphinx extension. Edit that module instead of this file."
    )
    return schema


def write_schema(output: Path) -> None:
    output.write_text(
        json.dumps(build_schema(), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
