"""The JSON Schema for case files, generated from the model definitions.

The hand-written base holds everything that is not per model. Model classes,
their ports, parameters, and monitorable variables come from the same records
the documentation tables use.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

from .doxygen import JSON_NAMES, Item, Model


def _annotations(item: Item) -> dict[str, str]:
    return {"description": item.description} if item.description else {}


def _closed(properties: dict[str, Any]) -> dict[str, Any]:
    return {
        "type": "object",
        "properties": properties,
        "additionalProperties": False,
    }


def _parameters(model: Model) -> dict[str, Any]:
    return _closed(
        {
            item.name: {"type": ["number", "boolean"], **_annotations(item)}
            for item in model.parameters
        }
    )


def _ports(model: Model) -> dict[str, Any]:
    names = [item.name for item in model.ports]
    if len(names) != len(set(names)):
        raise ValueError(f"{model.name} has duplicate port names")
    return _closed(
        {
            item.name: {
                "$ref": "#/$defs/nonnegative_id",
                **_annotations(item),
            }
            for item in model.ports
        }
    )


def _monitors(model: Model) -> dict[str, Any]:
    schema: dict[str, Any] = {"type": "array", "uniqueItems": True}
    if model.monitors:
        schema["items"] = {"enum": [item.name for item in model.monitors]}
    else:
        schema["maxItems"] = 0
    return schema


def _device(model: Model) -> dict[str, Any]:
    schema = _closed(
        {
            "class": {"const": model.json_name},
            "ports": _ports(model),
            "id": {"type": "string", "minLength": 1},
            "params": _parameters(model),
            "mon": _monitors(model),
            "extension": {"type": "object"},
        }
    )
    schema.update(title=model.name, required=["class", "ports", "id", "params"])
    return schema


def _bus(model: Model) -> dict[str, Any]:
    schema = _closed(
        {
            "number": {"$ref": "#/$defs/nonnegative_id"},
            "class": {"const": model.json_name},
            "name": {"type": "string"},
            "init": {"$ref": "#/$defs/bus_init"},
            "params": _parameters(model),
            "mon": _monitors(model),
            "freq_base": {"type": "number", "exclusiveMinimum": 0},
            "va_base": {"type": "number", "exclusiveMinimum": 0},
            "extension": {"type": "object"},
        }
    )
    schema.update(
        title=model.name,
        required=["number", "class", "name", "params"],
    )
    return schema


def _dispatch(title: str, models: list[Model]) -> dict[str, Any]:
    return {
        "title": title,
        "oneOf": [{"$ref": f"#/$defs/{model.name}"} for model in models],
    }


def build_schema(models: dict[str, Model], base: Path) -> dict[str, Any]:
    buses = [model for model in models.values() if model.kind == "bus"]
    devices = [model for model in models.values() if model.kind == "device"]

    unnamed = {model.name for model in buses} - JSON_NAMES.keys()
    if unnamed:
        raise ValueError(f"missing JSON names for bus models: {sorted(unnamed)}")

    schema = json.loads(base.read_text(encoding="utf-8"))
    if canonical := os.environ.get("READTHEDOCS_CANONICAL_URL"):
        schema["$id"] = canonical.rstrip("/") + "/case.schema.json"
    schema["$comment"] = (
        "Generated from Doxygen XML by the GridKit Sphinx extension. "
        "Edit case.schema.base.json or the C++ model documentation instead."
    )
    definitions = schema["$defs"]
    definitions["bus"] = _dispatch("Bus models", buses)
    definitions["device"] = _dispatch("Device models", devices)
    definitions.update((model.name, _bus(model)) for model in buses)
    definitions.update((model.name, _device(model)) for model in devices)
    return schema


def write_schema(models: dict[str, Model], base: Path, output: Path) -> None:
    output.write_text(
        json.dumps(build_schema(models, base), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
