from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

from .inventory import Item, Model, read_models


BUS_CLASS_NAMES = {"Bus": "bus", "BusInfinite": "infinite_bus"}


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
    items = (*model.buses, *model.signal_inputs, *model.signal_outputs)
    names = [item.name for item in items]
    if len(names) != len(set(names)):
        raise ValueError(f"{model.name} has duplicate port names")
    return _closed(
        {
            item.name: {
                "$ref": "#/$defs/nonnegative_id",
                **_annotations(item),
            }
            for item in items
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
            "class": {"const": model.name},
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
            "class": {"const": BUS_CLASS_NAMES[model.name]},
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


def build_schema(xml: Path, base: Path) -> dict[str, Any]:
    models = list(read_models(xml).values())
    buses = [model for model in models if model.kind == "bus"]
    devices = [model for model in models if model.kind == "device"]

    unknown_buses = {model.name for model in buses} - BUS_CLASS_NAMES.keys()
    if unknown_buses:
        raise ValueError(f"missing JSON names for bus models: {sorted(unknown_buses)}")

    schema = json.loads(base.read_text(encoding="utf-8"))
    if canonical := os.environ.get("READTHEDOCS_CANONICAL_URL"):
        schema["$id"] = canonical.rstrip("/") + "/case.schema.json"
    schema["$comment"] = (
        "Generated from Doxygen XML by the GridKit Sphinx schema extension. "
        "Edit case.schema.base.json or the C++ model documentation instead."
    )
    definitions = schema["$defs"]
    definitions["bus"] = _dispatch("Bus models", buses)
    definitions["device"] = _dispatch("Device models", devices)
    definitions.update((model.name, _bus(model)) for model in buses)
    definitions.update((model.name, _device(model)) for model in devices)
    return schema


def write_schema(xml: Path, base: Path, output: Path) -> None:
    schema = build_schema(xml, base)
    output.write_text(
        json.dumps(schema, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
