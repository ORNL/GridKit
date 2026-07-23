"""Sphinx directive that embeds an interactive oneline diagram.

Usage in a MyST page::

    ```{oneline} ../../../../../examples/PhasorDynamics/Large/WECC/wecc.json
    :height: 480px
    :colormap: viridis
    ```

The directive reads the GridKit case file at build time, converts buses and
branches into the latkit topology arrays, and inlines the result into the
page as JSON. The client-side loader (``_static/js/oneline.js``) renders it
with WebGPU; when WebGPU is unavailable the container hides itself and the
page's static figures stand alone.
"""

import json
import math
import re
from pathlib import Path

from docutils import nodes
from docutils.parsers.rst import directives
from sphinx.util.docutils import SphinxDirective

_HEIGHT_RE = re.compile(r"^\d+(\.\d+)?(px|em|rem|vh|%)$")
_COLORMAP_RE = re.compile(r"^[a-z0-9_-]+$")


def _height(argument):
    value = argument.strip()
    if not _HEIGHT_RE.match(value):
        raise ValueError(f"invalid height {value!r} (use e.g. 480px, 32em, 60vh)")
    return value


def _colormap(argument):
    value = argument.strip()
    if not _COLORMAP_RE.match(value):
        raise ValueError(f"invalid colormap name {value!r}")
    return value


def _number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _topology_payload(case):
    """Convert a GridKit case dict into the payload the loader consumes.

    Mirrors the bus/branch reading used by the lattice tooling: buses become
    vertices in declaration order, ``Branch`` devices become edges, and
    coordinates are taken from ``extension.longitude``/``latitude`` only when
    every bus has them.
    """
    header = case.get("header")
    if not isinstance(header, dict):
        header = {}
    buses = case["buses"]
    index = {}
    names = []
    kv = []
    vm = []
    coords = []
    for i, bus in enumerate(buses):
        index[bus["number"]] = i
        names.append(str(bus.get("name", bus["number"])))
        params = bus.get("params") or {}
        kv.append(params.get("kv") if _number(params.get("kv")) else None)
        init = bus.get("init") or {}
        vr, vi = init.get("Vr"), init.get("Vi")
        vm.append(round(math.hypot(vr, vi), 5) if _number(vr) and _number(vi) else None)
        ext = bus.get("extension") or {}
        lon, lat = ext.get("longitude"), ext.get("latitude")
        if _number(lon) and _number(lat):
            coords.extend((lon, lat))

    edges = []
    for device in case.get("devices", ()):
        if device.get("class") != "Branch":
            continue
        ports = device.get("ports") or {}
        try:
            edges.extend((index[ports["bus1"]], index[ports["bus2"]]))
        except KeyError as error:
            raise ValueError(
                f"branch {device.get('id')!r} references unknown bus {error}"
            ) from None

    return {
        "name": header.get("case_name") or "",
        "vertexCount": len(buses),
        # Omitted when any bus lacks extension.{longitude,latitude}; the
        # renderer then falls back to its generated layout.
        "coords": coords if buses and len(coords) == 2 * len(buses) else None,
        "edges": edges,
        "names": names,
        "kv": kv,
        "vm": vm if vm and all(v is not None for v in vm) else None,
    }


class OnelineDirective(SphinxDirective):

    required_arguments = 1
    final_argument_whitespace = True
    has_content = False
    option_spec = {
        "height": _height,
        "colormap": _colormap,
    }

    def run(self):
        _, case_path = self.env.relfn2path(self.arguments[0])
        self.env.note_dependency(case_path)
        try:
            case = json.loads(Path(case_path).read_text(encoding="utf-8"))
            payload = _topology_payload(case)
        except (OSError, ValueError, KeyError, TypeError) as error:
            raise self.error(
                f"oneline: cannot build topology from {self.arguments[0]}: {error}"
            )

        height = self.options.get("height", "480px")
        colormap = self.options.get("colormap", "viridis")
        # "</" would terminate the inline script element early.
        data = json.dumps(payload, separators=(",", ":")).replace("</", "<\\/")
        html = (
            f'<div class="oneline" data-colormap="{colormap}" style="height:{height}">'
            f'<script type="application/json">{data}</script>'
            f"</div>"
        )
        return [nodes.raw("", html, format="html")]


def setup(app):
    app.add_directive("oneline", OnelineDirective)
    app.add_js_file("js/oneline.js", type="module")
    app.add_css_file("css/oneline.css")
    return {
        "version": "0.1.0",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
