# sphinx-latkit

A pip package, published from the latkit repo, that lets any Sphinx site embed
interactive `<latkit-network>` diagrams. The wheel ships the compiled
`embed.bundle.js` and its assets built from the same commit, so one pin in
`docs/requirements.txt` pins the Python mechanics and the exact JavaScript
together. latkit is a general network tool; this package knows nothing about
any domain — domain projects subclass one directive class.

## Design

- **One custom doctree node carries each diagram's payload** from read to
  translation. Sphinx's own pipeline handles caching, incremental rebuilds,
  parallelism, and purging — no environment state, no `env-purge-doc` /
  `env-merge-info` / `build-finished` handlers.
- **Sidecars are written at HTML translation time**, content-addressed
  (`_latkit/<sha1>.json`), straight into the output directory. Writes are
  idempotent, so parallel writers and re-runs are safe by construction, and a
  sidecar exists exactly when a page referencing it was built.
- **The `src` URI is computed by the builder that is actually writing the
  page** (`relative_uri` against `get_target_uri`), so `html`, `dirhtml`, and
  `epub` are all correct from one code path, and cached doctrees stay
  builder-agnostic.
- **The extension API is stock Sphinx:** subclass `NetworkDirective`, override
  `compile`, register with `app.add_directive`. No bespoke registration
  function, no dynamic class creation.

## Layout (latkit repo)

```
sphinx/
  pyproject.toml
  README.md
  sphinx_latkit/
    __init__.py
    static/latkit/host.css        # bundle + assets force-included at wheel build
  tests/
    test_build.py
```

## `sphinx_latkit/__init__.py`

```python
"""Sphinx integration for @latkit/embed: directives that emit <latkit-network>.

The wheel ships the element bundle and its assets. Each diagram travels through
the doctree as a `network` node and is written at translation time as a
content-addressed sidecar under _latkit/, fetched lazily by the element.
"""

import hashlib
import html
import json
from importlib.metadata import version
from importlib.resources import files
from pathlib import Path

from docutils import nodes
from docutils.parsers.rst import directives
from sphinx.util.docutils import SphinxDirective
from sphinx.util.osutil import relative_uri

#: Live view attributes of <latkit-network>; mirrors VIEW_ATTRIBUTES in @latkit/embed.
VIEW_ATTRIBUTES = ("colormap", "projection", "vertex-color", "vertex-height",
                   "vertex-size", "edge-color", "edge-dash")

FALLBACK = "<p>Interactive diagram requires a WebGPU browser.</p>"


class network(nodes.General, nodes.Element):
    """Carries one diagram's serialized payload and options until translation."""


class NetworkDirective(SphinxDirective):
    """Embed a network-data JSON file. Subclasses override `compile` for domain formats."""

    required_arguments = 1
    final_argument_whitespace = True
    option_spec = {"height": directives.unchanged,
                   **dict.fromkeys(VIEW_ATTRIBUTES, directives.unchanged)}

    def compile(self, source):
        return source

    def run(self):
        _, path = self.env.relfn2path(self.arguments[0])
        self.env.note_dependency(path)
        try:
            data = self.compile(json.loads(Path(path).read_text("utf-8")))
        except Exception as error:
            raise self.error(f"{self.name}: {self.arguments[0]}: {error}")

        node = network(payload=json.dumps(data, separators=(",", ":")),
                       options=dict(self.options))
        self.set_source_info(node)
        return [node]


def visit_network_html(translator, node):
    builder = translator.builder
    payload = node["payload"]
    digest = hashlib.sha1(payload.encode()).hexdigest()[:16]

    sidecar = Path(builder.outdir, "_latkit", f"{digest}.json")
    if not sidecar.exists():
        sidecar.parent.mkdir(parents=True, exist_ok=True)
        sidecar.write_text(payload, "utf-8")

    src = relative_uri(builder.get_target_uri(builder.current_docname),
                       f"_latkit/{digest}.json")
    options = node["options"]
    attributes = "".join(f' {name}="{html.escape(options[name])}"'
                         for name in VIEW_ATTRIBUTES if name in options)
    height = html.escape(options.get("height", "480px"))
    translator.body.append(f'<latkit-network src="{src}"{attributes} '
                           f'style="height:{height}">{FALLBACK}</latkit-network>')
    raise nodes.SkipNode


def skip_network(translator, node):
    raise nodes.SkipNode


def add_static_path(app):
    app.config.html_static_path.append(str(files(__name__) / "static"))


def setup(app):
    app.add_node(network, html=(visit_network_html, None),
                 **dict.fromkeys(("latex", "text", "man", "texinfo"),
                                 (skip_network, None)))
    app.add_directive("latkit-network", NetworkDirective)
    app.add_js_file("latkit/embed.bundle.js", type="module")
    app.add_css_file("latkit/host.css")
    app.connect("builder-inited", add_static_path)
    return {"version": version("sphinx-latkit"),
            "parallel_read_safe": True,
            "parallel_write_safe": True}
```

## `sphinx_latkit/static/latkit/host.css`

```css
/* Host defaults for <latkit-network>. display:block applies before element
   upgrade and without JavaScript, so height and fallback render correctly. */
latkit-network {
    display: block;
    margin: 1em 0 1.5em;
    border: 1px solid rgba(0, 0, 0, 0.12);
    border-radius: 6px;
    background: #fdfdfd;
    overflow: hidden;
}
latkit-network > p {
    display: flex;
    align-items: center;
    justify-content: center;
    height: 100%;
    margin: 0;
    color: #555;
}
```

## `sphinx/pyproject.toml`

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "sphinx-latkit"
version = "0.0.0"  # stamped from packages/embed/package.json by the release job
description = "Sphinx directives for interactive @latkit/embed network diagrams."
readme = "README.md"
license = "MIT"
requires-python = ">=3.10"
dependencies = ["sphinx>=7"]

[tool.hatch.build.targets.wheel]
packages = ["sphinx_latkit"]

# The element bundle and assets come from the pnpm build; wheel builds fail
# fast when packages/embed/dist is missing.
[tool.hatch.build.targets.wheel.force-include]
"../packages/embed/dist/embed.bundle.js" = "sphinx_latkit/static/latkit/embed.bundle.js"
"../packages/embed/dist/assets" = "sphinx_latkit/static/latkit/assets"
```

## Prerequisite in `@latkit/embed`

The self-contained bundle entry (`dist/register.js` has bare specifiers a
browser cannot resolve):

```ts
// tsup.config.ts — second config alongside the existing one
{
  entry: { 'embed.bundle': 'src/register.ts' },
  format: ['esm'],
  noExternal: [/^@latkit\//],
  minify: true,
  sourcemap: true,
  target: 'es2022',
  platform: 'browser',
}
```

## Release (appended to the existing latkit release job)

PyPI trusted publishing supports pending publishers, so the OIDC connection is
configured before the first release; no token exists at any point.

```yaml
      - name: Build sphinx-latkit
        if: steps.changesets.outputs.published == 'true' && contains(steps.changesets.outputs.publishedPackages, '@latkit/embed')
        run: |
          VERSION=$(node -p "require('./packages/embed/package.json').version")
          sed -i "s/^version = .*/version = \"$VERSION\"/" sphinx/pyproject.toml
          pipx run build --wheel sphinx/

      - name: Publish sphinx-latkit to PyPI
        if: steps.changesets.outputs.published == 'true' && contains(steps.changesets.outputs.publishedPackages, '@latkit/embed')
        uses: pypa/gh-action-pypi-publish@release/v1
        with:
          packages-dir: sphinx/dist
```

The wheel version equals the `@latkit/embed` version by construction.
latkit's own docs install `./sphinx` as a path dependency and embed a demo via
the plain `latkit-network` directive, so every docs build exercises the wheel.

## Tests (`sphinx/tests/test_build.py`, pytest + `sphinx.testing`)

- Emitted element markup: `src` path, passthrough attributes, escaping, height.
- Sidecar `_latkit/<digest>.json` written and byte-identical to the payload.
- `dirhtml` resolves `src` one level deeper than `html` from the same sources.
- Incremental rebuild after deleting the output directory restores sidecars.
- A failing `compile` surfaces as a directive error naming the source file.
- `latex` build succeeds and emits nothing for the directive.

## Consumer: GridKit

`docs/requirements.txt` gains `sphinx-latkit==0.1.0`. Deleted: `docs/js/`,
`docs/_static/js/oneline.js`, `docs/_static/css/oneline.css`. Reverted:
`.readthedocs.yaml` (no Node, no npm) and `.gitignore`. Case pages and
`conf.py` are unchanged. The feature is one file:

```python
# docs/_ext/oneline.py — the entire GridKit-side feature
"""`oneline` directive: embed a GridKit case as an interactive network diagram."""

import math

from sphinx_latkit import NetworkDirective


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def param(element, key):
    return (element.get("params") or {}).get(key)


def voltage(bus):
    init = bus.get("init") or {}
    vr, vi = init.get("Vr"), init.get("Vi")
    return round(math.hypot(vr, vi), 5) if finite(vr) and finite(vi) else None


def bus_label(bus):
    name = f"Bus {bus.get('name', bus['number'])}"
    kv = param(bus, "kv")
    return f"{name} — {kv:g} kV" if finite(kv) else name


def field(fid, label, unit, scope, values):
    # @latkit/embed rejects non-finite values, so columns are all-or-nothing.
    if values and all(finite(v) for v in values):
        return {"id": fid, "label": label, "unit": unit, "scope": scope, "values": values}


def compile_case(case):
    """GridKit case -> latkit network data (topology, labels, vm/kv/r/x fields)."""
    buses = case["buses"]
    rows = {bus["number"]: row for row, bus in enumerate(buses)}
    if len(rows) != len(buses):
        raise ValueError("duplicate bus numbers")
    branches = [dev for dev in case.get("devices", ()) if dev.get("class") == "Branch"]

    def endpoint(branch, port):
        number = (branch.get("ports") or {}).get(port)
        if number not in rows:
            raise ValueError(f"branch {branch.get('id')!r}: port {port!r} references unknown bus {number}")
        return rows[number]

    topology = {"vertexCount": len(buses),
                "edges": [endpoint(branch, port) for branch in branches for port in ("bus1", "bus2")]}
    coords = [(bus.get("extension") or {}).get(axis)
              for bus in buses for axis in ("longitude", "latitude")]
    if coords and all(finite(coord) for coord in coords):
        topology["vertexCoords"] = [round(coord, 6) for coord in coords]

    fields = [column for column in (
        field("vm", "Voltage", "pu", "vertex", [voltage(bus) for bus in buses]),  # first => default vertex-color
        field("kv", "Base voltage", "kV", "vertex", [param(bus, "kv") for bus in buses]),
        field("r", "Resistance", "pu", "edge", [param(branch, "R") for branch in branches]),
        field("x", "Reactance", "pu", "edge", [param(branch, "X") for branch in branches]),
    ) if column]

    data = {"topology": topology, "labels": {"vertex": [bus_label(bus) for bus in buses]}}
    return data | ({"fields": fields} if fields else {})


class OnelineDirective(NetworkDirective):
    def compile(self, case):
        return compile_case(case)


def setup(app):
    app.setup_extension("sphinx_latkit")
    app.add_directive("oneline", OnelineDirective)
    return {"version": "2.0", "parallel_read_safe": True, "parallel_write_safe": True}
```

Authoring is unchanged, with every element attribute available as an option:

````markdown
```{oneline} ../../../../../examples/PhasorDynamics/Large/WECC/wecc.json
:height: 480px
```
````

emits:

```html
<latkit-network src="../../_latkit/3fbb90d2c1a44e19.json" style="height:480px">
  <p>Interactive diagram requires a WebGPU browser.</p>
</latkit-network>
```

`vm` is the first vertex field, so voltage coloring applies with zero
configuration; upgrades reach GridKit as a one-line pip version bump.
