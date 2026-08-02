# Overhead Line Inputs

Complete overhead-line input documents in the `line` 1.0 format defined by
[`line.schema.json`](../../../GridKit/Model/EMT/line.schema.json). Each file
describes one homogeneous line section: conductor data, tower attachment
geometry, route, and earth return. All values are SI. Conductor materials,
tensions, and earth data are illustrative, not normative specifications.

File | Structure
---- | ---------
[`69kv-wood-pole.line.json`](69kv-wood-pole.line.json) | Vertical wood pole with underbuilt neutral, no shield wire
[`138kv-delta.line.json`](138kv-delta.line.json) | Single-circuit delta, one shield wire, GIS route
[`345kv-horizontal.line.json`](345kv-horizontal.line.json) | Horizontal H-frame, twin bundles, tension supplied
[`500kv-double-circuit.line.json`](500kv-double-circuit.line.json) | Double-circuit vertical lattice, twin bundles
[`765kv-horizontal.line.json`](765kv-horizontal.line.json) | Horizontal, four-conductor square bundles

## Tower Cross-Sections

![69 kV wood pole cross-section](figures/69kv-wood-pole.png)

![138 kV delta cross-section](figures/138kv-delta.png)

![345 kV horizontal cross-section](figures/345kv-horizontal.png)

![500 kV double-circuit cross-section](figures/500kv-double-circuit.png)

![765 kV horizontal cross-section](figures/765kv-horizontal.png)

## Rendering

[`render_tower.py`](render_tower.py) regenerates the cross-section figures:

```sh
python render_tower.py *.line.json
```

## Validation

Any JSON Schema 2020-12 validator works, for example:

```sh
check-jsonschema --schemafile ../../../GridKit/Model/EMT/line.schema.json *.line.json
```

The loader additionally checks what the schema cannot express: `tower.x`,
`tower.height`, and optional `tower.tension` have one entry per conductor, in
conductor order.
