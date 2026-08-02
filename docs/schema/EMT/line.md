(schema-emt-line)=
# Overhead Line Input

The `line` 1.0 format describes one homogeneous overhead-line section:
`conductors` hangs physical conductors on the named attachment points of a
tower type, `path` gives the span and route, and `earth` gives the Carson
return path. Field names and units match the
[EMT parameter model documentation](../../GridKit/Model/EMT/Parameters/README.md);
all values are SI. Complete instances are in the
[overhead line examples](../../examples/EMT/Lines/README.md).

Conductor data and tower geometry are always referenced by type name.
Types live in catalogs: the document's own `catalog` section for
self-contained files, or shared catalog files listed in `include`, so many
line documents can draw on one library.

## Catalogs

A catalog names reusable types in two sections: `conductors` (dimensions,
material, and weight) and `towers` (the attachment points of a
representative structure, one per physical conductor). Catalog files are
YAML with a `catalog` format header, so libraries stay readable,
commentable, and easy to extend:

```yaml
catalog: "1.0"
name: Common North American overhead types

conductors:
  drake-acsr:
    description: Drake 795 kcmil 26/7 ACSR
    radius: { outer: 0.01407, inner: 0.00514 } # m, aluminum annulus
    conductivity: 3.5e+7                       # S/m
    weight: 16.0                               # N/m

towers:
  345kv-h-frame:
    attachments:
      left-1:      { x: -7.4285, h: 24.0 }
      left-2:      { x: -6.9715, h: 24.0 }
      # ...
      shield-left: { x: -4.3,    h: 29.5 }
```

A line document names its tower type and maps each physical conductor onto
an attachment point:

```json
"tower": "345kv-h-frame",
"conductors": [
  { "at": "left-1",      "phase": "a", "type": "drake-acsr", "tension": 28000.0 },
  { "at": "shield-left", "phase": "g", "type": "ehs-3-8-steel" }
]
```

The same `conductors` and `towers` sections can be embedded in a line
document's `catalog` object, which keeps single-file documents
self-contained. Write YAML exponents with an explicit sign (`3.5e+7`, not
`3.5e7`): YAML 1.1 loaders read unsigned exponents as strings, and
catalog validation rejects them.

Material properties with a vacuum reference are relative: conductor
`permeability` and earth `permittivity` are dimensionless multiples of the
vacuum values, default 1, so nonmagnetic conductors and vacuum-like earth
simply omit them. The loader scales them to absolute SI values during
resolution, and the schema bounds reject absolute values entered by
mistake.

## Generated Schema

The schema is authored in `docs/_ext/gridkit/line_schema.py` and generated
at documentation build time; the repository does not commit the JSON
artifact. Use the raw JSON Schema at {{ line_schema_url }}
as a validation endpoint. Editors can enable completion and inline
validation for an individual document:

{{ line_schema_example }}

Standalone catalog files validate against the catalog document definition
inside the same schema, `#/$defs/catalog-document`.

## Semantic Validation

The loader enforces what JSON Schema cannot express:

- Every `include` path resolves, relative to the document, to a valid
  YAML catalog file.
- Names of each kind are unique across the document's `catalog` section
  and all included catalogs; collisions are rejected rather than shadowed.
- The `tower` reference and every conductor `type` reference name a
  defined type.
- Every conductor `at` names an attachment point of the tower type, and
  no attachment point carries more than one conductor.
- Every conductor type's inner radius is below its outer radius.
- When `tension` is supplied, the conductor weight and `path.span` give a
  positive minimum conductor height.

Resolution replaces every reference with its catalog data and produces one
flat record per physical conductor — attachment coordinates, phase,
circuit, tension, and the conductor type data with relative material
properties scaled to absolute SI, in `conductors` order. These records are
the quantities the
[parameter models](../../GridKit/Model/EMT/Parameters/README.md) consume;
type names do not survive resolution.

## Schema Reference

```{jsonschema} ../../../GridKit/Model/EMT/line.schema.json
:lift_description:
:lift_definitions:
:auto_reference:
:auto_target:
```
