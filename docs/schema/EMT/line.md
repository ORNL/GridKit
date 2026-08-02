(schema-emt-line)=
# Overhead Line Input

The `line` 1.0 format describes one homogeneous overhead-line section as four
sections: `conductors` (physical conductor data), `tower` (cross-sectional
attachment geometry), `path` (span and route), and `earth` (Carson return
path). Field names and units match the
[EMT parameter model documentation](../../GridKit/Model/EMT/Parameters/README.md);
all values are SI. Complete instances are in the
[overhead line examples](../../examples/EMT/Lines/README.md).

Use the raw JSON Schema at {{ line_schema_url }}
as a validation endpoint, or
{download}`download a copy <../../../GridKit/Model/EMT/line.schema.json>`.
Editors can enable completion and inline validation for an individual document:

{{ line_schema_example }}

```{jsonschema} ../../../GridKit/Model/EMT/line.schema.json
:lift_description:
:lift_definitions:
:auto_reference:
:auto_target:
```
