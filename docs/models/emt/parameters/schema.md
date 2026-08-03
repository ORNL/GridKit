(models-emt-parameters-schema)=
# JSON Schema Reference

This reference documents the exact JSON Schema accepted by GridKit for
overhead-line input documents. See the
[parameter models](#models-emt-parameters) for the quantities a resolved
document supplies, the [catalog](#catalog) for the shipped conductor and
tower types, and the [overhead line examples](#examples-emt-lines) for
complete inputs.

Use the raw JSON Schema at {{ line_schema_url }}
as a validation endpoint, or
{download}`download a copy <../../../../GridKit/Model/EMT/line.schema.json>`.
Editors can enable validation for an individual document without workspace
settings:

{{ line_schema_example }}

Standalone catalog files validate against the catalog document definition
inside the same schema, `#/$defs/catalog-document`.

:::{div} gk-table

```{jsonschema} ../../../../GridKit/Model/EMT/line.schema.json
:lift_description:
:lift_definitions:
:auto_reference:
:auto_target:
```

:::
