(cases-phasor-dynamics-schema)=
# JSON Schema Reference

This reference documents the exact JSON Schema accepted by GridKit. See the
[case format](format.md) for conventions and the
[example inputs](format.md#example-file-for-a-2-bus-system) for complete inputs.

Use the raw JSON Schema at {{ schema_url }}
as a validation endpoint, or
{download}`download a copy <../../../GridKit/Model/case.schema.json>`.
Editors can enable validation for an individual case without workspace settings:

{{ schema_example }}

```{jsonschema} ../../../GridKit/Model/case.schema.json
:lift_description:
:lift_definitions:
:auto_reference:
:auto_target:
```
