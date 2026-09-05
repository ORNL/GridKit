# EMT case format specification

## Overview

This document describes the JSON data format for EMT cases. The format
follows the grid dynamics case format described in
[PhasorDynamics INPUT_FORMAT](../PhasorDynamics/INPUT_FORMAT.md), specialized
to instantaneous phase-coordinate models: every quantity is in SI units, and
there is no system-wide power or frequency base.

## Format

The root object contains `header` and `devices`, and may also contain
`signals` and `monitors`. `header` contains information about the case and
`devices` is an array of system components, including buses and nested
containers. A Container recursively has its own `devices` and `signals`, so
the same model-body format is used at every level.

### Header

Contained in the `header` key is an object with the following items:

   Name              | Value
 --------------------|-------------------------------------------------------
  `format_version`   | Optional number identifying the case format version
  `format_revision`  | Optional integer identifying the case format revision
  `case_name`        | A string containing the name of the case
  `case_date_time`   | Optional string in the ISO 8601 format indicating a datetime associated with the case
  `case_description` | A string with more specific description of what is modeled in the case
  `case_comments`    | A string with additional notes as needed

### Monitors

Contained in the `monitors` key is an array of objects, each of which
describes an output for monitored variables (those listed in the `mon` field
of a [device](#devices)). The following fields are supported:

  Name               | Description
  -------------------|------------------------------------------------------
  `file_name`        | Optional string indicating output file name. If omitted, `stdout` is used.
  `format`           | One of { "CSV", "JSON", "YAML" } (case-insensitive)
  `delim`            | Optional string specifying delimiter to use for CSV output (default is `","`).

### Signals

Contained in the optional `signals` key is an array of objects, each of which
represents a scalar signal connecting device signal ports:

  Name               | Description
  -------------------|------------------------------------------------------
  `id`               | Nonempty string uniquely identifying the signal in its Container

### Devices

Contained in the `devices` key is an array of objects, each of which
represents a device and has the following fields:

  Name               | Description
  -------------------|------------------------------------------------------
  `class`            | A string indicating the device class. See the table below
  `id`               | Nonempty string uniquely identifying the component in its Container
  `params`           | Optional object mapping parameter names to values, using the JSON names from the device model's Model Parameters table
  `inputs`           | Optional object mapping supported input keys to component or signal IDs
  `outputs`          | Optional object mapping supported output keys to signal IDs
  `mon`              | Optional array of variables to monitor, from the device model's Monitors table

Electrical inputs such as `bus`, `bus1`, and `bus2` refer to Bus component
IDs. Scalar control inputs refer to signal IDs, and scalar outputs assign the
component output to a signal ID. Device and signal IDs share one local
namespace and must be unique within their Container. The same local ID may be
reused in a different Container. Local IDs must not contain `.`, which is
reserved for a child boundary reference such as `plant.speed`.

Electrical current injection and bus-voltage sharing are established together
when an electrical input is connected. They do not require a second case-file
output entry. Consequently, a Bus has no `inputs` or `outputs` entries in case
JSON even though its model equations exchange voltage and current with the
connected devices.

A Bus is an ordinary device. It may additionally contain an optional `init`
object with the initial instantaneous phase voltages `va`, `vb`, and `vc` in
volts; missing entries default to zero.

Legacy top-level `buses` arrays and per-device `ports` objects are not
supported. Move each Bus into `devices`, then split every old `ports` entry
between the `inputs` and `outputs` objects according to the table below.

#### Containers

A Container is a Component whose state is the concatenation of its child
Components. It owns no physical variables or equations itself. Its optional
`inputs` and `outputs` objects define the names visible to its parent:

- `inputs` maps a public name to its source endpoint in the parent scope. The
  public name becomes an endpoint inside the Container and may be used by its
  devices as a scalar Signal or three-phase electrical port.
- `outputs` maps a public name to a Signal or Bus inside the Container. A
  component in the parent may consume it with an input reference such as
  `child.speed` or `child.terminal`.

An input source may be a local Signal or Bus in the parent, an input imported
by the parent, or an explicitly exposed child output. An output may expose a
local Signal or Bus or forward an explicitly exposed child output. References
do not search ancestor scopes and cannot access a child's private devices or
signals. The root SystemModel cannot have inputs because it has no parent
scope.

For example, these two independently scoped systems each own a Bus called
`bus`, expose those buses as terminals, and are wired by a line in the parent:

```json
{
  "devices": [
    {
      "class": "Container",
      "id": "left",
      "outputs": { "terminal": "bus" },
      "devices": [
        { "class": "Bus", "id": "bus" },
        {
          "class": "VoltageSource",
          "id": "source",
          "inputs": { "bus": "bus" }
        }
      ]
    },
    {
      "class": "Container",
      "id": "right",
      "outputs": { "terminal": "bus" },
      "devices": [
        { "class": "Bus", "id": "bus" },
        {
          "class": "LoadZ",
          "id": "load",
          "inputs": { "bus": "bus" }
        }
      ]
    },
    {
      "class": "LineLumped",
      "id": "tie",
      "inputs": {
        "bus1": "left.terminal",
        "bus2": "right.terminal"
      }
    }
  ]
}
```

`left.terminal` is the exact three-phase port owned by `left`'s Bus. A
boundary adds no state or residual equation. Two independently owned Bus
variables therefore cannot be identified by an alias; connect them with an
explicit Line, Switch, transformer, or constraint Component.

Inputs likewise preserve endpoint identity. This controller reads the exact
Signal exported by `plant`; no relay Signal or copy equation is introduced:

```json
{
  "class": "Container",
  "id": "controller",
  "inputs": { "speed": "plant.speed" },
  "outputs": { "pmech": "pmech" },
  "signals": [{ "id": "pmech" }],
  "devices": [
    {
      "class": "Tgov1",
      "id": "governor",
      "inputs": { "speed": "speed" },
      "outputs": { "pmech": "pmech" }
    }
  ]
}
```

The same rule applies to electrical inputs. A child may bind
`"inputs": {"terminal": "parent_bus"}` and connect an internal LoadZ to
`"terminal"`; both refer to the exact Port3 owned by `parent_bus`.

An inline Container may contain only `id`, `class`, `inputs`, `outputs`,
`signals`, and `devices`. It does not have its own header, monitor sinks,
parameters, or monitored variables. The recursive body is also the basis for
future file-backed Containers; file inclusion is not part of this revision.

#### Case connections

  Class                    | Key     | Direction | Target kind   | Required
  -------------------------|---------|-----------|---------------|---------
  `Bus`                    | —       | —         | —             | —
  `VoltageSource`          | `bus`   | Input     | Bus component | Yes
  `DependentVoltageSource` | `bus`   | Input     | Bus component | Yes
  `DependentVoltageSource` | `ea`    | Input     | Signal        | Yes
  `DependentVoltageSource` | `eb`    | Input     | Signal        | Yes
  `DependentVoltageSource` | `ec`    | Input     | Signal        | Yes
  `Machine`                | `bus`   | Input     | Bus component | Yes
  `Machine`                | `pm`    | Input     | Signal        | No
  `Machine`                | `efd`   | Input     | Signal        | No
  `Machine`                | `speed` | Output    | Signal        | No
  `LineLumped`             | `bus1`  | Input     | Bus component | Yes
  `LineLumped`             | `bus2`  | Input     | Bus component | Yes
  `LoadZ`                  | `bus`   | Input     | Bus component | Yes
  `Tgov1`                  | `speed` | Input     | Signal        | No
  `Tgov1`                  | `pref`  | Input     | Signal        | No
  `Tgov1`                  | `pmech` | Output    | Signal        | Yes
  `Switch`                 | `bus1`  | Input     | Bus component | Yes
  `Switch`                 | `bus2`  | Input     | Bus component | Yes
  `Container`              | user-defined | Input/output | Public boundary | No

Declaring a signal creates a named connection. Its value is supplied by a
component output or by the embedding program. Omit an optional input to use
the model's internal default or latched value.

#### Device classes

  Class                 | Model
  ----------------------|------------------------------------------------------
  `Bus`                 | [Bus](Component/Bus/README.md)
  `DependentVoltageSource` | [DependentVoltageSource](Component/Source/DependentVoltageSource/README.md)
  `VoltageSource`       | [VoltageSource](Component/Source/VoltageSource/README.md)
  `LineLumped`          | [LineLumped](Component/Line/LineLumped/README.md)
  `LoadZ`               | [LoadZ](Component/Load/LoadZ/README.md)
  `Switch`              | [Switch](Component/Switch/README.md)
  `Machine`              | [Machine](Component/Source/Machine/README.md)
  `Tgov1`                | [TGOV1](Component/Controller/TGOV1/README.md)
  `Container`            | Recursive collection of devices and signals

#### Parameter values

Parameter values are typed by their JSON representation:

- A boolean, a number with a decimal point, or an integer maps to a Boolean,
  real, or integer parameter respectively; real-valued scalar parameters must
  be written with a decimal point.
- A length-3 array maps to a three-phase vector parameter, integer-valued
  when every element is an integer.
- A 3x3 nested array maps to a real three-phase matrix parameter.
