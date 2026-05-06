# PhasorDynamics Input Format

A PhasorDynamics input file is a JSON object that describes one dynamic phasor
case: metadata, buses, network branches, devices, non-network signals, and
monitor output sinks.

This document defines the portable file-level schema and the canonical component
class strings. Component-specific `params`, `init`, `ports`, and `mon` entries
are defined by the linked model README files in the class catalog.

## General Rules

- Required root arrays must be present. Use `[]` when a category is empty.
- `class` values are case-sensitive and must match a class in the catalog.
- Unknown fields are invalid outside `extension`.
- `extension` objects are implementation-defined. Portable readers may ignore
  them.
- Object identifiers are unique only within their declared category.

## Case Object

The root JSON object has these fields:

| Field | Required | Type | Description |
| --- | --- | --- | --- |
| `header` | Yes | object | Case metadata and system base values. |
| `buses` | Yes | array | Bus objects. |
| `branches` | Yes | array | Network branch objects. |
| `devices` | Yes | array | Non-branch component objects. |
| `signals` | Yes | array | Non-network signal objects. |
| `monitors` | No | array | Monitor sink objects. |
| `extension` | No | object | Implementation-defined root data. |

```json
{
  "header": {
    "format_version": 1,
    "format_revision": 0,
    "case_name": "Example",
    "case_description": "Minimal canonical layout",
    "case_comments": "",
    "freq_base": 60.0,
    "va_base": 100000000.0
  },
  "buses": [],
  "branches": [],
  "devices": [],
  "signals": []
}
```

## Header Object

| Field | Required | Type | Description |
| --- | --- | --- | --- |
| `format_version` | Yes | integer | Major input format version. |
| `format_revision` | Yes | integer | Revision within the major input format version. |
| `case_name` | Yes | string | Case name. |
| `case_date_time` | No | string | ISO 8601 date or datetime. Include a UTC offset when a time is present. |
| `case_description` | Yes | string | Case description. |
| `case_comments` | Yes | string | Additional notes. Use `""` when there are no notes. |
| `freq_base` | Yes | number | System frequency base in Hz. |
| `va_base` | Yes | number | System apparent power base in VA. |
| `extension` | No | object | Implementation-defined header data. |

`format_version` and `format_revision` are non-negative integers. Fractional
versions such as `0.1` are invalid.

## Component Objects

The model topology is split across three component arrays:

- `buses` contains bus objects identified by integer bus number.
- `branches` contains network branch objects identified by string ID.
- `devices` contains non-branch components identified by string ID.

### Bus Object

| Field | Required | Type | Description |
| --- | --- | --- | --- |
| `number` | Yes | integer | Positive bus number, unique within `buses`. |
| `class` | Yes | string | Bus class from the class catalog. |
| `name` | No | string | Label. Names may be empty or non-unique. |
| `params` | Class-defined | object | Bus parameters. |
| `init` | Class-defined | object | Bus initialization values. |
| `mon` | No | array | Ordered monitor variable names. Duplicate entries are invalid. |
| `extension` | No | object | Implementation-defined bus data. |

### Branch and Device Object

| Field | Required | Type | Description |
| --- | --- | --- | --- |
| `id` | Yes | string | Identifier unique within the object's root array. |
| `class` | Yes | string | Component class from the class catalog. |
| `name` | No | string | Label. Names may be empty or non-unique. |
| `params` | Class-defined | object | Component parameters. |
| `init` | Class-defined | object | Component initialization values. |
| `ports` | Class-defined | object | External references as `{ "port_name": target }`. |
| `mon` | No | array | Ordered monitor variable names. Duplicate entries are invalid. |
| `extension` | No | object | Implementation-defined component data. |

Class documentation defines which `params`, `init`, `ports`, and `mon` entries
are required, optional, defaulted, or invalid for each class.

## References and Ports

References are scoped by target category:

| Target category | Identifier |
| --- | --- |
| Bus | Bus `number` integer. |
| Branch | Branch `id` string. |
| Device | Device `id` string. |
| Signal | Signal `signal_id` integer. |

Each port has a class-defined target category. A port value must resolve to an
object in that category. References are not globally unique because the port
contract determines the namespace.

Physical network attachment is represented by ports that target buses, branches,
or devices. Bus voltage and current balance variables are intrinsic model
variables and are not authored as root `signals`.

## Signal Objects

Signals are first-class non-network nets used by signal-targeting ports. Signal
objects do not have `class`, `params`, `init`, `ports`, or `mon` fields.

| Field | Required | Type | Description |
| --- | --- | --- | --- |
| `signal_id` | Yes | integer | Positive signal ID, unique within `signals`. |
| `name` | No | string | Label. |
| `extension` | No | object | Implementation-defined signal data. |

Signals are initialized through connected component models, class-defined
defaults, or explicit source devices.

Unless a component README defines different behavior, signal nets obey these
rules:

- Multiple consumers may share a signal.
- A signal with required consumers must have one producer or a documented
  default.
- Multiple producers for the same signal are invalid unless aggregation is
  class-defined.
- Connected signal ports must agree on dimension.

## Monitor Sinks

Component `mon` entries select class-defined variables for output. The listed
order is the requested output order. Composite monitor variables expand as
defined by the component class.

`monitors` declares output sinks. If `monitors` is omitted, monitor output uses
CSV to standard output with comma delimiter.

| Field | Required | Type | Description |
| --- | --- | --- | --- |
| `file_name` | No | string | Output file path. If omitted, write to standard output. |
| `format` | Yes | string | Output format: `csv`, `json`, or `yaml`. |
| `delim` | No | string | CSV delimiter. Default is `","`. |
| `extension` | No | object | Implementation-defined sink data. |

## Class Catalog

The catalog assigns each component class to the root array where instances must
appear. Classes not listed here are non-canonical unless documented by an
extension.

| Root array | Type | Class | Description | Documentation |
| --- | --- | --- | --- | --- |
| `buses` | Bus | `bus` | Positive-sequence AC phasor bus. | [Bus](Bus/README.md) |
| `buses` | Bus | `infinite_bus` | Positive-sequence AC phasor bus with fixed voltage. | [Bus](Bus/README.md) |
| `branches` | Branch | `Line` | Two-terminal line model. | [Branch](Branch/README.md) |
| `branches` | Branch | `Transformer` | Transformer branch model. Reserved; unimplemented. | [Branch](Branch/README.md) |
| `devices` | Load | `Load` | Static impedance load. | [Load](Load/README.md) |
| `devices` | Load | `LoadZIP` | Static ZIP load. | [LoadZIP](LoadZIP/README.md) |
| `devices` | Machine | `Genrou` | Sixth-order round-rotor synchronous machine. | [GENROU](SynchronousMachine/GENROUwS/README.md) |
| `devices` | Machine | `GenClassical` | Classical synchronous machine. | [GenClassical](SynchronousMachine/GenClassical/README.md) |
| `devices` | Controller | `Tgov1` | TGOV1 turbine governor. | [TGOV1](Governor/Tgov1/README.md) |
| `devices` | Controller | `Ieeet1` | IEEE Type 1 excitation system. | [IEEET1](Exciter/IEEET1/README.md) |
| `devices` | Controller | `SexsPti` | SEXS-PTI simplified excitation system. | [SEXS-PTI](Exciter/SEXS-PTI/README.md) |
| `devices` | Controller | `Ieeest` | IEEEST stabilizer. | [IEEEST](Stabilizer/IEEEST/README.md) |
| `devices` | Event | `BusFault` | Impedance-based fault at a bus. | [BusFault](BusFault/README.md) |

Root `signals` entries declare signal nets and do not use component classes.

## Versioning

Adding a component class normally extends the class catalog and the corresponding
model README without changing the file-level schema. New root fields or
incompatible validation rules require a new `format_revision` or
`format_version`.
