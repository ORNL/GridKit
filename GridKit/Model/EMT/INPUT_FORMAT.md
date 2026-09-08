# EMT case format specification

## Overview

This document describes the JSON data format for EMT cases. The format
follows the grid dynamics case format described in
[PhasorDynamics INPUT_FORMAT](../PhasorDynamics/INPUT_FORMAT.md), specialized
to instantaneous phase-coordinate models. Quantities use each model's documented
units; there is no system-wide power or frequency base.

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
  `value`            | Optional finite number supplying a constant signal; must not also have a component producer

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

All model inputs are scalar signals. Electrical models read `va`, `vb`, and
`vc`; two-terminal models read `v1a`, `v1b`, `v1c`, `v2a`, `v2b`, and `v2c`.
Voltage signals expose the bus values and derivatives. The bus registers
connected-device current signals and owns their KCL contributions; electrical
connections do not require separate current-output wiring.

The case parser accepts `"bus": "b"` as shorthand for
`"va": "b.va", "vb": "b.vb", "vc": "b.vc"`. Likewise `bus1` and `bus2`
expand to the corresponding terminal phases. The shortcut and an explicit
phase for that terminal cannot appear together. Only scalar phase mappings
remain in model data. A Container can export each phase, for example
`plant.va`, as an ordinary scalar output.

Bus `outputs` can publish `va`, `vb`, and `vc` on named signals. Its optional
`ia`, `ib`, and `ic` inputs add externally supplied injections to the same KCL
rows. Physical components already inject through their voltage connections;
do not add those contributions a second time with explicit Bus inputs.

Bus `shunts` contains named VectorFit admittances with zero incident current:

```json
{
  "class": "Bus", "id": "pcc",
  "shunts": {
    "filter": {"E": [[2e-5,0,0],[0,2e-5,0],[0,0,2e-5]]}
  }
}
```

Assembly splits line `Yp` (or `Gp`/`Cp`) between two bus-owned Norton sources
at scale `dx/2`. Each terminal preserves its scalar aliases and phase order.
The bus exposes `<name>_Ish_a`, `b`, `c` outputs and `<name>_inc_a`, `b`, `c`
inputs. `LineLumped` supplies `i21 = -i12` to terminal 1 and `i12` to terminal 2.

`LineDistributed` similarly creates one bus-owned characteristic admittance
per terminal from `submodels.Yc`, at unit scale. Its `submodels.H` contains
`K` and a nonempty `modes` array, each entry holding a positive `tau` and a
proper square VectorFit `H`. The line supplies each bus's incident current
from the opposite terminal's reflected current through an independent
propagation instance. See [Propagation](Operators/Shift/Propagation/README.md).

Initial values belong exclusively in the [state file](STATE.md), keyed by
component path and existing output names.
Bus voltages default to zero. Model initialization reconstructs
internal variables from known outputs and attached inputs.

Device and signal IDs share one local namespace and must be unique within
their Container. Local IDs cannot contain `.`, which separates component,
boundary, and voltage-phase references.

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
      "outputs": { "va": "bus.va", "vb": "bus.vb", "vc": "bus.vc" },
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
      "outputs": { "va": "bus.va", "vb": "bus.vb", "vc": "bus.vc" },
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
        "bus1": "left",
        "bus2": "right"
      }
    }
  ]
}
```

`left.va`, `left.vb`, and `left.vc` reference the phase signals owned by `left`'s Bus. A
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

Electrical boundaries use scalar signals too. A child binds
`"inputs": {"va": "parent_bus.va", "vb": "parent_bus.vb", "vc": "parent_bus.vc"}`.
An internal LoadZ connects with `"inputs": {"va": "va", "vb": "vb", "vc": "vc"}`.

An inline Container may contain only `id`, `class`, `inputs`, `outputs`,
`signals`, and `devices`. It does not have its own header, monitor sinks,
parameters, or monitored variables. The recursive body is also the basis for
future file-backed Containers; file inclusion is not part of this revision.

#### Case connections

  Class                    | Key     | Direction | Target kind   | Required
  -------------------------|---------|-----------|---------------|---------
  `Bus`                    | `va`, `vb`, `vc` | Output | Signal | No
  `Bus`                    | `ia`, `ib`, `ic` | Input | Signal | No
  `VoltageSource`          | `va`, `vb`, `vc` | Input | Voltage signal | Yes
  `DependentVoltageSource` | `va`, `vb`, `vc` | Input | Voltage signal | Yes
  `DependentVoltageSource` | `ea`    | Input     | Signal        | Yes
  `DependentVoltageSource` | `eb`    | Input     | Signal        | Yes
  `DependentVoltageSource` | `ec`    | Input     | Signal        | Yes
  `PWM`                    | `m`     | Input     | Three Signal IDs | No
  `PWM`                    | `s`     | Output    | Three Signal IDs | No
  `OuterPowerControl`      | `i`, `ilim` | Input | Two Signal IDs | Yes
  `OuterPowerControl`      | `icmd` | Output | Two Signal IDs | No
  `PLL`                    | `va`, `vb`, `vc` | Input | Voltage signal | Yes
  `PLL`                    | `theta`, `omega` | Output | Signal | No
  `Park`                   | `input` | Input     | Three Signal IDs | Yes
  `Park`                   | `theta` | Input     | Signal | Yes
  `Park`                   | `out` | Output | Three Signal IDs | No
  `Angle`                  | `omega` | Input | Signal | Yes
  `Angle`                  | `theta` | Output | Signal | No
  `Modulation`             | `u` | Input | Three Signal IDs | Yes
  `Modulation`             | `vdc` | Input | Signal | Yes
  `Modulation`             | `m` | Output | Three Signal IDs | No
  `InnerCurrentControl`    | `v`, `i`, `icmd` | Input | Two Signal IDs | Yes
  `InnerCurrentControl`    | `omega`, `vdc` | Input | Signal | Yes
  `InnerCurrentControl`    | `ilim`, `u` | Output | Two Signal IDs | No
  `OuterVoltageControl`    | `vref`, `v`, `ig`, `ilim` | Input | Two Signal IDs | Yes
  `OuterVoltageControl`    | `omega` | Input | Signal | Yes
  `OuterVoltageControl`    | `icmd` | Output | Two Signal IDs | No
  `DCLink`                 | `isrc`, `idc` | Input | Signal | Yes
  `DCLink`                 | `vdc` | Output | Signal | No
  `Converter`              | `s`     | Input     | Three Signal IDs | Yes
  `Converter`              | `vdc`   | Input     | Signal        | Yes
  `Converter`              | `i`     | Input     | Three Signal IDs | Yes
  `Converter`              | `idc`   | Output    | Signal        | No
  `Converter`              | `vo`    | Output    | Three Signal IDs | No
  `Machine`                | `va`, `vb`, `vc` | Input | Voltage signal | Yes
  `Machine`                | `pm`    | Input     | Signal        | No
  `Machine`                | `efd`   | Input     | Signal        | No
  `Machine`                | `speed`, `ia`, `ib`, `ic` | Output | Signal | No
  `REGFMA`                 | `v` | Input | Three voltage signal IDs | Yes
  `REGFMA`                 | `pref`, `qref`, `vref` | Input | Signal | No
  `REGFMA`                 | `i` | Output | Three Signal IDs | No
  `VoltageSource`, `DependentVoltageSource`, `LoadZ` | `ia`, `ib`, `ic` | Output | Signal | No
  `LineLumped`, `Switch` | `i12a`, `i12b`, `i12c` | Output | Signal | No
  `LineLumped`             | `v1a`, `v1b`, `v1c` | Input | Voltage signal | Yes
  `LineLumped`             | `v2a`, `v2b`, `v2c` | Input | Voltage signal | Yes
  `LineDistributed`        | `v1a`, `v1b`, `v1c`, `v2a`, `v2b`, `v2c` | Input | Voltage signal | Yes
  `LineDistributed`        | `i_ref1a`, `i_ref1b`, `i_ref1c`, `i_ref2a`, `i_ref2b`, `i_ref2c` | Output | Signal | No
  `LineDistributed`        | `i_inc1a`, `i_inc1b`, `i_inc1c`, `i_inc2a`, `i_inc2b`, `i_inc2c` | Output | Signal | No
  `LoadZ`                  | `va`, `vb`, `vc` | Input | Voltage signal | Yes
  `SexsPti`                | `va`, `vb`, `vc` | Input | Voltage signal | Yes
  `SexsPti`                | `vref`, `vs`, `vuel`, `voel` | Input | Signal | No
  `SexsPti`                | `efd` | Output | Signal | Yes
  `Tgov1`                  | `speed` | Input     | Signal        | No
  `Tgov1`                  | `pref`  | Input     | Signal        | No
  `Tgov1`                  | `pmech` | Output    | Signal        | Yes
  `Ieeet1`                 | `va`, `vb`, `vc` | Input | Voltage signal | Yes
  `Ieeet1`                 | `speed` | Input     | Signal        | No
  `Ieeet1`                 | `vref`  | Input     | Signal        | No
  `Ieeet1`                 | `vs`    | Input     | Signal        | No
  `Ieeet1`                 | `vuel`  | Input     | Signal        | No
  `Ieeet1`                 | `voel`  | Input     | Signal        | No
  `Ieeet1`                 | `efd`   | Output    | Signal        | Yes
  `Switch`                 | `v1a`, `v1b`, `v1c` | Input | Voltage signal | Yes
  `Switch`                 | `v2a`, `v2b`, `v2c` | Input | Voltage signal | Yes
  `Transformer`            | `v1a`, `v1b`, `v1c` | Input | Voltage signal | Yes
  `Transformer`            | `v2a`, `v2b`, `v2c` | Input | Voltage signal | Yes
  `Transformer`            | `i1a`, `i1b`, `i1c`, `i2a`, `i2b`, `i2c` | Output | Signal | No
  `GastPti`                | `speed`, `pref` | Input | Signal | No
  `GastPti`                | `pmech` | Output | Signal | Yes
  `Ieeest`                 | `input`, `speed` | Input | Signal | Exactly one
  `Ieeest`                 | `vct` | Input | Signal | When voltage cutout is enabled
  `Ieeest`                 | `output` | Output | Signal | Yes
  `Container`              | user-defined | Input/output | Public boundary | No

Declaring a signal creates a named connection. Its value is supplied by the
optional constant `value`, a component output, or the embedding program.
Declared constants can be updated by the application's `signal_step` events;
their finite values are preserved when connected controllers initialize.
Omit an optional input to use the model's internal default or latched value.

`IEEET1` is also accepted as the class name for `Ieeet1`. Its rated
line-to-line RMS voltage parameter `V` converts bus voltages to per unit;
connect its `efd` output to the Machine field-voltage input.

PWM, Converter, and Modulation vector ports use arrays of three scalar signal
IDs in phase order `a`, `b`, `c`. The scalar keys `ma`, `mb`, `mc`, `sa`, `sb`,
`sc`, and `voa`, `vob`, `voc` address individual phases. Vector monitors `m`,
`s`, and `vo` expand to these three scalar columns.

```json
{
  "class": "PWM",
  "id": "pwm",
  "params": { "M": 0.8, "fm": 60, "fc": 900 },
  "outputs": { "s": ["sa", "sb", "sc"] },
  "mon": ["s"]
}
```

To drive PWM from a controller, connect all three modulation inputs. Only
`fc` is required in this mode; `M` and `fm` apply to the unconnected sinusoidal
mode. The switching function uses the current modulation input and propagates
its derivatives through the connected signals.

```json
{
  "class": "PWM",
  "id": "pwm",
  "params": { "fc": 10000 },
  "inputs": { "m": ["ma", "mb", "mc"] },
  "outputs": { "s": ["sa", "sb", "sc"] },
  "mon": ["s"]
}
```

```json
{
  "class": "Converter",
  "id": "bridge",
  "inputs": { "s": ["sa", "sb", "sc"], "vdc": "dc", "i": ["ia", "ib", "ic"] },
  "outputs": { "vo": ["ea", "eb", "ec"], "idc": "idc" },
  "mon": ["vo", "idc"]
}
```

For a constant DC link, declare `{"id": "dc", "value": 1000.0}`.
Alternatively, the embedding program or another component supplies `dc`. A
DependentVoltageSource can consume `ea`, `eb`, and `ec` and publish its phase
currents to `ia`, `ib`, and `ic`. The bridge publishes the current drawn from
the DC link as `idc`, with `vdc * idc = vo · i`. Computed signals
are evaluated when read, including through Container boundaries. These
connections introduce no DAE variables.

For a dynamic DC link, declare `dc` and `idc` without constant values, and a
source-current signal such as `{"id": "isrc", "value": 80.0}`. Connect the capacitor
to the bridge above:

```json
{
  "class": "DCLink",
  "id": "capacitor",
  "params": { "C": 0.02 },
  "inputs": { "isrc": "isrc", "idc": "idc" },
  "outputs": { "vdc": "dc" },
  "mon": ["vdc", "isrc", "idc", "energy"]
}
```

The capacitor adds one differential voltage. Set its initial value with
`"capacitor": {"vdc": 600.0}` in the state file's `devices` object.

The current and voltage controllers use two-signal vector ports in power-invariant
`d`, `q` order. The same angle and angular frequency must be used for all
connected Park transforms and controller inputs. For example:

```json
{
  "class": "InnerCurrentControl",
  "id": "current",
  "params": { "L": 0.002, "Kp": 4.0, "Ki": 200.0, "Kaw": 2000.0, "Imax": 30.0, "Mmax": 0.95 },
  "inputs": { "v": ["vd", "vq"], "i": ["id", "iq"], "icmd": ["icmdd", "icmdq"], "omega": "omega", "vdc": "dc" },
  "outputs": { "ilim": ["ilimd", "ilimq"], "u": ["ud", "uq"] },
  "mon": ["xi", "ilim", "u"]
}
```

Use `Park` with `input` in `a`, `b`, `c` order and `out` in `d`, `q`, `0` order.
Set `params: {"inverse": true}` to reverse the transformation. Its scalar keys
are `u1`, `u2`, `u3` and `y1`, `y2`, `y3`. Connect `Angle.theta` to each Park
operator and supply the same `omega` to the Angle and controllers.

For the switching bridge, inverse-transform `[ud, uq, 0]`, connect the resulting
three-phase voltage command to `Modulation.u`, and connect `Modulation.m` to
`PWM.m`. Both the current controller and Modulation use the bridge's DC-link
voltage. Modulation requires a finite positive DC voltage.

For cascaded grid-forming control, connect `OuterVoltageControl.icmd` to the
current controller's `icmd`, and return the limited `ilim` to the voltage
controller. The `v` input is the filter-capacitor voltage; the current loop's
`i` input is the inverter-side current and the voltage loop's `ig` input is the
grid-side current. Vector monitors expand to scalar `d` and `q` columns.

#### Device classes

  Class                 | Model
  ----------------------|------------------------------------------------------
  `PWM`                 | [PWM](Component/Controller/PWM/README.md)
  `Park`                | [Park](Operators/Reference/Park/README.md)
  `PLL`                 | [PLL](Operators/Reference/PLL/README.md)
  `Angle`               | [Angle](Operators/Reference/Angle/README.md)
  `Modulation`          | [Modulation](Operators/Modulation/README.md)
  `InnerCurrentControl` | [InnerCurrentControl](Component/Controller/InnerCurrentControl/README.md)
  `OuterPowerControl`   | [OuterPowerControl](Component/Controller/OuterPowerControl/README.md)
  `OuterVoltageControl` | [OuterVoltageControl](Component/Controller/OuterVoltageControl/README.md)
  `DCLink`              | [DCLink](Component/Controller/DCLink/README.md)
  `Converter`           | [Converter](Operators/Converter/README.md)
  `Bus`                 | [Bus](Component/Bus/README.md)
  `DependentVoltageSource` | [DependentVoltageSource](Component/Source/DependentVoltageSource/README.md)
  `VoltageSource`       | [VoltageSource](Component/Source/VoltageSource/README.md)
  `LineLumped`          | [LineLumped](Component/Line/LineLumped/README.md)
  `LineDistributed`     | [LineDistributed](Component/Line/LineDistributed/README.md)
  `LoadZ`               | [LoadZ](Component/Load/LoadZ/README.md)
  `Switch`              | [Switch](Component/Switch/README.md)
  `Transformer`         | [Transformer](Component/Transformer/README.md)
  `Machine`              | [Machine](Component/Source/Machine/README.md)
  `REGFMA`, `Regfma`      | [REGFMA/REGFM_A1](Component/Source/REGFMA/README.md)
  `Tgov1`                | [TGOV1](Component/Controller/TGOV1/README.md)
  `Ieeet1`              | [IEEET1](Component/Controller/IEEET1/README.md)
  `SexsPti`             | [SEXS-PTI](Component/Controller/SEXS-PTI/README.md)
  `GastPti`             | [GASTPTI](Component/Controller/GASTPTI/README.md)
  `Ieeest`              | [IEEEST](Component/Controller/IEEEST/README.md)
  `Container`            | Recursive collection of devices and signals

#### Parameter values

Parameter types follow each model's declarations:

- Real parameters accept finite numeric values written as integer or floating
  literals. Boolean values are distinct from numbers.
- Index parameters require nonnegative integral values within the model's index
  range; fractional values and out-of-range conversions are rejected.
- Three-phase vectors require exactly three entries, and matrices require
  exactly three rows of three numeric entries.
- Unknown fields, parameter names, ports, monitored variables, and submodels are
  errors. Invalid input reports the component and offending field.

`SexsPti`, `SEXS-PTI`, and `SEXS` name the EMT [SEXS-PTI controller](Component/Controller/SEXS-PTI/README.md).
Its required `V` is the terminal line-to-line RMS voltage in volts; optional
`Tr` adds terminal-voltage measurement lag and defaults to zero.

`GastPti`, `GASTPTI`, and `GAST` name the EMT [GASTPTI controller](Component/Controller/GASTPTI/README.md).
Inputs `speed` and `pref` are optional, and output `pmech` is required.
`speed` is absolute rotor speed (one at synchronous). Required parameter
`S` is the connected machine rating in VA; `pref` and `pmech` use that base.
The optional `Trate` is in MW and defines the internal turbine base.

### IEEEST stabilizer

`Ieeest` and `IEEEST` select the complete PhasorDynamics stabilizer cascade
with EMT signal connections and compensated-voltage cutout. Connect
exactly one of `inputs.input` (generic per-unit signal) or `inputs.speed`
(absolute per-unit rotor speed), and connect `outputs.output` to the
exciter `vs` signal. Nonzero `Vcl` or `Vcu` requires `inputs.vct` in per
unit. Zero disables a threshold. The PSLF delay extension is unsupported;
`Tdelay` must be zero. See [IEEEST](Component/Controller/IEEEST/README.md)
for the parameter, bypass, initialization, and limiter contracts.
