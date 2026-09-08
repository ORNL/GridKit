# EMT state format specification

## Overview

This document describes the JSON data format for EMT operating points. A
case file carries parameters and topology only; a state file carries the
operating point. Fields use each model's documented initial-state names.
Missing records or null fields use model defaults. Every supplied component
path and output name must exist in the case, including those with null values.
Unknown fields and duplicate component paths across sections are errors.

State fields use each model's documented quantities and units. Phase voltages
and currents are instantaneous SI values. Initialization helpers may derive
them from RMS or phasor data before applying the state.

## Format

The root object may contain `header`, `buses`, `devices`, and `history`.

### Header

Contained in the optional `header` key is an object with the following items:

   Name          | Value
 ----------------|-------------------------------------------------------
  `version`      | Optional integer identifying the state format version
  `time`         | Optional floating-point model time of the state
  `created`      | Optional string with the wall-clock creation time
  `description`  | Optional string describing the state

The application initializes at time zero. A supplied `time` must be finite and
zero; a nonzero time cannot be used as a restart time. Null metadata fields
are treated as omitted.

### Buses

Contained in the `buses` key is an object mapping Bus component paths from the
case file to bus states. A root Bus uses its local ID; a nested Bus uses its
dot-qualified path, for example `left.network.bus`:

   Name  | Value
 --------|-------------------------------------------------------
  `va`   | Optional instantaneous phase a voltage in volts
  `vb`   | Optional instantaneous phase b voltage in volts
  `vc`   | Optional instantaneous phase c voltage in volts

### Devices

Contained in the `devices` key is an object mapping component paths from the
case file to device states. Nested paths are qualified by each containing
Container, for example `plant.machine`:

   Name   | Value
 ---------|-------------------------------------------------------
  `open`  | Optional Boolean switch command, true is open
  `theta`, `omega` | Optional `PLL` angle in radians and angular frequency in radians per second; defaults are the inferred voltage angle and nominal frequency
  `theta` | Optional `Angle` reference angle in radians, default zero
  `ud`, `uq`, `ilimd`, `ilimq` | Optional `InnerCurrentControl` outputs; voltage commands determine the integral, limited currents must match the limiter
  `icmdd`, `icmdq` | Optional `OuterVoltageControl` current-command outputs in amperes; determine the integral contributions
  `icmdd`, `icmdq` | Optional `OuterPowerControl` current commands in amperes; defaults give zero integral contribution
  `i12a`, `i12b`, `i12c` | Optional instantaneous `LineLumped` or `Switch` series currents from terminal 1 to terminal 2, in amperes
  `i12a`, `i12b`, `i12c`, `psi1a`, `psi1b`, `psi1c`, `psi2a`, `psi2b`, `psi2c` | Optional `Transformer` series leakage current and magnetizing flux linkages in per unit, default zero
  `ia`, `ib`, `ic` | Optional instantaneous `Machine`, `LoadZ`, `VoltageSource`, or `DependentVoltageSource` current injections into the bus, in amperes

Other outputs use their model output names and units. All state values must
be finite; missing or null values use model defaults. The current and voltage
controllers preserve their supplied integral states; consistent initialization
resolves their derivatives. Their computed control outputs cannot be prescribed
in the state file.

### History

Each `LineDistributed` requires a record in `history`, keyed by its qualified
component path. All five fields are required:

Name | Units | Value
---- | ----- | -----
`omega` | [rad/s] | Finite nonnegative prehistory angular frequency
`i_ref1`, `i_ref2` | [A] | Three instantaneous reflected-current values per terminal
`d_i_ref1`, `d_i_ref2` | [A/s] | Three reflected-current derivatives per terminal

For `omega: 0`, both derivatives must be zero and the prehistory is constant.
For positive frequency the values and derivatives define the harmonic
prehistory described by [LineDistributed](Component/Line/LineDistributed/README.md#initialization).
These are explicit time-domain history data, separate from topology and fit
coefficients. A zero history is suitable for a line energized at time zero:

```json
"history": {
  "tie": {
    "omega": 0,
    "i_ref1": [0, 0, 0], "i_ref2": [0, 0, 0],
    "d_i_ref1": [0, 0, 0], "d_i_ref2": [0, 0, 0]
  }
}
```

## Application

The application reads the state file and passes its values to
`SystemModel::initialize(state)` after allocation. The integrator then solves
consistent algebraic variables and derivatives.

The system validates component paths and output names, resolves initialization
dependencies, and reconciles machine operating-point requirements with
prescribed controller outputs before changing state. Conflicting requirements
or constants are rejected; components initialize only their own variables.
