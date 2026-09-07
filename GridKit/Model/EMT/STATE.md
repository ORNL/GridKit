# EMT state format specification

## Overview

This document describes the JSON data format for EMT operating points. A
case file carries parameters and topology only; a state file carries the
operating point. Fields use each component's existing `Outputs` enum.
Missing records or null fields use model defaults. Every supplied component
path and output name must exist in the case, including those with null values.
Unknown fields and duplicate component paths across sections are errors.

The state carries instantaneous SI values at the initialization instant.
Synthesizing instantaneous values from an RMS or phasor description is an
upstream tool concern and is outside the EMT model specification.

## Format

The root object may contain `header`, `buses`, and `devices`.

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
  `i12a`, `i12b`, `i12c` | Optional instantaneous `LineLumped` or `Switch` series currents from terminal 1 to terminal 2, in amperes
  `ia`, `ib`, `ic` | Optional instantaneous `Machine`, `LoadZ`, `VoltageSource`, or `DependentVoltageSource` current injections into the bus, in amperes

Other outputs use their model output names and units. All output values must
be finite; missing or null values use model defaults.

## Application

The application reads the state file and passes its values to
`SystemModel::initialize(state)` after allocation. The integrator then solves
consistent algebraic variables and derivatives.

The system validates component paths and output names, resolves initialization
dependencies, and reconciles machine operating-point requirements with
prescribed controller outputs before changing state. Conflicting requirements
or constants are rejected; components initialize only their own variables.
