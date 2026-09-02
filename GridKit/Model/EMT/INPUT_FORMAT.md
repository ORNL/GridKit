# EMT case format specification

## Overview

This document describes the JSON data format for EMT cases. The format
follows the grid dynamics case format described in
[PhasorDynamics INPUT_FORMAT](../PhasorDynamics/INPUT_FORMAT.md), specialized
to instantaneous phase-coordinate models: every quantity is in SI units, and
there is no system-wide power or frequency base.

## Format

The root object contains `header`, `buses`, and `devices`, and may also
contain `signals` and `monitors`. `header` contains information about the
case, `buses` is an array of electrical buses, and `devices` is an array of
system devices.

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
of a [bus](#buses) or [device](#devices)). The following fields are
supported:

  Name               | Description
  -------------------|------------------------------------------------------
  `file_name`        | Optional string indicating output file name. If omitted, `stdout` is used.
  `format`           | One of { "CSV", "JSON", "YAML" } (case-insensitive)
  `delim`            | Optional string specifying delimiter to use for CSV output (default is `","`).

### Buses

Contained in the `buses` key is an array of objects, each of which represents
a [Bus](Bus/README.md) and has the following fields:

  Name               | Description
  -------------------|------------------------------------------------------
  `number`           | Unique non-negative integer identifying the bus
  `class`            | The string `"Bus"`
  `name`             | String containing the name of the bus
  `init`             | Optional object with the initial instantaneous phase voltages `va`, `vb`, and `vc` in volts. Missing entries default to zero
  `params`           | Optional object mapping bus parameter names to values
  `mon`              | Optional array of variables to monitor, from the bus model's Monitors table

### Signals

Contained in the optional `signals` key is an array of objects, each of which
represents a scalar signal node connecting device signal ports:

  Name               | Description
  -------------------|------------------------------------------------------
  `name`             | String containing the name of the signal node
  `signal_id`        | Unique non-negative integer identifying the signal node

### Devices

Contained in the `devices` key is an array of objects, each of which
represents a device and has the following fields:

  Name               | Description
  -------------------|------------------------------------------------------
  `class`            | A string indicating the device class. See the table below
  `id`               | A string disambiguating this device
  `params`           | Object mapping parameter names to values, using the JSON names from the device model's Model Parameters table
  `ports`            | Object mapping port names to bus numbers and signal ids, using the names from the device model's Model Parameters and Model Ports tables
  `mon`              | Optional array of variables to monitor, from the device model's Monitors table

#### Device classes

  Class                 | Model
  ----------------------|------------------------------------------------------
  `VoltageSource`       | [VoltageSource](Component/Source/VoltageSource/README.md)
  `LineLumped`          | [LineLumped](Component/Line/LineLumped/README.md)
  `LoadZ`               | [LoadZ](Component/Load/LoadZ/README.md)
  `Switch`              | [Switch](Component/Switch/README.md)
  `SynchronousMachine`  | [SynchronousMachine](Component/Machine/SynchronousMachine/README.md)
  `Tgov1`               | [Tgov1](Governor/Tgov1/README.md)

#### Parameter values

Parameter values are typed by their JSON representation:

- A boolean, a number with a decimal point, or an integer maps to a Boolean,
  real, or integer parameter respectively; real-valued scalar parameters must
  be written with a decimal point.
- A length-3 array maps to a three-phase vector parameter, integer-valued
  when every element is an integer.
- A 3x3 nested array maps to a real three-phase matrix parameter.
