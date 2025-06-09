# Grid Dynamics case format specification

Adam Birchfield, Texas A&M University, 2/28/2025
Version 0.1

## Overview

This document describes a data format for grid dynamics cases intended to
be used in the SciDAC-OE project "Next-Generation Grid Simulations". The
format is designed first for implementation as UTF-8 encoded JSON but may
also be encoded as [MessagePack](https://msgpack.org).

### Goals

-   Simplicity both of user understanding and of software
    implementation.

-   Flexible to handle various kinds of power system dynamics models
    including phasor dynamics (PD), electromechanical transients (EMT)
    and hybrid models of the two.

-   Conformity, as much as possible, to the style and formulations of
    commercial software formats.

-   Extensible to be able to add many different kinds of nodes and
    device classes.

-   Self-describing, by explicitly labeling parameters.

-   Backward compatibility. Newer software should always be able to read
    older files. When possible, also have forward compatibility where
    older software can read newer files as long as they do not contain
    newly added or modified devices.

## Format

The root element in the format is an object with three keys: `header`,
`buses`, and `devices`. `header` contains information about the case,
`buses` is an array of buses (which can include electrical buses of 
different types as well as signal buses), and `devices` is an array 
of system devices.

### Header

Contained in the `header` key is an object with the following items:

   Name              | Value
 --------------------|-------------------------------------------------------
  `format_version`   | Non-negative integer indicating the format version
  `format_revision`  | Non-negative integer indicating the format revision
  `case_name`        | A string containing the name of the case
  `case_date`        | A string in standard date/time format
  `case_description` | A string with more specific description of what is modeled in the case
  `case_comments`    | A string with additional notes as needed
  `freq_base`        | A floating point value indicating the system frequency base in hertz (Hz). This is commonly 60 Hz
  `va_base`          | A floating point value indicating the system power base in volt-amperes (VA). This is commonly 100e6 VA

### Buses

Contained in the `buses` key is an array of objects, each of which represent
a bus and has the following fields:

  Name               | Description
  -------------------|------------------------------------------------------
  `number`           | Unique positive (> 0) integer identifying the node
  `class`            | A string indicating the class of node. See the table below for more information
  `name`             | Optional string containing the name of the node. This may be empty or non-unique
  `init`             | Optional object mapping string variable names to floating point values, specifying default voltages or signal values. The available initialization variables are dependent upon the node class. Any variables missing will be given default values. If this object is missing, all variables will be given default values. See the table below for more information
  `vbase`            | Optional floating point value giving the voltage base in volts (V). If omitted, default value of 1 V is assumed (common for signal buses)
  `mon`              | Optional field, which is an array specifying variables to monitor the value of in an output channel. Available variables include all the initialization variables, along with others as determined by the node class. See the table below for more information
  `freq_base`        | Optional uncommon field to override the system frequency base at this bus
  `va_base`          | Optional uncommon field to override the system power base at this bus
  `extension`        | Optional field containing an object with implementation-defined keys

#### Bus classes

As of the current version and revision, the following bus classes are
specified:

  Bus class          | Description                                                | Initialization variables | Other variables available to monitor
  -------------------|------------------------------------------------------------|------------------------- | -------------------------
  `bus`              | Positive-sequence, AC phasor domain bus                    | `Vr`, `Vi`               | `Vm`, `Va`
  `infinite_bus`     | Positive-sequence, AC phasor domain bus with fixed voltage | `Vr`, `Vi`               | `Vm`, `Va`
  `emt_bus`          | 3-phase bus with instantaneous voltages                    | `Va`, `Vb`, `Vc`         | 
  `infinite_emt_bus` | 3-phase bus with instantaneous voltages                    | `Va`, `Vb`, `Vc`         |
  `control`          | A single control signal                                    | `x`                      |

This list is subject to change.

### Devices

Contained in the `devices` section is an array of objects, each of which
represent a device and has the following fields:

  Name              | Description
  ------------------|------------------------------------------------------
  `class`           | A string indicating the class of device. See the table below for more information
  `ports`           | An object mapping the object's port names (depending on the device class as specified in the table below) to the associated bus number to which it is connected. Any field listed under variables available to monitor can also be added here as a read-only port
  `id`              | A string disambiguating the device from others. Each device in a class must have a unique combination of required port bus numbers and this string. This string should be 1 or 2 characters long.
  `params`          | An object mapping initialization parameters to numerical values, depending on the class. See the table below for more information
  `mon`             | Optional field, which is an array specifying variables to record the value of in an output channel. Available variables are determined by the device class, as specified in the table below
  `va_base`         | Optional field to override the system power base for this device
  `freq_base`       | Optional uncommon field to override the system frequency base for this device
  `extension`       | Optional field containing an object with implementation-defined keys

#### Device classes

As of the current version and revision, the following device classes
are specified:

  Device class  | Description                                          | Ports                                         | Initialization parameters | Variables available to monitor
  --------------|------------------------------------------------------|-----------------------------------------------|---------------------------- | -------------------------
  `branch`      | a basic algebraic pi model for a line or transformer | `bus1`, `bus2`                             | `R`, `X`, `G`, `B` | `ir1`, `ii1`, `im1`, `p1`, `q1`, `ir2`, `ii2`, `im2`, `p2`, `q2`
  `static_load` | a basic static ZIP load                              | `bus`                                      | `Pz`, `Qz`, `Pi`, `Qi`, `Pp`, `Qp` | `ir`, `ii`, `p`, `q`
  `GENROU`      | 6th order machine model                              | `bus`, `exciter_signal`\*, `governor_signal`\* | `p0`, `q0`, `H`, `D`, `Ra`, `Tdop`, `Tdopp`, `Tqopp`, `Tqop`, `Xd`, `Xdp`, `Xdpp`, `Xq`, `Xqp`, `Xqpp`, `Xl`, `S10`, `S12` | `ir`, `ii`, `p`, `q`, `delta`, `omega`
  `bus_fault`   | simple impedance-based fault at a bus                | `bus`, `control_signal`\*                    | `state0`, `R`, `X` | `state`, `ir`, `ii`

Ports marked with \* are optional and, if missing, will be assumed to be connected to a constant value.


This list is subject to change.

## Example File for a 2-Bus System

```json
{
    "header": {
        "format_version": 0,
        "format_revision": 1,
        "case_name": "Two-bus test case 1",
        "case_description": "A two-bus test case for demonstrating the dynamics format",
        "comments": "This case is set up to monitor the voltage at both buses and the machine angle and speed",
        "freq_base": 60,
        "va_base": 100e6
    },
    "buses": [
        { "number": 1, "class": "bus", "name": "Bus 1", "init": {"Vr":0.994988, "Vi":0.099997}, "vbase": 115e3, "mon": ["Vr", "Vi"] },
        { "number": 2, "class": "infinite_bus", "name": "Bus 2", "init": {"Vr":1.0, "Vi":0}, "vbase": 115e3 }
    ],
    "devices": [
        { "class": "branch", "ports": {"bus1":1, "bus2":2}, "id": "1", "params": {"R":0, "X":0.1, "G":0, "B":0} },
        { "class": "GENROU", "ports": {"bus":1}, "id": "1", "params": {"p0":1, "q0":0.05013, "H":3, "D":0, "Ra":0, "Tdop":7, "Tdopp":0.04, "Tqopp":0.05, 
               "Tqop":0.75, "Xd":2.1, "Xdp":0.2, "Xdpp":0.18, "Xq":0.5, "Xqpp":0.18, "Xl":0.15, "S10":0, "S12":0}, "mon": ["delta", "omega"] }
        { "class": "bus_fault", "ports": {"bus":1}, "id": "1", "params": {"state0":0, "R":0, "X":1e-3} }
    ]
}
```

## Appendix: Output file format

There could be multiple output file formats, but the simplest output
file format is a comma-separated text file table. The first row is a
list of headers for the monitor channels. Then each subsequent row
represents a time point and gives the value for each channel. The first
channel is always "Time [s]". The second channel is always "Solver
Status", an integer where 0 is success, positive integers are warnings
and negative integers are terminal errors. Then the other channels come
from any nodes or devices which have the "monitor" property specified in
their extra data, with one channel for each item in the monitor array.
The channels must be in the order they appear in the input file. In the
example above, 4 channels are identified: "Vr" and "Vi" for Bus 1 and
"delta" and "omega" for the GENROU device at Bus 1. A sample output file
might look like this:

```
Time [s], Solver Status, Bus 1 (1) Vr, Bus 1 (1) Vi, GENROU Bus 1 (1-1) delta, GENROU Bus 1 (1-1) omega

0, 0, 0.994988, 0.099997, 0.553983, 0
0.00416667, 0, 0.994988, 0.099999, 0.553983, 6.84E-09
0.00833333, 0, 0.994988, 0.0999991, 0.553983, 1.32E-08
0.0125, 0, 0.994988, 0.0999992, 0.553984, 1.92E-08
0.0166667, 0, 0.994988, .0999992, 0.553984, 2.49E-08
0.0208333, 0, 0.994988, 0.0999993, 0.553984, 3.02E-08
0.025, 0, 0.994988, 0.0999993, 0.553984, 3.51E-08
```
