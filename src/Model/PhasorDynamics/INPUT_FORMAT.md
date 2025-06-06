# Grid Dynamics case format specification

Adam Birchfield, Texas A&M University, 2/28/2025
Version 0.1

## Overview

This document describes a data format for grid dynamics cases intended to
be used in the SCIDAC-OE project "Next-Generation Grid Simulations". The
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

-   Non-repetitive to keep file size concise and enhance human
    readibility. For example, if we have 300 GENROU devices, we do not
    need to say "Tdopp" for every one of them. Just putting "GENROU" for
    each lets us know what the parameters mean and their order.

-   Backward compatibility. Newer software should always be able to read
    older files. When possible, also have forward compatibility where
    older software can read newer files as long as they do not contain
    newly added or modified devices.

## Format

The root element in the format is an object with three keys: `header`,
`nodes`, and `devices`. `header` contains information about the case,
`nodes` is an array of nodes, and `devices` is an array of system devices.

### Header

Contained in the `header` key is an object with the following items:

   Name              | Value
 --------------------|-------------------------------------------------------
  `format_version`   | Non-negative integer indicating the format version
  `format_revision`  | Non-negative integer indicating the format revision
  `case_name`        | A string containing the name of the case
  `case_description` | A string describing what the case models
  `case_comments`    | A string containing any additional notes about the case
  `freq_base`        | A floating point value indicating the system frequency base in hertz (Hz). This is commonly 60 Hz
  `va_base`          | A floating point value indicating the system power base in volt-amperes (VA). This is commonly 100e6 VA

### Nodes

Contained in the `nodes` key is an array of objects, each of which represent
a node and have the following fields:

  Name               | Description
  -------------------|------------------------------------------------------
  `id`               | Unique positive (> 0) integer identifying the node
  `class`            | A string indicating the class of node. See the table below for more information
  `name`             | A string containing the name of the node. This may be empty or non-unique
  `default_voltages` | An array of floating point values specifying default voltages or signal values. The length and meaning of this field is dependent upon the node class. See the table below for more information
  `voltage_base`     | A floating point value functioning as a scaling factor for the voltages in volts (V)
  `freq_base`        | Optional field to override the system frequency base
  `va_base`          | Optional field to override the system power base
  `monitor`          | Optional field, which is an array specifying variables to record the value of in an output channel. Available variables are determined by the node class
  `extension`        | Optional field containing an object with implementation-defined keys

#### Node classes

As of the current version and revision, the following node classes are
specified:

  Node class         | Description                                                | Available variables
  -------------------|------------------------------------------------------------|----------------------
  `bus`              | Positive-sequence, AC phasor domain bus                    | `["Vr", "Vi"]`
  `infinite_bus`     | Positive-sequence, AC phasor domain bus with fixed voltage | `["Vr", "Vi"]`
  `emt_bus`          | 3-phase bus with instantaneous voltages                    | `["Va", "Vb", "Vc"]`
  `infinite_emt_bus` | 3-phase bus with instantaneous voltages                    | `["Va", "Vb", "Vc"]`
  `control`          | A single control signal                                    | `["x"]`

This list is subject to change.

### Devices

Contained in the `devices` section is an array of objects, each of which
represent a device and have the following fields:

  Name              | Description
  ------------------|------------------------------------------------------
  `class`           | A string indicating the class of device. See the table below for more information
  `connected_nodes` | An array of node ids indicating the nodes connected to the device. The length and meaning of the indices of this array is dependent upon the device class. An index may be `0` to indicate that no node is connected to the device at that index
  `disambiguator`   | A string disambiguating the device from others. Each device must have a unique combination of connected nodes and this string. It is recommended that this string is within the range of 1-2 characters (inclusive)
  `init_parameters` | An array of initialization parameters. The contents of this field are dependent upon the device class. See the table below for more information
  `freq_base`       | Optional field to override the system frequency base
  `va_base`         | Optional field to override the system power base
  `monitor`         | Optional field, which is an array specifying variables to record the value of in an output channel. Available variables are determined by the device class
  `extension`       | Optional field containing an object with implementation-defined keys

#### Device classes

As of the current version and revision, the following device classes
are specified:

  Device class  | Description                                          | Nodes                                         | Initialization parameters
  --------------|------------------------------------------------------|-----------------------------------------------|----------------------------
  `branch`      | a basic algebraic pi model for a line or transformer | 2: `Bus1`, `Bus2`                             | 4: `R`, `X`, `G`, `B`
  `static_load` | a basic static ZIP load                              | 1: `Bus`                                      | 6: `Pz`, `Qz`, `Pi`, `Qi`, `Pp`, `Qp`
  `GENROU`      | 6th order machine model                              | 3: `Bus`, `exciter_signal`, `governor_signal` | 18: `p0`, `q0`, `H`, `D`, `Ra`, `Tdop`, `Tdopp`, `Tqopp`, `Tqop`, `Xd`, `Xdp`, `Xdpp`, `Xq`, `Xqp`, `Xqpp`, `Xl`, `S10`, `S12`
  `bus_fault`   | simple impedance-based fault at a bus                | 2: `Bus`, `control signal`                    | 3: `state0`, `R`, `X`

NOTE: we should add the variables above too as done for the nodes

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
    "nodes": [
        {
            "id": 1,
            "class": "bus",
            "name": "Bus 1",
            "default_voltages": [0.994988, 0.099997],
            "voltage_base": 115e3,
            "monitor": ["Vr", "Vi"]
        },
        {
            "id": 2,
            "class": "infinite_bus",
            "name": "Bus 2",
            "default_voltages": [1.0, 0],
            "voltage_base": 115e3
        }
    ],
    "devices": [
        {
            "class": "branch",
            "connected_nodes": [1, 2],
            "disambiguator": "1",
            "init_parameters": [0, 0.1, 0, 0]
        },
        {
            "class": "GENROU",
            "connected_nodes": [1, 0, 0],
            "disambiguator": "1",
            "init_parameters": [1, 0.05013, 3, 0, 0, 7, 0.04, 0.05, 0.75, 2.1, 0.2, 0.18, 0.5, 0.5, 0.18, 0.15, 0, 0],
            "monitor": ["delta", "omega"]
        }
        {
            "class": "bus_fault",
            "connected_nodes": [1],
            "disambiguator": "1",
            "init_parameters": [0, 0, 1e-3]
        }
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
