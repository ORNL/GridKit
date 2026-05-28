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
  `case_date_time`   | Optional string in the ISO 8601 format indicating a datetime associated with the case. Including a time is optional, but if one is specified, it is recommended that a UTC offset be included.
  `case_description` | A string with more specific description of what is modeled in the case
  `case_comments`    | A string with additional notes as needed
  `freq_base`        | A floating point value indicating the system frequency base in hertz (Hz). This is commonly 60 Hz
  `va_base`          | A floating point value indicating the system power base in volt-amperes (VA). This is commonly 100e6 VA

### Monitors

Contained in the `monitors` key is an array of objects, each of which describes
an output for monitored variables (those listed in the `mon` field of a
[bus](#buses) or [device](#devices). The following fields are supported:

  Name               | Description
  -------------------|------------------------------------------------------
  `file_name`        | Optional string indicating output file name. If omitted, `stdout` is used.
  `format`           | One of { "CSV", "JSON", "YAML" } (case-insensitive)
  `delim`            | Optional string specifying delimiter to use for CSV output (default is `","`).

__NOTE__: If `monitors` entry is omitted entirely, you will get the default CSV
output to `stdout`. If you wish to change the console output format, use an
entry here without the `file_name` field. For example:
```json
  "monitors": [ { "format": "YAML" } ]
```

### Buses

Contained in the `buses` key is an array of objects, each of which represent
a bus and has the following fields:

  Name               | Description
  -------------------|------------------------------------------------------
  `number`           | Unique positive (> 0) integer identifying the node
  `class`            | A string indicating the class of node. See the table below for more information
  `name`             | Optional string containing the name of the node. This may be empty or non-unique
  `init`             | Optional object mapping string variable names to floating point values, specifying default voltages or signal values. The available initialization variables are dependent upon the node class. Any variables missing will be given default values, which are specified beneath the table below. If this object is missing, all variables will be given default values. See the table below for more information
  `v_base`           | Optional floating point value giving the voltage base in volts (V).
  `mon`              | Optional field, which is an array specifying variables to monitor the value of in an output channel. Available variables include all the initialization variables, along with others as determined by the node class. See the table below for more information
  `extension`        | Optional field containing an object with implementation-defined keys

#### Bus classes

As of the current version and revision, the following bus classes are
specified:

  Bus class          | Description                                                | Initialization variables | Other variables available to monitor
  -------------------|------------------------------------------------------------|------------------------- | -------------------------
  `bus`              | Positive-sequence, AC phasor domain bus                    | `Vr`, `Vi`               | `Vm`, `Va`
  `infinite_bus`     | Positive-sequence, AC phasor domain bus with fixed voltage | `Vr`, `Vi`               | `Vm`, `Va`
  `emt_bus`          | 3-phase bus with instantaneous voltages                    | `Va`, `Vb`, `Vc`         |

For fields named `Vr` or `Va`, the default value is `1.0`, otherwise it is
`0.0`. This list is subject to change.

### Signals

Contained in the `signals` key is an array of objects, each of which represent
a signal and has the following fields:

  Name               | Description
  -------------------|------------------------------------------------------
  `signal_id`        | Unique positive (> 0) integer identifying the signal
  `name`             | Optional string containing the name or description of the signal. This may be empty or non-unique

All interactions with signals are handled by `devices`, which reference signals via `signal_id`


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
  `extension`       | Optional field containing an object with implementation-defined keys

For more information, see the detailed documentation for each device class containing equations, electrical bus configuration, and signal inlets/outlets specifications.

#### Device classes

As of the current version and revision, the following device classes
are specified:

  Device class         | Description                                          | Ports                            | Initialization parameters   | Variables available to monitor
  ---------------------|------------------------------------------------------|----------------------------------|---------------------------- | -------------------------
  `Branch`             | algebraic pi model for a line or off-nominal transformer branch | `bus1`, `bus2`                   | `R`, `X`, `G`, `B`, `tap`, `phase` | `ir1`, `ii1`, `im1`, `p1`, `q1`, `ir2`, `ii2`, `im2`, `p2`, `q2`
  `Load`               | a basic static impedence load model                  | `bus`                            | `R`, `X` | `p`, `q`
  `Genrou`             | 6th order machine model                              | `bus`, `pmech`\*, `speed`\*, `efd`\*    | `p0`, `q0`, `H`, `D`, `Ra`, `Tdop`, `Tdopp`, `Tqop`, `Tqopp`, `Xd`, `Xdp`, `Xdpp`, `Xq`, `Xqp`, `Xqpp`, `Xl`, `S10`, `S12`, `mva`  | `ir`, `ii`, `p`, `q`, `delta`, `omega`, `speed`
  `Gensal`             | 5th order salient-pole machine model                 | `bus`, `pmech`\*, `speed`\*, `efd`\*    | `p0`, `q0`, `H`, `D`, `Ra`, `Tdop`, `Tdopp`, `Tqopp`, `Xd`, `Xdp`, `Xdpp`, `Xq`, `Xl`, `S10`, `S12`, `mva`  | `ir`, `ii`, `p`, `q`, `delta`, `omega`, `speed`, `Eqp`, `psidp`, `psiqpp`, `psidpp`, `vd`, `vq`, `te`, `id`, `iq`
  `GenClassical`       | the classical machine model                          | `bus`, `pmech`\*, `speed`\*, `efd`\*  | `p0`, `q0`, `H`, `D`, `Ra`, `Xdp`, `mva` | `ir`, `ii`, `p`, `q`, `delta`, `omega`
  `Tgov1 `             | the TGOV1 governor model                             | `pmech`, `speed`                 | `R`, `T1`, `T2`, `T3`, `Pvmax`, `Pvmin`, `Dt` | `none`
  `GastPti`            | the GASTPTI gas turbine-governor model               | `pmech`, `speed`\*, `pref`\*   | `R`, `T1`, `T2`, `T3`, `At`, `Kt`, `Vmax`, `Vmin`, `Dturb`, `Trate` | `pmech`, `fuelvalve`, `fuelflow`, `exhausttemp`, `vload`, `vtemp`, `vlv`
  `Ieeet1`             | the IEEET1 exciter model                             | `bus`, `speed`, `efd`, `vs`\*    | `Tr`, `Ka`, `Ta`, `Ke`, `Te`, `Kf`, `Tf`, `Vrmin`, `Vrmax`, `E1`, `E2`, `Se1`, `Se2`, `Ispdlim` | `efd`, `ksat`
  `SexsPti`            | the SEXS-PTI simplified exciter model                | `bus`, `efd`, `vs`\*             | `Ta`, `Tb`, `Te`, `K`, `Efdmax`, `Efdmin` | `efd`
  `Ieeest`      | the IEEEST stabilizer model                          | `input`, `output`                | `A1`, `A2`, `A3`, `A4`, `A5`, `A6`, `T1`, `T2`, `T3`, `T4`, `T5`, `T6`, `Ks`, `Lsmin`, `Lsmax`, `Vcl`, `Vcu`, `Tdelay` | `vss`
  `BusFault`           | simple impedance-based fault at a bus                | `bus`, `status`\*                | `state0`, `R`, `X` | `state`, `ir`, `ii`
  `BusToSignalAdapter` | signal adapter component for a bus                   | `bus`, `vr`, `vi`, `ir`, `ii`    |                             |

Ports marked with \* are optional and, if missing, will be assumed to be
connected to a constant value. This list is subject to change.


For `Branch`, `tap` and `phase` are optional parameters. If omitted, `tap`
defaults to `1.0` and `phase` defaults to `0.0` radians. Bus `bus1` is the tap
side for off-nominal transformer branches.

## Example File for a 2-Bus System

```json
{
   "header": {
       "format_version": 0,
       "format_revision": 1,
       "case_name": "Two-bus test case 1",
       "case_description": "A two-bus test case for demonstrating the dynamics format",
       "case_comments": "This case is set up to monitor the voltage at both buses and the machine angle and speed",
       "freq_base": 60.0,
       "va_base": 100e6
   },
   "buses": [
       { "number": 1, "class": "bus", "name": "Bus 1", "init": {"Vr":0.994988, "Vi":0.099997}, "v_base": 115e3, "mon": ["Vr", "Vi"] },
       { "number": 2, "class": "infinite_bus", "name": "Bus 2", "init": {"Vr":1.0, "Vi":0.0}, "v_base": 115e3 }
   ],
   "signals": [
       { "signal_id": 1, "name": "Machine Speed Deviation"},
       { "signal_id": 2, "name": "Mechanical Power"},
       { "signal_id": 3, "name": "Excitation Field"}
   ],
   "devices": [
       { "class": "Branch", "ports": {"bus1":1, "bus2":2}, "id": "BR1", "params": {"R":0.0, "X":0.1, "G":0.0, "B":0.0} },
       { "class": "Genrou", "ports": {"bus":1, "speed": 1, "pmech":2, "efd":3}, "id": "DV1", "params": {"p0":1.0, "q0":0.05013, "H":3.0, "D":0.0, "Ra":0.0, "Tdop":7.0, "Tdopp":0.04, "Tqopp":0.05, "Tqop":0.75, "Xd":2.1, "Xdp":0.2, "Xdpp":0.18, "Xq":0.5, "Xqp": 0.0, "Xqpp":0.18, "Xl":0.15, "S10":0.0, "S12":0.0}, "mon": ["delta", "omega"] },
       { "class": "Tgov1", "ports": {"bus":1, "speed": 1, "pmech":2}, "id": "DV2", "params": {"R":0.05, "T1":0.5,"T2":2.5, "T3":7.5, "Pvmax":0, "Pvmin":1, "Dt":0}},
       { "class": "Ieeet1", "ports": {"bus":1, "speed": 1, "efd":3}, "id": "DV3", "params": {"Tr":0.001, "Ka":50.0, "Ta":0.04, "Ke":-0.06, "Te":0.6, "Kf":0.09, "Tf":1.46, "Vrmin":-1, "Vrmax":1, "E1":2.8, "E2":3.373, "Se1":0.04, "Se2":0.33, "Ispdlim":0}},
       { "class": "BusFault", "ports": {"bus":1}, "id": "EVT1", "params": {"state0": false, "R":0.0, "X":1e-3} }
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
