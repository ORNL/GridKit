## Grid Dynamics Case Format

Adam Birchfield, Texas A&M University, 2/28/2025
Version 0.1

### Overview

This document describes the data format for grid dynamics cases which
will be used in the SCIDAC-OE project "Next-Generation Grid Simulations"
and can also be used to inform future formats and software efforts.

#### Overall goals of the format

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

#### Overview of the input format

-   Using JSON in UTF-8 as the base format style

-   Root is an object with three main sections as key-value pair items:

    -   "header" section is itself an object with various system-wide
        information about the case as key-value pairs: name, date, base
        frequency, etc.

    -   "nodes" section is an array of nodes: primarily electrical buses
        and control signals

    -   "devices" section is an array of system devices, including
        generators, loads, branches, controllers, and events

```json
{
    "header":{... },
    "nodes":[... ],
    "devices":[... ]
}
```


### Header Section

```json
"header": { "key": value, ...}
```

The header section is an object, with the following standardized items,
allowing for extensibility as well.

   Name             | Value                                                  
 -------------------|---------------------------------------------------------
  format_version    | Integer of main format version (0 for this version 0.1)
  format_revision   | Integer of format revision (1 for this version 0.1)
  case_name         | Short string for the case name 
  case_date         | String in a standard date/time format 
  case_description  | Longer string giving a more specific description of what is modeled in the case 
  comments          | Even longer string with notes as needed 
  freq_base         | System frequency base, floating point in Hz (commonly 60)
  va_base           | System power base, floating point in VA (commonly 100e6)

### Nodes Section

```json
"nodes": [[..., ], [..., ], [..., ]]
```

The nodes section is an array of arrays. Each array represents one node,
and has exactly 6 items as follows.

  Index |  Node item
  ------|---------------------------------------------------------
  0     | Number, a unique integer \> 0 to identify the node 
  1     |  Node class, string, for what type of node it is. See table below for supported classes 
  2     | String name of the node, not necessarily unique. If desired, for simplicity it can be an empty string. Nodes are identified by number, not name, but this can be a helpful label for debugging.
  3     |  Array of default voltages or signal values. Length of array and meaning depends on node class (see table below).
  4     | Voltage base for per-unit, floating point value in V. Set to 1 to indicated that voltages or signals given are already in actual units.
  5     |  Extra, object with key-value pairs for any extra information to provide with the node. This can be an empty object. If you want to specify a non-default frequency or power base, do it here with `"freq_base"` or `"va_base"`. If you add the item `"monitor":[]`, the voltage or signal values specified in the array for this node will be added to an output channel. Other uses of this section would be coordinates for one-line diagram drawing, membership in components or areas, etc.
  

### Table of Supported Node Classes

For now the following node classes are specified (though not all
implemented yet)

  Node Class   | Description                                 | Variables
  -------------|---------------------------------------------|--------------
  bus          | positive-sequence, AC phasor domain bus     | [Vr, Vi]
  infinite_bus | positive-sequence, AC phasor domain bus with fixed voltage | [Vr, Vi]
  emt_bus      | 3-phase bus with instantaneous voltages     | [Va, Vb, Vc]
  infinite_emt_bus | 3-phase bus with instantaneous voltages | [Va, Vb, Vc]
  control      | single control signal                       | [x]


### Devices Section

```json
"devices": [[..., ], [..., ], [..., ]]
```

The devices section is also an array of arrays. Each array represents
one device, and has exactly 5 items as follows.

  Index | Node item
  ------|----------------------------------------------------------------
  0     | Device class, a string identifier for what type of device this is. See table below for currently supported devices. A goal is to continually increase the number of supported devices.
  1     | Node numbers, an array of nodes connected with this device. The length of this array is fixed depending on the device class, but some nodes may be 0 if that functionality is not used. For example, GENROU class devices must have 3 nodes: bus, exciter signal, and governor signal. But the bus is the only one which must be connected to an actual node.
  2     | String ID disambiguator. Each device must have a unique combination of node numbers plus this string. It is recommended for most devices for this to only be 1-2 characters to facilitate converting to industry formats.
  3     | Array of initialization parameters. The length and meaning of these values is fixed and specified by the device class (see table below).
  4     | Extra, object with key-value pairs for any extra information to provide with the device. This can be an empty object {}. If you want to specify a non-default frequency or power base for this device, do it here with `"freq_base"` or `"va_base"`. If you add the item `"monitor":[]`, the variable values specified in the array will be recorded to an output channel.


### Table of Supported Device Classes

For now the following device classes are supported, with more to be
added in future versions. Note that the format version number and
revision number are extremely important, as the device list and the
exact lists of nodes and parameters may change with updates in different
versions.


  Class Name    | Description  | Nodes                  | Initialization Parameters
  --------------|--------------| -----------------------|-----------------------
  branch        | a basic algebraic pi model for a line or transformer | 2: Bus1, Bus2 | 4: R, X, G, B                         
  static_load   | a basic static ZIP load   | 1: Bus        | 6: Pz, Qz, Pi, Qi, Pp, Qp
  GENROU        | 6th order machine model   | 3: Bus, exciter_signal, governor_signal | 18: p0, q0, H, D, Ra, Tdop, Tdopp, Tqopp, Tqop, Xd, Xdp, Xdpp, Xq, Xqp, Xqpp, Xl, S10, S12
  bus_fault     | simple impedance-based fault at a bus | 2: Bus, control signal | 3: state0, R, X                                       


### Example File for 2-Bus System

```json
{"header": {
"format_version": 0, "format_revision": 1, "case_name": "two-bus test case 1", "case_data": "2/28/2025",
"case_description": "A two-bus test case for demonstrating the dynamics format.",
"comments":"The case is set up to monitor the voltage at both buses and the machine angle and speed",
"freq_base":60, "va_base": 100e6
},

"nodes": [

[1, "bus", "Bus 1", [0.994988, 0.099997], 115e3, {"monitor": [ "Vr", "Vi"]}],
[2, "infinite_bus", "Bus 2", [1.0, 0], 115e3, {}]

],

"devices": [

["branch", [1, 2], "1", [0, 0.1, 0, 0], {}],
["GENROU", [1, 0, 0], "1", [1, 0.05013, 3, 0, 0, 7, 0.04, 0.05, 0.75, 2.1, 0.2, 0.18, 0.5, 0.5, 0.18, 0.15, 0, 0], {"monitor": ["delta", "omega"]}],
["bus_fault": [1], "1", [0, 0, 1e-3], {}]

] }
```

### Appendix: Output File Format

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
