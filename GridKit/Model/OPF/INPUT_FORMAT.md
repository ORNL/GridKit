# OPF case format specification

Version 0.1

## Overview

This document describes a data format for steady state network cases used in
power flow and optimal power flow studies. The format is implemented as UTF-8
encoded JSON and uses the file extension `.opf.json`.

The format stores no operating state. Bus voltages, device dispatch, load
demand, service statuses, and tap positions belong to the companion state
format with the extension `.state.json`. State records are matched to this
format through bus numbers and device id strings.

All electrical quantities are given in per unit on the system power base
`va_base`. Bus voltage bases are given by the bus parameter `kv`. Missing
limit parameters are treated as unbounded.

## Format

The root element is an object with four keys `header`, `params`, `buses`,
and `devices`.

### Header

Contained in the `header` key is an object with the following items.

  Name               | Value
  -------------------|------------------------------------------------------
  `format_version`   | Non-negative integer indicating the format version
  `format_revision`  | Non-negative integer indicating the format revision
  `case_name`        | A string containing the name of the case
  `case_date_time`   | Optional string in the ISO 8601 format
  `case_description` | Optional string describing what is modeled in the case
  `case_comments`    | Optional string with additional notes as needed

### Params

Contained in the `params` key is an object with the following items.

  Name        | Value
  ------------|---------------------------------------------------------
  `freq_base` | System frequency base in hertz
  `va_base`   | System power base in volt amperes

### Buses

Contained in the `buses` key is an array of objects, each of which represents
a bus and has the following fields.

  Name     | Description
  ---------|--------------------------------------------------------------
  `class`  | A string indicating the class of bus
  `id`     | A string naming the bus, unique among all buses
  `params` | An object mapping parameter names to numerical values

Bus classes.

  Bus class | Description                             | Parameters
  ----------|-----------------------------------------|------------------------------
  `Bus`     | A standard bus                          | `number`, `kv`, `vmin`, `vmax`
  `Slack`   | The slack bus giving the angle reference | `number`, `kv`, `vmin`, `vmax`

Exactly one bus should have the class `Slack`.

Bus parameters.

  Name     | Description
  ---------|--------------------------------------------------------------
  `number` | Unique non-negative integer referenced by device bus maps
  `kv`     | Bus voltage base in kilovolts
  `vmin`   | Minimum voltage magnitude
  `vmax`   | Maximum voltage magnitude

### Devices

Contained in the `devices` key is an array of objects, each of which
represents a device and has the following fields.

  Name     | Description
  ---------|--------------------------------------------------------------
  `class`  | A string indicating the class of device
  `buses`  | An object mapping the device bus names to bus numbers
  `id`     | A string identifying the device, unique among all devices
  `params` | An object mapping parameter names to numerical values

Device classes.

  Device class | Description                          | Buses          | Parameters
  -------------|--------------------------------------|----------------|--------------------------------------------------------
  `Branch`     | Pi model line or transformer branch  | `from`, `to`   | `R`, `X`, `G`, `B`, `smax`
  `Generator`  | Dispatchable injection               | `bus`          | `pmin`, `pmax`, `qmin`, `qmax`, `mva`, `c0`, `c1`, `c2`
  `Load`       | Demand                               | `bus`          | `pmin`, `pmax`, `qmin`, `qmax`
  `Shunt`      | Fixed shunt compensation             | `bus`          | `G`, `B`

Device parameters.

  Name           | Description
  ---------------|-------------------------------------------------------
  `R`            | Series resistance
  `X`            | Series reactance
  `G`            | Total shunt conductance
  `B`            | Total shunt susceptance
  `smax`         | Apparent power rating applied at each terminal
  `pmin`, `pmax` | Active power limits
  `qmin`, `qmax` | Reactive power limits
  `mva`          | Machine power base in megavolt amperes
  `c0`, `c1`, `c2` | Generation cost coefficients

Branch parameters follow the same pi model conventions as the PhasorDynamics
`Branch` class. Generator cost is `c0 + c1 p + c2 p^2` in currency per hour
with active power `p` in per unit.

## Example File for a 2-Bus System

```json
{
    "header": {
        "format_version": 0,
        "format_revision": 1,
        "case_name": "Basic 2-bus OPF case",
        "case_description": "A two-bus case for demonstrating the OPF format"
    },
    "params": {
        "freq_base": 60.0,
        "va_base": 100e6
    },
    "buses": [
        {
            "class": "Bus",
            "id": "Bus_1",
            "params": {"number": 0, "kv": 115.0, "vmin": 0.95, "vmax": 1.05}
        },
        {
            "class": "Slack",
            "id": "Bus_2",
            "params": {"number": 1, "kv": 115.0, "vmin": 0.95, "vmax": 1.05}
        }
    ],
    "devices": [
        {
            "class": "Branch",
            "buses": {"from": 0, "to": 1},
            "id": "BR1",
            "params": {"R": 0.0, "X": 0.1, "G": 0.0, "B": 0.0, "smax": 2.5}
        },
        {
            "class": "Generator",
            "buses": {"bus": 0},
            "id": "DV1",
            "params": {
                "pmin": 0.0,
                "pmax": 1.2,
                "qmin": -0.6,
                "qmax": 0.6,
                "mva": 120.0,
                "c0": 0.0,
                "c1": 14.0,
                "c2": 0.02
            }
        },
        {
            "class": "Generator",
            "buses": {"bus": 1},
            "id": "EQ1",
            "params": {
                "pmin": -99.0,
                "pmax": 99.0,
                "qmin": -99.0,
                "qmax": 99.0,
                "mva": 1000.0,
                "c1": 25.0
            }
        },
        {
            "class": "Load",
            "buses": {"bus": 0},
            "id": "LD1",
            "params": {}
        }
    ]
}
```
