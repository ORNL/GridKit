# PDSim — Phasor Dynamics Study Runner

## Overview

PDSim is a command-line application for running parameterized phasor dynamics
simulations. It reads a JSON solver file that specifies the case (system
model) to load, the simulation parameters, an optional schedule of runtime
cues, and output configuration.

Components that respond to runtime control (e.g. `BusFault`, `Branch`) are declared in
the case file as ordinary devices. The solver file's `schedule` then drives
their state changes by referencing each component's `id` from the case file.

### Usage

```
PDSim <study-file.solver.json>
```

## Solver file format

The solver file configures a simulation run. All paths are resolved relative
to the solver file's directory.

### Top-level fields

  Name               | Required | Description
  -------------------|----------|------------------------------------------------------
  `format_version`   | No       | Format version integer (default: `1`)
  `case_file`        | Yes      | Path to the case (system model) JSON file
  `dt`               | Yes      | Simulation time step in seconds
  `tmax`             | Yes      | Simulation end time in seconds
  `schedule`         | No       | Ordered array of timed cues (see below)
  `output_file`      | No       | Path to monitor output file (CSV). If omitted, writes to a default sink declared in the case file
  `reference_file`   | No       | Path to reference CSV for validation
  `error_tolerance`  | No       | Maximum allowed error vs reference (default: `1e-4`)

### Output

Which variables are monitored — and in what format — is configured by the
case file (the `mon` fields on buses and devices, and the monitor sink
declarations). The solver file's `output_file` only chooses the destination
path for the study's CSV sink; everything else about the output (format,
delimiter, monitored variables) is a model property.

### Schedule

The schedule is an ordered list of timed **cues**. Each cue addresses a
device declared in the case file by its `id` and asks it to perform an
action at the given simulation time.

  Name     | Description
  ---------|------------------------------------------------------
  `time`   | Simulation time to fire the cue, in seconds (must be ≤ `tmax`)
  `target` | Device `id` from the case file (e.g. a `BusFault` device's `id`)
  `action` | Action verb interpreted by the target device's class

The schedule must be non-decreasing in time; ties are allowed for
simultaneous cues (e.g. trip line A and B at the same instant). PDSim
runs to each distinct cue time, applies every cue scheduled at that
instant, re-initializes for a consistent IC, and continues. After the
last cue, it runs to `tmax`.

#### Devices that accept cues

  Class       | Action vocabulary | Effect
  ------------|-------------------|----------------------------------------------
  `BusFault`  | `on`, `off`       | Engage or clear the fault impedance at the bus
  `Branch`    | `on`, `off`       | Close or open the branch (in service / tripped)

A device is cue-targetable only if its `id` is referenced by the schedule.
Devices whose class does not accept cues (or unknown actions) cause
PDSim to fail with a clear error.

### Validation

When both `output_file` and `reference_file` are provided, PDSim compares
the simulation output against the reference and reports whether the
maximum error is within `error_tolerance`. The process exit code is `0`
on pass, `1` on failure.

## Example

Solver file:

```json
{
    "format_version":  1,
    "case_file":       "ThreeBusBasic.case.json",
    "dt":              0.00416666666666,
    "tmax":            10,
    "schedule": [
        { "time": 1.0, "target": "fault_1", "action": "on"  },
        { "time": 1.1, "target": "fault_1", "action": "off" }
    ],
    "output_file":     "ThreeBus_six_cycle_fault_bus1.csv",
    "reference_file":  "ThreeBusBasic.ref.csv",
    "error_tolerance": 1e-4
}
```

The corresponding case file declares the `BusFault` device whose `id` is
referenced as `target`:

```json
{
    "devices": [
        {
            "class":  "BusFault",
            "id":     "fault_1",
            "ports":  { "bus": 2 },
            "params": { "R": 0.0, "X": 1e-5, "state0": false }
        }
    ]
}
```
