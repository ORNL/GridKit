# PDSim — Phasor Dynamics Study Runner

## Overview

PDSim is a command-line application for running parameterized phasor dynamics
simulations. It reads a JSON study file that specifies the system model,
simulation parameters, disturbance events, and output configuration.

Event definitions (e.g. bus faults) are specified in
the study file rather than in the system model JSON, keeping the model data
purely descriptive.

### Usage

```
PDSim <study-file.study.json>
```

## Study file format

The study file configures a simulation run. All paths are resolved relative
to the study file's directory.

### Top-level fields

  Name               | Required | Description
  -------------------|----------|------------------------------------------------------
  `format_version`   | No       | Format version integer (default: `1`)
  `study_name`       | No       | Short name for the study
  `study_description`| No       | Free-text description of the study
  `case_file`        | Yes      | Path to the system model JSON file
  `dt`               | Yes      | Simulation time step in seconds
  `tmax`             | Yes      | Simulation end time in seconds
  `events`           | No       | Array of event definitions (see below). Required when simulation includes disturbance events
  `schedule`         | No       | Ordered array of timed cues (see below). Required when simulation includes disturbance events
  `output`           | No       | Monitor output configuration object (see [Output](#output) below)
  `reference_file`   | No       | Path to reference CSV for validation
  `error_tol`        | No       | Maximum allowed error vs reference (default: `1e-4`)

### Output

The study file controls where and how monitored variables are written.
Which variables are monitored is determined by the `mon` fields on buses
and devices in the model JSON.

The `output` object has the following fields:

  Name               | Description
  -------------------|------------------------------------------------------
  `file_name`        | Output file path. If omitted, output is written to `stdout`
  `format`           | One of { `"CSV"`, `"JSON"`, `"YAML"` } (case-insensitive, default: `"CSV"`)
  `delim`            | Delimiter string for CSV output (default: `","`)

If the `output` key is omitted entirely, the default is CSV output to
`stdout`.

The CSV format has a header row followed by one row per time step.
The first column is always `Time [s]`, the second is `Solver Status`
(0 = success, positive = warning, negative = error), followed by one
column per monitored variable in the order they appear in the model file.

> In a future version, `output` should become an array of sink objects to
> support multiple simultaneous outputs (e.g., CSV to file and YAML to
> stdout).

### Events

Each event defines a disturbance that can be activated during the simulation.
Events are referenced by `id` in the schedule.

Each event has the following common fields:

  Name     | Description
  ---------|------------------------------------------------------
  `id`     | Unique string identifier for the event
  `type`   | Event type string (see supported types below)
  `params` | Object containing type-specific parameters

#### `bus_fault` — fault-to-ground at a bus

  Name  | Description
  ------|------------------------------------------------------
  `bus` | Bus number to fault
  `R`   | Fault resistance in per-unit (default: `0.0`)
  `X`   | Fault reactance in per-unit (default: `1e-5`)

#### `trip_line` — trip a transmission line or transformer *(not yet implemented)*

  Name   | Description
  -------|------------------------------------------------------
  `bus1` | Bus number at the first end of the branch
  `bus2` | Bus number at the second end of the branch
  `id`   | Device disambiguation string (required when multiple branches connect the same two buses)

### Schedule

The schedule is an ordered list of timed cues that trigger event actions:

  Name     | Description
  ---------|------------------------------------------------------
  `time`   | Simulation time to execute the action in seconds
  `event`  | References an event `id`
  `action` | `"on"` to apply, `"off"` to clear

Each entry in the schedule array is a **cue** — a timed trigger for an event
action. The simulator runs to each cue's time, executes the action,
reinitializes, and continues. After all cues, it runs to `tmax`.

### Validation

When both `output.file_name` and `reference_file` are provided, PDSim compares the
simulation output against the reference and reports whether the maximum error
is within `error_tol`. The process exit code is `0` on pass, `1` on
failure.

(NOTE: I think in future versions we should generalize the reference and error
specification with 'post-process' or something similar.)

## Example

```json
{
    "format_version": 1,
    "case_file": "ThreeBusBasic.json",
    "dt": 0.00416666666666,
    "tmax":  10,
    "events": [
        {"id": "fault_1", "type": "bus_fault", "params": {"bus": 2, "R": 0.0, "X": 1e-5}},
        {"id": "trip_1",  "type": "trip_line", "params": {"bus1": 1, "bus2": 2, "id": "BR1"}}
    ],
    "schedule": [
        {"time": 1.0,  "event": "fault_1", "action": "on"},
        {"time": 1.1,  "event": "fault_1", "action": "off"},
        {"time": 1.1,  "event": "trip_1",  "action": "on"}
    ],
    "output": {"file_name": "ThreeBus_fault_results.csv"},
    "reference_file": "ThreeBusBasic.ref.csv",
    "error_tol":  1e-4
}
```
