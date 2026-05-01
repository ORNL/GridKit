# PDSim (Phasor Dynamics Simulation)


PDSim runs a phasor-dynamics simulation defined by a *solver file*
(`*.solver.json`) and a *case file* (`*.case.json`). The case file specifies
the system model, and the solver file specifies the simulation enviornment, including
any runtime events scheduled during the simulation.

The runtime event system (action vocabulary, dispatch model, schedule
semantics) is documented in [`Model/PhasorDynamics/EVENTS.md`](../../GridKit/Model/PhasorDynamics/EVENTS.md).

## Usage

```
pdsim <solver-file.solver.json>
```

Relative paths in `case_file`, `output_file`, and `reference_file` are
resolved against the solver file's directory. Absolute paths are used as
given.

## Solver file format

### Top-level fields

Name              | Required | Description
------------------|----------|------------------------------------------------------
`format_version`  | no       | Format version integer. Default `1`. Parser rejects unknown versions.
`case_file`       | yes      | Path to the case file
`dt`              | yes      | Reporting time step in seconds
`tmax`            | yes      | End time in seconds
`schedule`        | no       | Ordered array of events (see below). May be empty or omitted.
`output_file`     | no       | Path to monitor output (CSV)
`reference_file`  | no       | Path to reference CSV for validation
`error_tolerance` | no       | Maximum allowed error vs reference. Default `1e-4`.

### Schedule

The `schedule` field is an array of events. Each event has the shape:

Name     | Required           | Description
---------|--------------------|------------------------------------------------------
`time`   | yes                | Simulation time at which the event fires, in seconds. Must be `≤ tmax`.
`target` | yes                | Routing key of the component receiving the event (bus `number` or device `id`)
`action` | yes                | Action name; one of `open`, `close`, `fault`, `clear`
`params` | iff action carries | Action payload (e.g. `{ R, X, percent }` for `fault`)

The action vocabulary, dispatch model, schedule canonicalization (stable
sort by time, last-listed wins on conflicts), and same-time semantics are
documented in
[`EVENTS.md`](../../GridKit/Model/PhasorDynamics/EVENTS.md).

## Output and validation

If `output_file` is given, it sets the destination for the case file's
default CSV monitor sink. Other monitor sinks declared in the case file
are unaffected. Output formatting (channels, delimiter) is controlled by
the case file's monitor declarations.

If both `output_file` and `reference_file` are given, PDSim compares the
output against the reference and exits with status `0` if the maximum
error is within `error_tolerance`, `1` otherwise.

## Example

```json
{
    "format_version":  1,
    "case_file":       "texas.json",
    "dt":              0.00416666666666,
    "tmax":            20.0,
    "schedule": [
        { "time": 10.0,  "target": "1001",           "action": "fault",
          "params": { "R": 0.0, "X": 1.0e-5 } },
        { "time": 10.15, "target": "1001",           "action": "clear" },
        { "time": 10.15, "target": "BR_1001_1064_1", "action": "open"  }
    ],
    "output_file":     "texas.csv",
    "reference_file":  "texas.ref.csv",
    "error_tolerance": 1.0e-4
}
```

This run faults bus `1001` at `t=10`, clears the fault and trips line
`BR_1001_1064_1` simultaneously at `t=10.15` (a 9-cycle fault followed by
isolation of the faulted equipment), and integrates to `t=20`. The two
events at `t=10.15` batch into one IDA re-init. Output is compared against
the reference file with tolerance `1e-4`.
