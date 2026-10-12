# Input file for GridKit phasor dynamics application

## Root elements

   Name                | Value
 ----------------------|-------------------------------------------------------
  `system_model_file`  | Path to the system model file[^1]
  `dt_monitor`         | Monitor output time interval for recorded simulation results (default: 0, no intermediate monitoring)
  `tmax`               | A floating-point value for max time
  `rel_tol`            | Relative solver tolerance (default: 1.0e-7)
  `abs_tol`            | Absolute solver tolerance override (default: 1.0e-9)
  `signal_history`     | Optional delayed-signal accuracy and storage settings (see below)
  `dt_fixed`           | Fixed solver time step size, or 0 for adaptive stepping (default: 0)
  `max_steps`          | Maximum number of solver time steps, 0 for the IDA default, or a negative number for unlimited steps (default: 0)
  `max_order`          | Maximum IDA integration method order from 1 to 5 (default: 5; fixed stepping is capped at 2)
  `consistent_ic_type` | IDA consistent initial condition calculation type; one of { "y", "ya_ydp" } (default: "ya_ydp")
  `events`             | An array of event groups (see [Events](#events) below)
  `output_file`        | Path to output (CSV) file (optional)
  `reference_file`     | A string containing the name of the case (optional)
  `error_type`         | One of { "relative" (default), "absolute" }
  `error_tolerance`    | A floating-point value for highest allowable total error (default: 1.0e-4)
  `abs_err_threshold`  | A floating-point value for the smallest value at which to scale relative error (default: machine epsilon for double-precision)

[^1]: See system model [case format](../../GridKit/Model/PhasorDynamics/INPUT_FORMAT.md)

## Signal history

Delayed signals are declared in the case's `signals` array. Both applications
enable history control when positive delays are present. An optional solver-file
object overrides its defaults:

```json
"signal_history": {
  "rel_tol": 1e-7,
  "abs_tol": 1e-9,
  "compress": true
}
```

History tolerances default to the solver tolerances; a nonpositive solver absolute
tolerance uses `1e-9` for history. Set `compress` to `false` to retain every accepted
interpolation segment within the required delay window. Accuracy controls the
step limit, independently of the shortest delay. A fixed step that fails the
history error check reports an error. See the
[model equations and initialization](../../GridKit/Model/PhasorDynamics/SignalDelay/README.md).

## Events

Each event group describes a system event that occurs at a given time point

   Name              | Value
 --------------------|-------------------------------------------------------
  `time`             | A floating point value for time event occurs
  `type`             | Event type (one of { "fault_on", "fault_off" })
  `bus`              | Number of the faulted bus (omitted for `ContingencyAnalysis`, which faults every bus in turn using the same `R` and `X`)
  `R`                | Fault resistance in per unit (required for "fault_on")
  `X`                | Fault reactance in per unit (required for "fault_on")
