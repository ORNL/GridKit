# Input file for GridKit phasor dynamics application

## Root elements

   Name                | Value
 ----------------------|-------------------------------------------------------
  `system_model_file`  | Path to the system model file[^1]
  `dt_monitor`         | Monitor output time interval for recorded simulation results (default: 0, no intermediate monitoring)
  `tmax`               | A floating-point value for max time
  `rel_tol`            | Relative solver tolerance (default: 1.0e-7)
  `abs_tol`            | Absolute solver tolerance override (default: 1.0e-9)
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
  Name               | Value
 ---------------------|-------------------------------------------------------
  `system_model_file` | Path to the system model file[^1]
  `state_file`        | Path to the initial operating State file
  `dt_monitor`        | Monitor output time interval for recorded simulation results (default: 0, no intermediate monitoring)
  `tmax`              | A floating-point value for max time
  `rel_tol`           | Relative solver tolerance (default: 1.0e-7)
  `abs_tol`           | Absolute solver tolerance override (default: 1.0e-9)
  `dt_fixed`          | Fixed solver time step size, or 0 for adaptive stepping (default: 0)
  `events`            | An array of named input events (see [Events](#events) below)
  `output_file`       | Path to output (CSV) file (optional)
  `reference_file`    | A string containing the name of the case (optional)
  `error_type`        | One of { "relative" (default), "absolute" }
  `error_tolerance`   | A floating-point value for highest allowable total error (default: 1.0e-4)
  `abs_err_threshold` | A floating-point value for the smallest value at which to scale relative error (default: machine epsilon for double-precision)

[^1]: See system model [case format](../../GridKit/Model/PhasorDynamics/INPUT_FORMAT.md)

## Events

Each event assigns a scalar value to a named model input at a simulation time.
Events must be listed in chronological order. Adjacent events at the same time
are applied together before the solver is reinitialized.

   Name              | Value
 --------------------|-------------------------------------------------------
  `time`             | A floating point value for time event occurs
  `device_id`        | Case device ID receiving the input change
  `input`            | Name of the model input to change
  `value`            | Scalar value assigned to the input

For example, a six-cycle bus fault is represented by changing its `active`
input:

```json
"events": [
  { "time": 1.0, "device_id": "fault_1", "input": "active", "value": 1.0 },
  { "time": 1.1, "device_id": "fault_1", "input": "active", "value": 0.0 }
]
```
