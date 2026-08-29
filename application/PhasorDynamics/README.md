# Input file for GridKit phasor dynamics application

## Root elements

   Name                | Value
 ----------------------|-------------------------------------------------------
  `system_model_file`  | Path to the system model file[^1]
  `dt_monitor`         | Monitor output time interval for recorded simulation results (default: 0, no intermediate monitoring)
  `tmax`               | A floating-point value for max time
  `ida`                | IDA solver options (optional; see [IDA options](#ida-options))
  `consistent_ic_type` | IDA consistent initial condition calculation type; one of { "y", "ya_ydp" } (default: "ya_ydp")
  `fault_bus`          | Bus id used to relocate the system model's first bus fault (optional)
  `events`             | An array of event groups (see [Events](#events) below)
  `output_file`        | Path to output (CSV) file (optional)
  `step_trace_file`    | Path to optional IDA internal-step trace CSV
  `reference_file`     | A string containing the name of the case (optional)
  `error_type`         | One of { "relative" (default), "absolute" }
  `error_tolerance`    | A floating-point value for highest allowable total error (default: 1.0e-4)
  `abs_err_threshold`  | A floating-point value for the smallest value at which to scale relative error (default: machine epsilon for double-precision)

[^1]: See system model [case format](../../GridKit/Model/PhasorDynamics/INPUT_FORMAT.md)

## IDA options

All IDA options are optional. Omitted options retain the application or IDA
default shown below.

 Name                        | Default
 ----------------------------|-----------------------
 `rel_tol`                   | `1.0e-7`
 `abs_tol`                   | `1.0e-9`; use `0` for model-specific tolerances
 `fixed_step`                | Adaptive stepping
 `init_step`                 | Estimated by IDA
 `min_step`                  | No minimum
 `max_step`                  | Unbounded
 `max_order`                 | `5`
 `max_num_steps`             | `500`
 `max_err_test_fails`        | `10`
 `suppress_alg`              | `false`
 `max_nonlin_iters`          | `4`
 `max_conv_fails`            | `10`
 `nonlin_conv_coef`          | `0.33`
 `max_num_steps_ic`          | `5`
 `max_num_jacs_ic`           | `4`
 `max_num_iters_ic`          | `10`
 `max_backs_ic`              | `100`
 `line_search_off_ic`        | `false`
 `nonlin_conv_coef_ic`       | `0.0033`
 `step_tolerance_ic`         | IDA default
 `linear_solution_scaling`   | `true`
 `delta_cj_lsetup`           | `0.25`
 `klu_ordering`              | `"amd"`; one of `"amd"`, `"colamd"`, or `"natural"`

`fixed_step` cannot be combined with `init_step`, `min_step`, or `max_step`.
It retains GridKit's existing fixed-step tolerance and nonlinear-solver
behavior.

```json
"ida": {
    "rel_tol": 1.0e-6,
    "max_num_steps": 2000,
    "klu_ordering": "amd"
}
```

For compatibility, the former top-level `rel_tol`, `abs_tol`, `dt_fixed`,
`max_steps`, and `klu_ordering` fields remain accepted. A study must use either
those fields or `ida`, not both.

## Internal-step trace

When `step_trace_file` is set, `DynamicSimulation` runs IDA in one-internal-step
mode and records each accepted step. Normal monitor and callback output remains
at the requested `dt_monitor` times. The CSV contains the event segment, step
time and size, next proposed step, current and proposed order, and cumulative
step, residual, Jacobian, error-test, and nonlinear-solver counters.

## Events

Each event group describes a system event that occurs at a given time point

   Name              | Value
 --------------------|-------------------------------------------------------
  `time`             | A floating point value for time event occurs
  `type`             | Event type (one of { "fault_on", "fault_off" })
  `element_id`       | An integer value referencing the element associated with the event (e.g., bus fault id)

When `fault_bus` is specified, the first bus fault (`element_id: 0`) is moved to
that bus before the system is constructed. Other event element ids retain their
existing meaning.
