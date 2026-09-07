# Input file for GridKit EMT application

`EMTDynamicSimulation` requires SUNDIALS KLU and a sparse model Jacobian
(enable Enzyme for the EMT components). It reports the selected solver and
Jacobian size and rejects configurations that would fall back to dense.

Before consistent initialization, the assembled model derives its differential
variables from `F_yp` and checks the initial-value Jacobian
`[F_yp(:, differential), F_y(:, algebraic)]`. A structural matching failure or
numerical singularity is reported with component paths and local indices.
This validates the current variable partition and operating point; it does
not perform index reduction. The check and IDA variable tags are refreshed
at each restart, including switch events.

## Root elements

   Name                | Value
 ----------------------|-------------------------------------------------------
  `system_model_file`  | Path to the system model file[^1]
  `state_file`         | Optional path to an EMT operating-point [state file](../../GridKit/Model/EMT/STATE.md), relative to the solver file; uses component output names and Boolean switch status
  `dt_monitor`         | Monitor output time interval for recorded simulation results (default: 0, no intermediate monitoring)
  `tmax`               | A floating-point value for max time
  `rel_tol`            | Relative solver tolerance (default: 1.0e-7)
  `abs_tol`            | Absolute solver tolerance override (default: 1.0e-9)
  `mu`                | Positive finite CommonMath smoothing scale (default: 240); configured before model construction
  `signal_values`     | Optional object overriding declared constant signals by qualified path, e.g. `{"dc_4": 28918.846170570516}`; cannot override component outputs
  `dt_fixed`           | Fixed solver time step size, or 0 for adaptive stepping (default: 0)
  `max_steps`          | Maximum number of solver time steps, 0 for the IDA default, or a negative number for unlimited steps (default: 0)
  `consistent_ic_type` | IDA consistent initial condition calculation type; one of { "y", "ya_ydp" } (default: "ya_ydp")
  `events`             | An array of event groups (see [Events](#events) below)
  `output_file`        | Path to output (CSV) file (optional)
  `state_output_file`  | Optional CSV of every DAE variable and derivative at the monitor times, with a companion `.csv.json` index map. Output paths are relative to the working directory. Event times include pre-event and post-event rows.
  `reference_file`     | A string containing the name of the case (optional)
  `error_type`         | One of { "relative" (default), "absolute" }
  `error_tolerance`    | A floating-point value for highest allowable total error (default: 1.0e-4)
  `abs_err_threshold`  | A floating-point value for the smallest value at which to scale relative error (default: machine epsilon for double-precision)

[^1]: See system model [case format](../../GridKit/Model/EMT/INPUT_FORMAT.md)

The `mu` option follows `lukel/mu-control-dev`: it sets the process-wide
`Math::MU<RealT>`, affecting all CommonMath primitives, not only PWM. Configure
it before constructing models or starting workers, and keep it fixed during
the run. Separate application processes can use different values.

## Events

Each event group describes a system event that occurs at a given time point.
A switch event changes the Jacobian sparsity pattern, so the driver
rediscovers the structure, rebuilds the linear solver, and reinitializes the
integrator at the event time.

   Name              | Value
 --------------------|-------------------------------------------------------
  `time`             | A floating point value for time event occurs
  `type`             | Event type (one of { "switch_open", "switch_close" })
  `element_id`       | String ID of the `Switch` component associated with the event
