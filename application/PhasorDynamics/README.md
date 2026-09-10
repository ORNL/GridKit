# Input file for GridKit phasor dynamics application

## Root elements

- `system_model_file`: Path to the system model file.[^1]
- `dt_monitor`: Monitor output time interval for recorded simulation results.
  Default: 0, with no intermediate monitoring.
- `tmax`: A floating-point value for max time.
- `rel_tol`: Relative solver tolerance. Default: 1.0e-7.
- `abs_tol`: Absolute solver tolerance override. Default: 1.0e-9.
- `cosim_rel_tol`: Relative co-simulation interface residual tolerance.
  Default: 1.0e-8.
- `cosim_step_tol`: KINSOL scaled step tolerance for the co-simulation
  interface. Default: 1.0e-12.
- `cosim_function_error`: Estimated relative error of co-simulation residual
  evaluations, used to choose finite-difference perturbations. Default: 1.0e-10.
- `dt_fixed`: Fixed solver time step size, or 0 for adaptive stepping.
  Default: 0.
- `max_steps`: Maximum number of solver time steps, 0 for the IDA default,
  or a negative number for unlimited steps. Default: 0.
- `max_order`: Maximum IDA integration method order from 1 to 5. Default: 5;
  fixed stepping is capped at 2.
- `consistent_ic_type`: IDA consistent initial condition calculation type;
  one of `"y"`, `"ya_ydp"`. Default: `"ya_ydp"`.
- `events`: An array of event groups (see [Events](#events) below).
- `output_file`: Path to output file. Optional.
- `reference_file`: A string containing the name of the case. Optional.
- `error_type`: One of `"relative"` (default), `"absolute"`.
- `error_tolerance`: A floating-point value for greatest allowable total error.
  Default: 1.0e-4.
- `abs_err_threshold`: A floating-point value for the smallest value at which
  to scale relative error. Default: machine epsilon for double-precision.

[^1]: See system model [case format].

[case format](../../GridKit/Model/PhasorDynamics/INPUT_FORMAT.md)

## Events

Each event group describes a system event that occurs at a given time point.

- `time`: A floating-point value for the time the event occurs.
- `type`: Event type; one of `"fault_on"`, `"fault_off"`, `"in_service"`,
  `"out_of_service"`.
- `bus`: Number of the faulted bus. Omitted for `ContingencyAnalysis`, which
  faults every bus in turn using the same `R` and `X`.
- `R`: Fault resistance in per unit. Required for `"fault_on"`.
- `X`: Fault reactance in per unit. Required for `"fault_on"`.
- `device`: `id` of the branch to put in service (`"in_service"`) or take out of
  service (`"out_of_service"`).
