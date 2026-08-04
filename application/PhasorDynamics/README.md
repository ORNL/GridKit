# Input file for GridKit phasor dynamics application

## Root elements

   Name                | Value
 ----------------------|-------------------------------------------------------
  `system_model_file`  | Path to the system model file[^1]
  `dt_monitor`         | Monitor output time interval for recorded simulation results, 0 for no intermediate monitoring, or a negative number to output at solver-selected steps (default: 0)
  `tmax`               | A floating-point value for max time
  `rel_tol`            | Relative solver tolerance (default: 1.0e-7)
  `abs_tol`            | Absolute solver tolerance override (default: 1.0e-9)
  `dt_fixed`           | Fixed solver time step size, or 0 for adaptive stepping (default: 0)
  `max_steps`          | Maximum number of solver time steps, 0 for the IDA default, or a negative number for unlimited steps (default: 0)
  `consistent_ic_type` | IDA consistent initial condition calculation type; one of { "y", "ya_ydp" } (default: "ya_ydp")
  `ida_max_order`      | Optional integer from 1 to 5 limiting IDA's maximum integrator order
  `ida_max_dt`         | Optional positive floating point value limiting IDA's internal max step
  `events`             | An array of event groups (see [Events](#events) below)
  `output_file`        | Path to monitor output (CSV) file; omitted means monitors are disabled (optional)
  `ida_stats`          | Optional path to IDA aggregate statistics JSON
  `ida_steps`          | Optional path to IDA accepted-step JSON
  `reference_file`     | A string containing the name of the case (optional)
  `error_type`         | One of { "relative" (default), "absolute" }
  `error_tolerance`    | A floating-point value for highest allowable total error (default: 1.0e-4)
  `abs_err_threshold`  | A floating-point value for the smallest value at which to scale relative error (default: machine epsilon for double-precision)

[^1]: See system model [case format](../../GridKit/Model/PhasorDynamics/INPUT_FORMAT.md)

## Events

Each event group describes a system event that occurs at a given time point

   Name              | Value
 --------------------|-------------------------------------------------------
  `time`             | A floating point value for time event occurs
  `type`             | Event type (one of { "fault_on", "fault_off" })
  `element_id`       | An integer value referencing the element associated with the event (e.g., bus fault id)
