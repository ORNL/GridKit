# Input file for GridKit phasor dynamics application

## Root elements

   Name                | Value
 ----------------------|-------------------------------------------------------
  `partition`          | Optional partitioned integration settings (see below)
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

[^1]: See system model [case format](../../GridKit/Model/PhasorDynamics/INPUT_FORMAT.md)

## Events

Each event group describes a system event that occurs at a given time point

   Name              | Value
 --------------------|-------------------------------------------------------
  `time`             | A floating point value for time event occurs
  `type`             | Event type (one of { "fault_on", "fault_off" })
  `element_id`       | An integer value referencing the element associated with the event (e.g., bus fault id)

## Partitioned integration

`partition` contains `file` (a bus assignment `.partition.json`), `dt` (the
initial macro step), and optional `tol` (endpoint prediction mismatch relative
to `1 + |input|`; default 0 uses fixed steps).

Regions are advanced by multicolor Gauss-Seidel: regions that share no tie
branch form a color, and coupled regions advance in color order. A region starts
as soon as the coupled regions before it finish. Its boundary inputs ramp toward
regions already advanced in the step and extrapolate the others. Complete macro
steps are accepted or rejected together, and input slopes reset at events. A
rejected step restarts each region from its accepted state. Consistent coupling
at the start and after events uses Anderson-accelerated Gauss-Seidel sweeps
(KINSOL).

Output at `dt_monitor` comes from each region's interpolant, so it does not
limit the macro step; `tol` alone sets it.

`threads` (default 1) sets how many regions advance at once; values above 1
need a build configured with `-DGridKit_ENABLE_OPENMP=ON`. Results do not depend
on `threads`. Regional solves are memory-bound, so more threads than physical
performance cores slow a run down. The reported `Complete in` duration is elapsed wall time,
excluding setup.
