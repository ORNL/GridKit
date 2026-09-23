# Input file for GridKit phasor dynamics application

## Root elements

   Name                | Value
 ----------------------|-------------------------------------------------------
  `system_model_file`  | Path to the system model file[^1]
  `state_file`         | Path to a state file[^2] that sets the initial operating point (optional)
  `dt_monitor`         | Monitor output time interval for recorded simulation results (default: 0, no intermediate monitoring)
  `tmax`               | A floating-point value for max time (required by `DynamicSimulation` and `ContingencyAnalysis`)
  `rel_tol`            | Relative solver tolerance (default: 1.0e-7)
  `abs_tol`            | Absolute solver tolerance override (default: 1.0e-9)
  `dt_fixed`           | Fixed solver time step size, or 0 for adaptive stepping (default: 0)
  `max_steps`          | Maximum number of solver time steps, 0 for the IDA default, or a negative number for unlimited steps (default: 0)
  `max_order`          | Maximum IDA integration method order from 1 to 5 (default: 5; fixed stepping is capped at 2)
  `consistent_ic_type` | IDA consistent initial condition calculation type; one of { "y", "ya_ydp" } (default: "ya_ydp")
  `events`             | An array of event groups (see [Events](#events) below, default: none)
  `output_file`        | Path to output (CSV) file (optional)
  `output_state_file`  | Path to the state file[^2] the study writes (required by `OptimalDispatch`)
  `dispatch_file`      | Path to the MATPOWER case[^3] with limits and costs (required by `OptimalDispatch`)
  `ipopt`              | Ipopt options by name (optional, `OptimalDispatch` only)
  `reference_file`     | A string containing the name of the case (optional)
  `error_type`         | One of { "relative" (default), "absolute" }
  `error_tolerance`    | A floating-point value for highest allowable total error (default: 1.0e-4)
  `abs_err_threshold`  | A floating-point value for the smallest value at which to scale relative error (default: machine epsilon for double-precision)

[^1]: See system model [case format](../../GridKit/Model/PhasorDynamics/INPUT_FORMAT.md)
[^2]: See the [state format](../../GridKit/Model/STATE.md)
[^3]: See the optimal power flow [data](../../GridKit/Model/OptimalPowerFlow/README.md#data)

## Events

Each event group describes a system event that occurs at a given time point

   Name              | Value
 --------------------|-------------------------------------------------------
  `time`             | A floating point value for time event occurs
  `type`             | Event type (one of { "fault_on", "fault_off" })
  `element_id`       | An integer value referencing the element associated with the event (e.g., bus fault id)

## Optimal dispatch

`OptimalDispatch` solves the AC optimal power flow of the case with the limits
and costs of its MATPOWER case and writes the solution to `output_state_file`.
`DynamicSimulation` starts from that state through `state_file`.

```shell
OptimalDispatch IEEE39.solver.json
```

 Case device                                  | Optimal power flow component
 ---------------------------------------------|-----------------------------
 `Bus`, `BusInfinite`                         | `Bus`, infinite for `BusInfinite`
 `Branch`                                     | `Branch` with the same parameters
 `Genrou`, `Gensal`, `GenClassical`, `Regca`  | `Generator`
 `LoadZIP`                                    | `Load` with its demand from the state
 `LoadZ`                                      | `Shunt` with $G + jB = 1 / (R + jX)$

Controllers, exciters, governors, stabilizers, and faults do not enter the
optimal power flow. PhasorDynamics initialization checks the dispatch against
their limits.

Ipopt uses exact sparse derivatives. The application sets `bound_relax_factor`
to 0, because PhasorDynamics initialization rejects limits that relaxed bounds
violate. Options under `ipopt` are applied after it.
