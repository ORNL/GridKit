# Input file for GridKit EMT application

## Root elements

   Name                | Value
 ----------------------|-------------------------------------------------------
  `system_model_file`  | Path to the system model file[^1]
  `state_file`         | Path to the state file[^1] (optional)
  `dt_monitor`         | Monitor output time interval (default: 0)
  `tmax`               | Final time
  `rel_tol`            | Relative solver tolerance (default: 1.0e-5)
  `abs_tol`            | Absolute solver tolerance override (default: 0, model tolerances)
  `dt_fixed`           | Fixed solver time step size (default: 0, adaptive)
  `max_steps`          | Maximum solver steps per output interval (default: unlimited)
  `max_order`          | Maximum IDA integration method order (default: 5)
  `events`             | Array of [events](#events)
  `output_file`        | Path to output (CSV) file (optional)
  `output_layout`      | One of { "wide", "long" } (default: "wide")

[^1]: See [EMT input format](../../GridKit/Model/EMT/INPUT_FORMAT.md)

## Events

   Name              | Value
 --------------------|-------------------------------------------------------
  `time`             | Event time
  `type`             | One of { "switch_open", "switch_close" }
  `device`           | `Switch` id
