# Input file for GridKit EMT application

## Root elements

   Name                         | Value
 -------------------------------|-------------------------------------------------------
  `system_model_file`           | Path to the system model file[^1]
  `state_file`                  | Path to the state file[^1] (optional)
  `dt_monitor`                  | Monitor output time interval (default: 0)
  `tmax`                        | Final time
  `rel_tol`                     | Relative solver tolerance (default: 1.0e-2)
  `abs_tol`                     | Absolute solver tolerance override (default: 0, model tolerances)
  `dt_fixed`                    | Fixed solver time step size (default: 0, adaptive)
  `max_steps`                   | Maximum solver steps per output interval (default: unlimited)
  `max_order`                   | Maximum IDA integration method order (default: 2)
  `suppress_algebraic_errors`   | Leave the algebraic variables out of the error test (default: true)
  `events`                      | Array of [events](#events)
  `output_file`                 | Path to output (CSV) file (optional)
  `output_layout`               | One of { "wide", "long" } (default: "wide")

[^1]: See [EMT input format](../../GridKit/Model/EMT/INPUT_FORMAT.md)

## Events

An event sets values of the [state file](../../GridKit/Model/EMT/INPUT_FORMAT.md#state-file)
at a time. List events in chronological order. Events at the same time take
effect together, in file order, with one IDA restart.

   Name              | Value
 --------------------|-------------------------------------------------------
  `time`             | Event time
  `buses`            | Values by bus number (optional)
  `devices`          | Values by device id (optional)
