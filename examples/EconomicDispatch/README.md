# Economic Dispatch

This example source is retained but is not built, installed, or registered
with CTest. PhasorDynamics initialization from an OPF state is deferred.

The source solves the optimal power flow of a case in `cases/PhasorDynamics/`
with the limits and costs of its MATPOWER case and writes the solution state.
Previously recorded OPF results:

 Case          | Buses | Generators | Ipopt iterations | Cost [$/h]
 --------------|-------|------------|------------------|-----------
 `IEEE39`      | 39    | 10         | 16               | 40086.9
 `Hawaii`      | 37    | 39         | 20               | 4151.17
 `ACTIVSg2000` | 2000  | 334        | 35               | 1.22901e6

## Input

Dispatch requires `system_model_file`, `dispatch_file`, and `output_state_file`;
`state_file` and `ipopt` are optional. See the
[IEEE39 input](IEEE39/IEEE39.solver.json).
Input paths are relative to the solver file; output is relative to the working
directory. Missing state values use the case operating point; supplied currents
must include both `ir` and `ii`.

Machines marked `online: false` in the input state are currently rejected.

 Case device                                  | Optimal power flow component
 ---------------------------------------------|-----------------------------
 `Bus`, `BusInfinite`                         | `Bus`, infinite for `BusInfinite`
 `Branch`                                     | `Branch` with the same parameters
 `Genrou`, `Gensal`, `GenClassical`, `Regca`  | `Generator`
 `LoadZIP`                                    | `Load` with its demand from the state
 `LoadZ`                                      | `Shunt` with $G + jB = 1 / (R + jX)$

Controllers, exciters, governors, stabilizers, and faults do not enter the
optimal power flow.

Ipopt uses exact sparse derivatives. The application sets `bound_relax_factor`
to 0 to enforce the original bounds. Options under `ipopt` are applied after it.

## Pass criterion

When enabled, `EconomicDispatch_<case>` passes when Ipopt converges and the
state is written.

## Data

- `IEEE39`: MATPOWER `case39`.
- `Hawaii`: the MATPOWER export of the Texas A&M Hawaii40 case.
- `ACTIVSg2000`: MATPOWER `case_ACTIVSg2000`. The static generators are
  `LoadZIP` devices with negative demand in the case, so they stay fixed.
