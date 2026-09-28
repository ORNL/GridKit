# Economic Dispatch

Each example solves the optimal power flow of a case in `cases/PhasorDynamics/`
with the limits and costs of its MATPOWER case, writes the solution state, and
checks that PhasorDynamics starts from it in steady state.

 Case          | Buses | Generators | Ipopt iterations | Cost [$/h]
 --------------|-------|------------|------------------|-----------
 `IEEE39`      | 39    | 10         | 16               | 40086.9
 `Hawaii`      | 37    | 39         | 20               | 4151.17
 `ACTIVSg2000` | 2000  | 334        | 35               | 1.22901e6

Build and run the IEEE39 example from the repository root:

```shell
cmake --build build --target EconomicDispatch -j 10
cd build/examples/EconomicDispatch/IEEE39
../EconomicDispatch IEEE39.solver.json
```

Use the written state in a dynamics study:

```json
{
  "system_model_file": "IEEE39.case.json",
  "state_file": "IEEE39.state.json",
  "tmax": 10.0,
  "events": []
}
```

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
optimal power flow. PhasorDynamics initialization checks the dispatch against
their limits.

Ipopt uses exact sparse derivatives. The application sets `bound_relax_factor`
to 0, because PhasorDynamics initialization rejects limits that relaxed bounds
violate. Options under `ipopt` are applied after it.

## Pass and fail criteria

- `EconomicDispatch_<case>` passes when Ipopt converges and the state is
  written.
- `EconomicDispatch_<case>_steady_start` applies the state to the case and
  passes when both of these hold:
  - the largest PhasorDynamics residual at $t = 0$ is below its tolerance;
  - the largest change of any variable over 1 s, relative to one plus its
    initial magnitude, is below its tolerance, with the `DynamicSimulation`
    default solver tolerances.

 Case          | Residual | Tolerance | Drift   | Tolerance
 --------------|----------|-----------|---------|----------
 `IEEE39`      | 1.0e-13  | 2e-13     | 9.3e-15 | 2e-14
 `Hawaii`      | 1.3e-13  | 2.5e-13   | 2.4e-14 | 5e-14
 `ACTIVSg2000` | 5.1e-11  | 1e-10     | 4.1e-12 | 8e-12

## Data

- `IEEE39`: MATPOWER `case39`.
- `Hawaii`: the MATPOWER export of the Texas A&M Hawaii40 case.
- `ACTIVSg2000`: MATPOWER `case_ACTIVSg2000`. The static generators are
  `LoadZIP` devices with negative demand in the case, so they stay fixed.
