# Optimal Dispatch

Each example solves the optimal power flow of a case in `cases/PhasorDynamics/`
with the limits and costs of its MATPOWER case, writes the solution state, and
checks that PhasorDynamics starts from it in steady state.

 Case          | Buses | Generators | Ipopt iterations | Cost [$/h]
 --------------|-------|------------|------------------|-----------
 `IEEE39`      | 39    | 10         | 16               | 40086.9
 `Hawaii`      | 37    | 39         | 20               | 4151.17
 `ACTIVSg2000` | 2000  | 334        | 35               | 1.22901e6

Run an example from its build directory, then start a dynamic study from the
solution:

```shell
OptimalDispatch IEEE39.solver.json
```

```json
{
  "system_model_file": "IEEE39.case.json",
  "state_file": "IEEE39.state.json",
  "tmax": 10.0
}
```

## Pass and fail criteria

- `OptimalDispatch_<case>` passes when Ipopt converges and the state is
  written.
- `OptimalDispatch_<case>_steady_start` applies the state to the case and
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
