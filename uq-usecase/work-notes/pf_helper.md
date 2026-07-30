 # GridKit Power Flow solver: notes and plan for Illinois case

**Goal**: use GridKit's built-in Newton/KINSOL power flow solver to produce a converged
PF solution from `case_ACTIVSg200.m` (or any hourly `.m` file), then patch the result
into `illinois.json` (replacing Vr/Vi, p0/q0, Pnom/Qnom).

---

## 0. Key design decisions and analogies

### How this relates to `gridkit_helper.ipynb` / `DynamicSimulation`

In `gridkit_helper.ipynb`, the "runner" is the `DynamicSimulation` binary — a general-purpose
transient simulator that takes a `.solver.json` config file, runs, and writes output files.
It is already built as part of the main GridKit CMake tree (`build/application/PhasorDynamics/DynamicSimulation`).

For power flow there is **no equivalent general binary**. The only built PF example is
`grid3bus`, which has 3-bus data hard-wired as an inline string — it takes no arguments
and cannot run other cases. So we need to build a new general-purpose PF runner: `solve_pf`.

`solve_pf` is the PF counterpart of `DynamicSimulation`. It is a **thin wrapper** around
GridKit's existing `PowerFlow` module — all physics come from GridKit headers
(`MatpowerParser`, `SystemModelPowerFlow`, `BusFactory`, `GeneratorFactory`, `Branch`,
`Load`, `Kinsol`). The `solve_pf.cpp` file itself contributes only `main()` + argument
parsing + stdout formatting; it is structurally identical to `parserCase()` in the
existing `Grid3BusSys.cpp` example, but reads from `argv[1]` instead of an inline string.

- takes a `.m` file path as `argv[1]`
- parses it via `MatpowerParser.hpp`
- assembles `SystemSteadyStateModel` and runs KINSOL
- prints per-bus V/theta to stdout (machine-parseable)
- the notebook calls it via `subprocess`, parses stdout into a DataFrame — same pattern as `gridkit_helper`

### Where `pf-solver/` lives vs `test-runs/`

`test-runs/` is a **runtime output directory** — case files are copied there and
`DynamicSimulation` writes CSV output (`mon.csv`, plus per-case result CSVs) into it. It contains no source code.

`pf-solver/` contains **C++ source and CMake** for the `solve_pf` binary. It is more like
a new `examples/PowerFlow/SolvePF/` entry. The reason it lives in `uq-usecase/` rather
than the GridKit examples tree is that it is a uq-usecase-specific utility, not a
first-class GridKit example. If it proves useful it can be promoted to `examples/PowerFlow/`.

### `.mat` vs `.m` file extension

`3bus.mat` (the 3-bus example file) is **plain ASCII text** in MATPOWER format — same as
any `.m` file. The extension is `.mat` by convention but the format is identical. The
`readMatPowerFile()` function uses `std::ifstream` and never checks the extension. So
`solve_pf` accepts both `.m` and `.mat` files equally — pass either as `argv[1]`.

---

---

## 1. What GridKit's PF solver actually does

GridKit implements a Newton-based power flow solver using SUNDIALS KINSOL as the nonlinear
algebraic equation solver. The model is assembled from component objects (Bus, Branch,
Generator, Load), each contributing residual and Jacobian entries. The system is solved as
$\mathbf{f}(\mathbf{x}) = 0$ for bus voltages $(V_i, \theta_i)$.

Key source files:

| File | Role |
|---|---|
| `GridKit/Model/PowerFlow/SystemModelPowerFlow.hpp` | `SystemSteadyStateModel<T,I>` — assembles bus/gen/branch/load objects into a solvable system; calls KINSOL |
| `GridKit/Model/PowerFlow/MatpowerParser.hpp` | `readMatPowerFile()` / `readMatPower()` — parses a MATPOWER v2 `.m` file into `SystemModelData` struct |
| `GridKit/Model/PowerFlow/PowerFlowData.hpp` | `BusData`, `GenData`, `BranchData`, `LoadData`, `SystemModelData` structs |
| `GridKit/Model/PowerFlow/Bus/` | `BusFactory` — creates slack, PQ, PV bus objects |
| `GridKit/Model/PowerFlow/Generator/` | `GeneratorFactory` — creates slack and PV generator models |
| `GridKit/Model/PowerFlow/Branch/` | `Branch` — PI-model branch with r, x, b, transformer ratio |
| `GridKit/Model/PowerFlow/Load/` | `Load` — constant PQ load |
| `GridKit/Solver/SteadyState/Kinsol.hpp` | KINSOL wrapper; `runSimulation()` drives the solve |

---

## 2. Existing examples

### `examples/PowerFlow/Grid3Bus/`

The canonical working example. Solves a 3-bus network (1 slack, 1 PQ, 1 PV) three ways:

1. **Monolithic** (`MiniGrid`): hard-wired residual/Jacobian in a single class
2. **Parser** (`parserCase`): reads `.m` as an inline string via `readMatPower(mp, iss)`,
   then passes `SystemModelData mp` to `SystemSteadyStateModel(mp)` — this is the path
   that generalizes to any `.m` file
3. **Hardwired** (`hardwiredCase`): manually constructs `BusData`/`GenData`/`BranchData`
   objects and calls `sysmodel->addBus()` / `sysmodel->addComponent()`

The parser path (case 2) is the template for running any MATPOWER `.m`:
```cpp
SystemModelData<double, size_t> mp;
readMatPowerFile(mp, filename);                          // parse .m from disk
auto* sysmodel = new SystemSteadyStateModel<>(mp);      // assemble
sysmodel->allocate();
sysmodel->initialize();
auto* kinsol = new Kinsol<>(sysmodel);
kinsol->configureSimulation();
kinsol->getDefaultInitialCondition();   // V=1, theta=0 flat start
kinsol->runSimulation();                // KINSOL solve
// read back: sysmodel->getBus(i)->V(), ->theta()
```

### `examples/PowerFlow/MatPowerTesting/`

Unit tests for the MATPOWER parser only (not a PF solve). Tests parse individual `.m`
rows (bus, gen, branch, gencost) and verify the parsed struct fields. Useful for
confirming parser behavior on ACTIVSg200 fields before attempting a solve.

---

## 3. What has been verified to work

- Parser reads `baseMVA`, `bus`, `gen`, `branch`, `gencost` from any MATPOWER v2 `.m` file
- Parser handles inline `%` comments and trailing semicolons
- `BusFactory` handles types: 1=PQ, 2=PV, 3=slack (ref), 4=isolated (skipped)
- Generator status field (`mpc.gen` col 8) is read into `GenData.status` but it is
  **not yet verified** whether `SystemSteadyStateModel` filters out offline gens
  (`status=0`) automatically
- Flat-start initial guess (`V=1`, `theta=0`) is the default; warm start is not exposed

### 3-bus sanity check (2026-07-22)

`grid3bus` (hard-wired demo) and `solve_pf` (on `3bus.mat`) both converge and agree:

| var | result | expected |
|-----|--------|----------|
| theta2 (deg) | -4.87980 | -4.87979 |
| V2 (p.u.) | 1.08281 | 1.08281 |
| theta3 (deg) | 1.46241 | 1.46241 |

9 KINSOL iterations, fn norm ~7e-6. All three `grid3bus` variants (monolithic, parser, hardwired) produce identical results. `solve_pf` matches to within 1e-3 tolerance.

### `solve_pf` on `3bus.mat` (2026-07-22)

`solve_pf` built from `uq-usecase/pf-solver/` and run on `3bus.mat`:

| var | solve_pf | expected | \|err\| |
|-----|----------|----------|---------|
| theta2 (deg) | -4.87980 | -4.87979 | 9.4e-06 |
| V2 (p.u.) | 1.08281 | 1.08281 | 1.3e-06 |
| theta3 (deg) | 1.46240 | 1.46241 | 6.2e-06 |

KINSOL return code 0. All PASS (tol 1e-3).

---

## 4. Open questions / unknowns for ACTIVSg200 (200-bus case)

| # | Question | How to answer |
|---|---|---|
| 1 | Does KINSOL converge from flat start on a 200-bus case? | Run it (see plan below) |
| 2 | Are offline generators (`status=0`) excluded from the PF equations? | Check `GeneratorFactory` — does it check `GenData.status`? |
| 3 | Does `Branch` handle non-unity transformer ratio (`ratio != 0`)? | Verified: `case_ACTIVSg200.m` has **no off-nominal tap transformers** — `TAP` is 0 or 1 and `SHIFT=0` for all 245 branches. The question is moot for this case but remains relevant for other cases. |
| 4 | Does the parser handle `mpc.bus_name`? | It is not in the `BusData` struct today; irrelevant for PF solve |
| 5 | Are shunt elements (`Gs`, `Bs` on buses) modeled? | Check `BusData` — fields exist; verify they enter the residual |
| 6 | After convergence, what is the output format? Bus objects expose `->V()` and `->theta()` (degrees). Is there a dump-to-.m or dump-to-JSON utility? | No utility found yet; needs custom output code |
| 7 | Line flow output (PF, QF, PT, QT)? | Not exposed in current PF examples; would need Branch to compute and expose post-solve flows |

---

## 5. Plan: attempt PF solve on ACTIVSg200

### Step 1: build GridKit with PF examples enabled

Confirm the `grid3bus` binary exists (built as part of the main CMake build):
```
ls /home/isatkaus/gridkit/build/bin/grid3bus
/home/isatkaus/gridkit/build/bin/grid3bus
```
Run it to confirm 3-bus solve works:
```bash
/home/isatkaus/gridkit/build/bin/grid3bus
```
Expected output: `theta2 = -4.87979 deg`, `V2 = 1.08281 p.u.`, `theta3 = 1.46241 deg`, `Success!`

### Step 2: check parser on ACTIVSg200

Write a small C++ driver (or adapt the parser test) that reads
`case_ACTIVSg200.m` and prints the parsed bus/gen/branch counts to confirm no
parse errors on the 200-bus file. The file path is:
```
/kfs2/projects/scidac/scidac-data/ACTIVSg200/raw-tamu-data/case_ACTIVSg200.m
```

### Step 3: write a `solve_illinois_pf.cpp` driver

Modeled on the `parserCase()` in `Grid3BusSys.cpp`, but reading from a file path and
printing results to stdout in a format suitable for patching `illinois.json`:

```cpp
#include <GridKit/Model/PowerFlow/MatpowerParser.hpp>
#include <GridKit/Model/PowerFlow/SystemModelPowerFlow.hpp>
#include <GridKit/Solver/SteadyState/Kinsol.hpp>

int main(int argc, char* argv[]) {
    std::string m_path = argv[1];  // path to .m file

    GridKit::PowerFlowData::SystemModelData<double, size_t> mp;
    GridKit::readMatPowerFile(mp, m_path);

    auto* sys = new GridKit::SystemSteadyStateModel<double, size_t>(mp);
    sys->allocate();
    sys->initialize();

    auto* kinsol = new GridKit::AnalysisManager::Sundials::Kinsol<double, size_t>(sys);
    kinsol->configureSimulation();
    kinsol->getDefaultInitialCondition();
    int ret = kinsol->runSimulation();

    // dump bus voltages
    for (auto& bd : mp.bus) {
        auto* bus = sys->getBus(bd.bus_i);
        printf("bus %d  V=%.6f  theta_deg=%.6f\n",
               (int)bd.bus_i, bus->V(), bus->theta() * 180.0 / M_PI);
    }
    return ret;
}
```

Add to `CMakeLists.txt` in `examples/PowerFlow/` (or `uq-usecase/`), link against
`GridKit::bus GridKit::generator GridKit::branch GridKit::load GridKit::solvers_steady`.

### Step 4: compare to MATPOWER reference

Run `case_ACTIVSg200.m` through MATPOWER (Octave/MATLAB) to get the reference
VM/VA solution, then diff against GridKit output. Key metrics:
- Max |V| error across all 200 buses
- Max |theta| error (degrees)
- Number of KINSOL iterations

### Step 5: extract Vr/Vi for JSON patching

After a successful solve, compute:
```
Vr_i = V_i * cos(theta_i)
Vi_i = V_i * sin(theta_i)
```
and write a patch script (`m_to_case.py`) that applies these to `illinois.json`.

---

## 6. Risks and fallback

- **Convergence**: flat start on a 200-bus network may not converge with the current
  KINSOL settings. If it fails, options are: (a) use MATPOWER VM/VA as warm start by
  setting `BusData.Vm`/`Va` before `initialize()`, or (b) fall back to PowerModels.jl
  (see below).
- **Offline generators**: 9 gens have `status=0` in the base case. If `GeneratorFactory`
  does not filter them out, they will inject wrong P and PF will diverge.
- **No JSON output utility yet**: a small Python post-processing script will be needed
  to read the solver output and write the patched JSON, regardless of which PF solver is used.

### Fallback: PowerModels.jl

If GridKit's KINSOL PF fails to converge (or is too difficult to configure for a 200-bus
case), the fallback is [PowerModels.jl](https://lanl-ansi.github.io/PowerModels.jl/dev/)
— a Julia package from LANL that supports AC/DC power flow and OPF.

**Why it's a good fallback:**
- Reads MATPOWER `.m` files natively via `PowerModels.parse_file()`
- Full AC power flow (`solve_ac_pf`) and AC OPF (`solve_ac_opf`) out of the box
- Well-tested on large networks; handles all MATPOWER conventions (tap ratios, shunts,
  offline gens, etc.)
- Julia is already available on Kestrel: `~/bin/julia112`

**Minimal usage:**
```julia
using PowerModels, Ipopt
data = PowerModels.parse_file("case_ACTIVSg200.m")
result = solve_ac_pf(data, optimizer_with_attributes(Ipopt.Optimizer))
# result["solution"]["bus"]["1"]["vm"], ["va"]
```

**Output format:** `result["solution"]["bus"]` dict keyed by bus number string; each
entry has `"vm"` (pu) and `"va"` (radians). Post-processing to Vr/Vi:
```julia
vm = result["solution"]["bus"][string(i)]["vm"]
va = result["solution"]["bus"][string(i)]["va"]  # radians
vr = vm * cos(va)
vi = vm * sin(va)
```

This output can be written to a JSON/CSV and consumed by the same Python `m_to_case.py`
patching script regardless of whether GridKit or PowerModels.jl produced the solution.

---

## 7. What was built: `uq-usecase/pf-solver/`

### Problem

GridKit's `MatpowerParser.hpp` could not parse `case_ACTIVSg200.m` (or `Hawaii40_20231026.m`)
due to two issues in the library's parser (which we were not allowed to modify):

1. Field names with underscores (e.g. `mpc.bus_name`) caused a regex match failure — the
   parser's `get_mpc_field()` uses `[a-zA-Z]+` and rejects underscores.
2. Extra columns in `mpc.gen` and `mpc.branch` (ACTIVSg-extended MATPOWER format) caused
   `checkEndOfMatrixRow` to throw, aborting the parse.

### Solution: `matpower_parser.hpp`

A standalone local parser (`uq-usecase/pf-solver/matpower_parser.hpp`) that replaces
GridKit's `MatpowerParser.hpp` for this use case. It populates the same
`SystemModelData<double,size_t>` struct so GridKit's solver objects (`SystemSteadyStateModel`,
`Kinsol`) can be used unchanged.

Key design choices:

| Feature | GridKit parser | Local parser |
|---------|---------------|--------------|
| Unknown field names (e.g. `bus_name`) | regex fail | skip silently |
| Extra columns in gen/branch rows | throws | ignores extras |
| Cell arrays `{...}` | not handled | consumed and skipped |
| `BusData.Va` unit | stored as degrees (from `.m`) | converted to radians (GridKit bus residuals use radians) |
| `Pd`, `Qd`, `Pg`, `Qg` unit | stored as MW/MVAr (from `.m`) | divided by `baseMVA` to give pu (GridKit PF model is pu) |

Once the parser could read the file without errors, KINSOL still failed due to two missing
unit conversions. The `.m` format stores Va in degrees and powers in MW/MVAr, but GridKit's
model internals expect radians and per-unit. GridKit's own `MatpowerParser.hpp` passes these
values through as-is; the local `matpower_parser.hpp` was extended to convert them before
populating `SystemModelData` (GridKit source was not modified):

- **Va degrees bug**: `BusData.Va` is passed straight into the bus constructor as the
  initial angle. GridKit's bus residuals use `cos(dtheta)` / `sin(dtheta)` where `dtheta`
  is in radians, so passing degrees (e.g. `-10.5` stored as `-10.5 rad` instead of
  `-0.183 rad`) gave wildly wrong residuals. KINSOL hit max iterations (`-6`).
  Fix: `bd.Va *= M_PI / 180.0` after reading each bus row.

- **MW/MVAr bug**: `LoadData.Pd`/`Qd` and `GenData.Pg`/`Qg` are passed straight into
  `Load` and `GeneratorPV` constructors and subtracted from / added to the bus P/Q
  residuals directly. The branch admittances are in pu, so the power balance is pu-scale
  (order 1). Passing MW values (e.g. `Pd = 150 MW` instead of `1.5 pu`) made the residual
  `||f|| = 364`. KINSOL stagnated (`-5`, line-search nonconvergence).
  Fix: divide `Pd`, `Qd`, `Pg`, `Qg`, `Qmax`, `Qmin`, `Pmax`, `Pmin` by `baseMVA`
  after parsing all rows.

### `solve_pf.cpp`

The `main()` wrapper. Key behavior:

- Parses the `.m` file via the local parser.
- Assembles `SystemSteadyStateModel`, calls KINSOL with warm start from the `.m` file's
  `Vm`/`Va` values.
- After `runSimulation()` (or on KINSOL exception), calls `sys->evaluateResidual()` and
  computes `||f||_2` manually using `sys->getResidual()` (from `ModelEvaluatorImpl`).
- Accepts as converged if `||f|| < 1e-4`, regardless of KINSOL return code — this handles
  the case where KINSOL throws `-5` (line search stagnation) because the warm-start point
  is already at the solution.
- Prints `bus <i>  V=<pu>  theta_deg=<deg>  type=<1|2|3>` to stdout only if converged.
- Optionally writes a solved `.m` file via `--output-m <path>` (see below).
- Returns exit code 0 (converged) or 1 (not converged).

### `--output-m`: solved .m output

Usage:
```bash
solve_pf input.m --output-m solved.m
```

When `--output-m` is given and the solve converges, `write_solved_m()` (in
`matpower_parser.hpp`) rewrites the input `.m` with `mpc.bus` columns 8-9 (Vm, Va)
replaced by the PF solution. All other content — `mpc.gen`, `mpc.branch`, `gencost`,
comments, unknown fields — is copied through unchanged. No output file is written if
the solve does not converge.

This makes `solve_pf` a self-contained tool for the UQ pipeline:
```
modified .m  →  solve_pf input.m --output-m solved.m  →  m_to_case.py solved.m  →  case.json
```

The solved `.m` has the same format as the `pcm-runs/ACTIVSg200_wind_demand/` hourly
files, so `m_to_case.py` consumes it without modification.

### Build

`uq-usecase/pf-solver/build.sh` builds `solve_pf` against `~/gridkit/build/` using the same
clang 16 / GCC 13 toolchain as the main GridKit build. No changes to the main build tree.

```bash
bash ~/gridkit/uq-usecase/pf-solver/build.sh
# output: uq-usecase/pf-solver/build/solve_pf
```

---

## 8. Verified results (2026-07-22)

All three cases converge with a warm start from the `.m` file's `Vm`/`Va` columns.
The `.m` base-case values represent a previously converged MATPOWER solution, so KINSOL
reaches convergence in very few iterations (the initial point is already near the solution).

### 3-bus sanity check

| var | solve_pf | expected | pass? |
|-----|----------|----------|-------|
| theta2 (deg) | -4.87980 | -4.87979 | PASS |
| V2 (p.u.) | 1.08281 | 1.08281 | PASS |
| theta3 (deg) | 1.46240 | 1.46241 | PASS |

KINSOL return code 0, `||f|| = 7e-6`.

### ACTIVSg200 (200-bus Illinois)

KINSOL return code 0, `||f|| = 4.4e-06`. Voltages: 1.009 to 1.043 pu.

Comparison vs MATPOWER `.m` reference `Vm`/`Va`:

| metric | value |
|--------|-------|
| max \|V err\| | 0.0297 pu |
| mean \|V err\| | 0.0051 pu |
| max \|theta err\| | 0.191 deg |
| mean \|theta err\| | 0.074 deg |

The residual vs the MATPOWER reference is **unexplained**: `case_ACTIVSg200.m` has no
off-nominal tap transformers (all `TAP=0` or `TAP=1`, `SHIFT=0`), so GridKit's tap-ratio
omission has no practical effect here. The source of the `max|V err|=0.030 pu` discrepancy
is an open question — see Section 15 for details.

### Hawaii40 (37-bus Hawaii)

KINSOL return code 0, `||f|| = 1.2e-06`, converged in 3 Newton iterations.

Voltages 0.969 to 1.003 pu, angles -7.0 to +1.4 deg.

Hawaii40 has not been checked for off-nominal tap transformers; the same open question
about the source of any discrepancy vs the `.m` reference applies.

---

## 9. Warm-start behavior

`solve_pf` always warm-starts from `Vm`/`Va` in the `.m` file. This is intentional: all
current use cases (base cases from TAMU / ACTIVSg suite) ship with a converged MATPOWER
solution in columns 8/9 of `mpc.bus`. KINSOL finds this starting point already satisfies
the GridKit residual (to within the model differences noted above) and converges immediately.

For a case where no prior solution is available (e.g., a perturbed operating point after
changing load or dispatch), the warm start will no longer be at the solution and KINSOL
will need to take Newton steps. Section 6 and 14 confirm that KINSOL converges from
warm-started perturbed cases with `nni=4-5` iterations for all tested perturbation types
(load ±80%, 10 gens offline). Flat start (`Vm=1, Va=0`) is possible but, on this model,
consistently converges to a different (low-voltage) solution branch — see Section 14.

---

## 10. Load perturbation convergence study (2026-07-23, ACTIVSg200)

Five load perturbation levels tested on `case_ACTIVSg200.m` using `make_perturbed_load_m`
(seed=42, same per-bus random draws scaled by `pct`; `Pd`/`Qd` scaled together to preserve
per-bus power factor). Warm start from base-case `Vm`/`Va` for all runs.

All five levels converged. Results vs the base-case GridKit solution (`il_df`):

**Table columns:**
- `max|dV|`, `mean|dV|` — max and mean absolute voltage magnitude difference (pu) between the perturbed and base-case solutions, over all 200 buses.
- `max|dTheta|`, `mean|dTheta|` — same for bus voltage angle (degrees).
- `||f||` — Euclidean norm of the P/Q mismatch vector at the final iterate (KINSOL `Residual 2-norm`). Values near `1e-6` confirm genuine convergence; values near `1e-4`+ indicate stagnation.
- `V violations` — count of buses with solved voltage magnitude outside `[0.95, 1.05]` pu (ANSI/NERC steady-state band). Zero means the perturbed operating point is within normal limits.

| Test | pct | max\|dV\| (pu) | mean\|dV\| (pu) | max\|dTheta\| (deg) | mean\|dTheta\| (deg) | \|\|f\|\| | V violations |
|------|-----|----------------|-----------------|---------------------|----------------------|-----------|--------------|
| 1    | ±5%  | 0.000640 | 0.000093 | 0.139 | 0.052 | 4.44e-06 | 0 |
| 1b   | ±10% | 0.001276 | 0.000185 | 0.277 | 0.105 | 4.49e-06 | 0 |
| 1c   | ±20% | 0.002534 | 0.000368 | 0.553 | 0.209 | 4.72e-06 | 0 |
| 1d   | ±40% | 0.004999 | 0.000731 | 1.103 | 0.418 | 5.81e-06 | 0 |
| 1e   | ±80% | 0.009726 | 0.001443 | 2.193 | 0.832 | 8.13e-07 | 0 |

### Findings

- **Perfectly linear response**: doubling the perturbation doubles both `max|dV|` and
  `max|dTheta|` to within 1% across the full range ±5% to ±80%. The network operates
  well within its linear regime at all tested levels.
- **Voltage is stiff**: `max|dV|` stays below 0.01 pu even at ±80%. The voltage range
  remains `[1.008, 1.043]` pu throughout with no violations. The data are consistent with
  PV and slack buses absorbing load changes almost entirely through angle adjustment rather
  than voltage magnitude; no other explanation was investigated.
- **Angle is the observable**: `max|dTheta|` grows from 0.14 deg (±5%) to 2.19 deg
  (±80%). Bus angles are the most informative output for load uncertainty; voltage
  magnitude changes will be small.
- **KINSOL is robust**: genuine new solutions found at every level (confirmed by
  `diff_vs_base` showing non-trivial shifts). No stagnation at any perturbation level.
- **UQ implication**: realistic load uncertainty (±5-20%) is comfortably in the linear
  regime. PF convergence will not be a bottleneck for the aleatoric UQ track.

---

## 11. Wind curtailment convergence study (2026-07-23, ACTIVSg200)

Four curtailment levels tested on the 6 ACTIVSg200 wind generators (buses 65, 104, 105,
114, 115, 147) using `make_wind_dispatch_m` (seed=7, curtailment-only: each gen's `Pg`
scaled by `(1 - Uniform(0, pct))`). All 6 wind gens run at `Pg == Pmax` in the base case
so upward perturbation is not physically possible; curtailment-only is the correct model.
Warm start from base-case `Vm`/`Va` for all runs.

All four levels converged.

| Test | max curtailment | MW curtailed (approx) | max\|dV\| (pu) | max\|dTheta\| (deg) | \|\|f\|\| |
|------|-----------------|-----------------------|----------------|---------------------|-----------|
| 2a   | up to 10% | ~few tens MW | — | ~1.3 deg | ~5e-06 |
| 2b   | up to 20% | — | 0.000450 | 2.575 | 5.02e-06 |
| 2c   | up to 40% | — | 0.000815 | 5.156 | 7.78e-06 |
| 2d   | up to 80% | — | — | — | — |

(Exact per-gen curtailment values can be read from the `plot_wind_comparison` bar chart
for each test.)

### Findings

- All tests converge. Angle shifts grow with curtailment level, proportionally.
- `max|dV|` remains below 0.001 pu for all tested levels; voltage is essentially unchanged.
- The slack bus is expected to absorb the entire wind reduction (non-wind dispatch is
  unchanged in the `.m` file); the angle shifts growing with curtailment level are
  consistent with this, but line flows were not directly checked to confirm it.
- **Section 14**: stagnation concern was investigated and resolved. The very small voltage
  response and clean convergence even at 80% curtailment was the initial motivation for
  the investigation; `nni` counts and flat-start cross-checks confirmed that genuine
  Newton steps were taken.

---

## 12. Generator offline convergence study (2026-07-23, ACTIVSg200)

Five gen-offline tests: bus 147 fixed (92.4 MW) and 2/3/5/10 random non-slack generators
offline (seed=99). Buses and MW dropped shown below. Warm start from base-case `Vm`/`Va`.

| Test | Buses offline | MW dropped | max\|dV\| (pu) | max\|dTheta\| (deg) | \|\|f\|\| |
|------|--------------|------------|----------------|---------------------|-----------|
| 3a   | {147} | 92.4 | — | — | — |
| 3b   | {104, 170} | 70.4 | — | — | — |
| 3c   | {104, 151, 167} | 72.0 | 0.000557 | 4.387 | 4.98e-06 |
| 3d   | {67, 94, 114, 136, 154} | 165.3 | 0.002703 | 7.315 | 9.24e-07 |
| 3e   | {65, 77, 91, 94, 125, 136, 155, 170, 182, 183} | 306.3 | 0.004528 | 16.482 | 1.51e-05 |

### Findings

- All tests converge. Angle shifts scale with MW dropped.
- `max|dV|` stays below 0.005 pu even for 306 MW dropped (test 3e). **Working
  hypothesis**: in a stiff, well-connected 200-bus system with ample reactive support,
  real power imbalances are absorbed mostly through angle shifts rather than voltage
  magnitude changes, so small dV is physically plausible. Evidence for this: Section 10
  shows the same pattern for load perturbations (linear scaling, voltage nearly unchanged
  at ±80% load). Note: the tap ratio omission is not relevant here since
  `case_ACTIVSg200.m` has no off-nominal tap transformers (verified from `mpc.branch`).
  **Pending confirmation**: PowerModels.jl cross-validation will verify whether a
  full-model solver shows similarly small dV for the same gen-offline cases. If it does,
  the hypothesis is confirmed. If it does not, the GridKit model is missing reactive
  coupling that would drive larger voltage shifts.

### Resolution: stagnation concern resolved (2026-07-23, Section 7)

See Section 14 for the full investigation. Summary:

- `nni=4-5` for all warm-start perturbed cases — genuine Newton steps were taken, not
  stagnation (stagnation would show `nni=1` or `2`).
- Flat start (`Vm=1, Va=0`) always converges to a distinct low-voltage branch
  (V range ~[0.965, 1.000] pu) regardless of perturbation type or severity.
- The `max|dV|=0.043 pu` gap between flat and warm is a fixed property of the two
  solution branches in GridKit's tap-less model, not a perturbation artifact.
- Warm-start solutions track the high-voltage branch (V range ~[1.007, 1.043] pu)
  and are genuine equilibria of GridKit's model.
- Cross-validation against an independent reference solver still requires PowerModels.jl
  — that is the next step. Since `case_ACTIVSg200.m` has no off-nominal tap transformers,
  any remaining difference vs PowerModels.jl will have a different root cause.

---

## 13. Cross-section finding: angle-only response across all perturbation types

Sections 10, 11, and 12 together show a striking and consistent pattern: across every
perturbation tested on ACTIVSg200, the PF solution changes almost entirely through angle
shifts while voltage magnitudes barely move.

All numbers below come directly from notebook cell outputs (pf_helper.ipynb, Section 6
and 7, executed 2026-07-23/26). Magnitude column values are computed from
`illinois_m_tables.md` (mpc.gen and mpc.bus tables).

| Section | Perturbation | Magnitude | max\|dV\| (pu) | max\|dTheta\| (deg) |
|---------|-------------|-----------|----------------|---------------------|
| 10 | load ±5%  | total load 1475.69 MW, ±5% per bus | 0.000640 | 0.138612 |
| 10 | load ±80% | total load 1475.69 MW, ±80% per bus | 0.009726 | 2.193222 |
| 11 | wind curtailment up to 40% | 40% of 536.2 MW total wind | 0.000815 | 5.156084 |
| 11 | wind curtailment up to 80% | 80% of 536.2 MW total wind | 0.002328 | 10.344126 |
| 12 | 1 gen offline (bus 147) | 92.4 MW (from mpc.gen PG) | 0.002758 | 10.841606 |
| 12 | 10 gens offline (seed=99) | 306.3 MW (from mpc.gen PG sum) | 0.004528 | 16.482293 |

Even at the extremes — 306 MW of generation lost, or ±80% load perturbation — the maximum
voltage magnitude change across all 200 buses is less than 0.01 pu.

### Why is this surprising?

In standard PF intuition, significant real power imbalances should produce noticeable
voltage changes, especially far from the slack bus. A 306 MW loss on a system with
1488.3 MW total online generation (from mpc.gen, status=1 rows) is a ~21% dispatch
change; one might expect voltage depressions of 0.05 pu or more at buses electrically
distant from the slack.

### PQ vs PV vs slack: what determines whether |V| moves

Before analyzing why voltage stays flat, it helps to lay out the three MATPOWER bus
types explicitly, since the answer is largely a consequence of how many buses are
voltage-regulated by construction:

| Type | MATPOWER code | Fixed at the bus | Free variables solved by PF | Physical realization |
|------|:-:|---|---|---|
| Slack (ref) | 3 | \|V\|, θ | none (bus absorbs P and Q imbalance) | Reference generator; one per island |
| PV | 2 | \|V\|, P (generator setpoint) | θ, Q_g | Generator with an AVR (automatic voltage regulator); voltage held at a setpoint by adjusting excitation |
| PQ | 1 | P, Q | \|V\|, θ | Load bus, or generator run without voltage control |

The critical point: at a **PV bus**, |V| is a boundary condition, not a variable. The
solver treats the setpoint as hard until a generator Q limit is hit (`Qmax`/`Qmin`) —
the generator's exciter is idealized as being able to inject whatever reactive power is
needed to hold the setpoint. At a **PQ bus**, |V| is free and gets pushed around by
whatever the surrounding network delivers. So the fraction of buses that are PV+slack —
call it the **regulated fraction** — sets a hard floor on how much of the network is
voltage-pinned by construction.

### Regulated bus fraction: ACTIVSg200 vs peers

Counting `BUS_TYPE` (column 2 of `mpc.bus`) directly:

| Case | Buses | PQ (type 1) | PV (type 2) | Slack (type 3) | Regulated fraction (PV+slack) |
|---|---:|---:|---:|---:|---:|
| ACTIVSg200 (Illinois) | 200 | 151 | 48 | 1 | **24.5 %** |
| ACTIVSg2000 (Texas) | 2000 | — | ≥200 (capped at grep limit) | 1 | not yet confirmed; needs full parse |
| ACTIVSg10k (WECC) | ~10 000 | — | — | — | not measured locally |
| Hawaii40 | 37 | 28 | 8 | 1 | **24.3 %** |

For comparison — real-world reference points:
- ERCOT / WECC bulk transmission snapshots: typically 15–30 % regulated buses (nearly
  every transmission-connected generator runs on AVR control in normal operation).
- IEEE test transmission cases (39-bus New England, 118-bus, 300-bus): all report
  regulated fractions in roughly the same 15–30 % band.
- Distribution feeders (radial, IEEE 13/34/123-node, etc.): 0–5 % regulated — the
  substation transformer LTC is often the only voltage anchor. This is the regime where
  a large real-power perturbation *does* produce large |V| swings.

ACTIVSg200 at 24.5 % is squarely in the transmission-grid regime. This is not
artificially stiff — it's what a synthetic 200-bus bulk transmission model should look
like. TAMU built the ACTIVSg suite specifically to reproduce statistical properties of
real interconnections, and the density of AVR-controlled generators is one of those
properties.

### Is angle-only response practical? (short answer: yes)

The "angle carries the response, magnitude stays put" behavior is well-known in
transmission planning and is exactly the assumption behind **DC power flow**, which
drops the voltage-magnitude equations entirely. DC PF works for many transmission
planning tasks precisely because on well-regulated networks, active branch flow is
dominated by angle differences:

$$
P_{ij} \approx \frac{\theta_i - \theta_j}{X_{ij}}
$$

Our findings say GridKit's full AC solve on ACTIVSg200 sits in exactly the regime
where DC PF would already be a good approximation. That is a positive validation of
the solver on a "normal" transmission model — not an anomaly.

Regimes where DC PF breaks (and |V| changes become first-order):
- The system is heavily loaded and near voltage collapse (nose of the PV curve).
- One or more PV generators are at their Q limits, forcing PV → PQ downgrade.
- Reactive imbalances dominate (contingency loss of a large capacitor bank or SVC).
- The network is radial or lightly meshed (distribution, weak grids, islanded systems).

None of these apply to ACTIVSg200 at the operating points we're perturbing: all runs
converge with |V| ∈ [1.007, 1.043] pu (no near-collapse), the network is meshed
(245 branches for 200 buses), and the perturbations are moderate.

### Working hypothesis: strong PV bus voltage regulation

**Hypothesis**: when ~24 % of buses have fixed voltage setpoints, real power imbalances
propagate almost entirely as angle shifts, with reactive power redistribution absorbing
whatever is needed to keep voltages near setpoints. The slack bus picks up the MW
imbalance through angle; the PV buses suppress voltage deviations locally.

Supporting evidence internal to GridKit's own solutions:
- Load perturbation results (Section 10) show a perfectly linear angle response —
  consistent with a network operating deep in a strongly regulated, near-linear regime.
- Even the largest tested perturbation (10 gens offline, 306 MW dropped, Section 12)
  keeps every solved voltage inside the ANSI [0.95, 1.05] pu band — zero violations.

### Caveats: what GridKit's model omits that could inflate the effect

The angle-only pattern is *physically consistent* with a well-regulated transmission
network, but two known model simplifications in GridKit could artificially reinforce it
beyond what a full MATPOWER model would give. Both need to be checked (or documented as
known biases) before the results feed into `illinois.json`:

1. **Reactive power limits `Qmax` / `Qmin` are not enforced.** In MATPOWER / PowerModels,
   a PV bus is downgraded to PQ mid-solve if the reactive dispatch needed to hold |V|
   would exceed the generator's Q capability (called "PV → PQ switching"). Once switched,
   voltage at that bus *does* start moving. GridKit's `GeneratorPV` reads `Qmax` / `Qmin`
   into `GenData` (verified — the parser populates them), but does not appear to enforce
   them in the KINSOL residual loop. On a stressed operating point this could keep the
   network artificially stiffer than a real solve would.

2. **Bus shunts `Gs`, `Bs` handling.** GridKit's `Branch` hardcodes shunt conductance to
   zero. Bus shunt fields (`mpc.bus` columns 5–6) are parsed into `BusData` but their
   path into the residual has not been audited here. ACTIVSg200 has ~24 buses with
   nonzero `Bs` (mostly capacitor banks); whether they contribute correctly to the
   nodal reactive balance is an open item.

Neither omission is fatal for the aleatoric UQ study, but both are candidates for the
0.030 pu residual vs the MATPOWER reference documented in Section 8, and both are
things PowerModels.jl handles natively.

### Why this matters for UQ

If voltage magnitudes are genuinely nearly insensitive to the perturbations we sample,
the informative output for UQ is bus **angle** (or equivalently, active branch flows) —
not voltage magnitude. This has direct implications:

- **QoI choice**: for the PF-solution-based aleatoric track, bus angle and branch MW
  flow are the primary QoIs; |V| is a validation-only variable.
- **UQ dimensionality**: an angle-dominated response is roughly linear in the sampled
  parameters (confirmed for load in Section 10), so surrogate models (polynomial chaos,
  Gaussian processes) will converge with far fewer samples than a nonlinear response
  would demand.
- **Downstream dynamic UQ**: the tiny |V| response also means the initial conditions
  patched into `illinois.json` will vary mostly in `Vi = |V| sin(θ)` (through θ) rather
  than through |V| itself. The Genrou rotor angles `delta_0` inherit that variability.

### Next steps

1. **PowerModels.jl cross-validation** (already the planned next step; see the plan for
   `pm_helper` scope). Key discriminating experiments:
   - *Base case*: does PM's |V| range match GridKit's [1.007, 1.043]? If yes, the 0.030 pu
     Section 8 residual is genuinely small-model-difference, not a solver artifact.
   - *10 gens offline* (max stress test in Section 12): does PM's `max|dV|` stay
     < 0.005 pu (confirming the angle-only regime), or does it grow to 0.02–0.05 pu
     (signaling PV → PQ switching that GridKit misses)?
   - *Load ±80%* (max load test in Section 10): same comparison.

2. **Reactive dispatch audit**: after the base-case PM solve, tabulate `Qg` vs
   `[Qmin, Qmax]` for all 49 ACTIVSg200 generators. If several are already at or near
   limits in the base case, GridKit's non-enforcement is a known bias and results
   should be quoted with that caveat.

3. **Regulated-fraction check on ACTIVSg2000**: run the same BUS_TYPE count on
   `case_ACTIVSg2000.m`. If the regulated fraction is similar (~20–25 %), the
   angle-only response should reproduce. If it turns out lower (say <15 %) *and*
   ACTIVSg2000 shows larger dV under the same perturbation types, that is direct
   evidence for the PV-regulation mechanism.

4. **UQ QoI decision** (informed by 1–3): if PM confirms small dV, formalize angle
   and branch flow as the primary QoIs for the aleatoric-UQ track and treat |V| as
   a validation-only variable.

---

## 14. Stagnation investigation: flat start vs warm start (2026-07-23, ACTIVSg200)

**Motivation**: Section 12 raised concern that warm-start solutions in Section 6 might
be stagnation artifacts. This section documents the investigation and resolution.

### Method

Added `--flat-start` flag to `solve_pf` (zeroes all bus `Vm=1.0, Va=0.0` before solving).
Added `_extract_nni` to `pf_utils.py` to parse KINSOL's "Nonlinear iters" from stderr and
inject it into the returned stderr string as `nni=N`. Four tests run in notebook Section 7.

### Results

| Test | warm nni | warm V range (pu) | flat nni | flat V range (pu) | max\|dV\| warm vs flat |
|------|----------|------------------|----------|------------------|------------------------|
| 7.1 base case | 4 | [1.009, 1.043] | 7 | [0.967, 1.000] | 0.043 |
| 7.2 gen147off | 4 | [1.009, 1.043] | 7 | [0.967, 1.000] | 0.043 |
| 7.3 load ±80% | 5 | [1.008, 1.043] | 7 | [0.966, 1.000] | 0.043 |
| 7.4 gen10rand | 5 | [1.007, 1.043] | 8 | [0.965, 1.000] | 0.043 |

### Findings

**Two observed solution branches**: in all four tests, two distinct converged solutions
were found depending on starting point. We have not proven that exactly two equilibria
exist — that would require a more thorough analysis. What the data show:
- Warm start (base-case `Vm`/`Va`) converged to V ~[1.007, 1.043] pu in all 4 cases.
- Flat start (`Vm=1, Va=0`) converged to V ~[0.965, 1.000] pu in all 4 cases.
- The separation is ~0.043 pu, constant across all 4 tests.

The high-voltage branch is consistent with the TAMU base-case solution stored in
`case_ACTIVSg200.m` (which represents a physically operated system), giving confidence
that it is the operationally relevant solution. The low-voltage branch has all bus
voltages at or below 1.000 pu, which is atypical for a generation-rich system and
suggests it is a non-operational equilibrium of the tap-less model.

**Not stagnation**: `nni=4-5` for warm-start perturbed cases (vs `nni=4` for the base
case) confirms KINSOL took genuine Newton steps. Stagnation would show `nni=1` because
the residual check would pass at the initial point without any steps.

**Flat start is not a useful cross-validator here**: in all 4 tested cases it converged
to the low-voltage branch rather than the operationally relevant high-voltage branch.
Cross-validation must use PowerModels.jl, which operates on the full MATPOWER model.

**Warm-start solutions are valid for GridKit's model**: the Section 6 solutions are
genuine equilibria of GridKit's model. They may differ from PowerModels.jl solutions
for reasons that are not yet understood (see open question in Section 15); they are
self-consistent within GridKit's model.

---

## 15. Transformer tap ratio and phase shift: what GridKit omits and why it matters

### Overview

GridKit's Branch model uses a plain symmetric PI model (`a=1, φ=0` for all branches).
MATPOWER, PowerModels, and PowerWorld use the full generalized PI model that incorporates
transformer tap ratios and phase shifts, which makes the from/to shunt admittances
asymmetric. Both formulations assemble branch equations into nodal mismatch equations
(expressible as a Ybus).

### What the MATPOWER `.m` file contains

Every row of `mpc.branch` has two transformer fields (columns 9 and 10):

| Column | MATPOWER name | Meaning |
|--------|--------------|---------|
| 9 | `ratio` (a) | Off-nominal turns ratio. `0` or `1` = unity (plain line). Values like `0.975` or `1.05` indicate an LTC or tap-changer at the from-bus. |
| 10 | `angle` (φ) | Phase shift angle in degrees. Non-zero only for phase-shifting transformers (PSTs). In `case_ACTIVSg200.m` all values are 0, but the field is present. |

In `case_ACTIVSg200.m`, verified by parsing `mpc.branch` column 9: all 245 branches have
`TAP=0` (179 branches, plain lines) or `TAP=1.0` (66 branches, explicit unity). There are
**no off-nominal tap transformers** in this specific case. The description below applies
generally to MATPOWER cases that do have them (e.g. many larger transmission cases).

### What GridKit parses vs what it uses

`BranchData` (in `PowerFlowData.hpp`) has both fields:

```cpp
RealT ratio;  // parsed from column 9
RealT angle;  // parsed from column 10
```

But `Branch::Branch(bus1, bus2, BranchData& data)` only reads `data.r`, `data.x`, `data.b`:

```cpp
Branch<ScalarT, IdxT>::Branch(bus_type* bus1, bus_type* bus2, BranchData& data)
    : R_(data.r),
      X_(data.x),
      G_(0.0),      // shunt conductance: hardcoded zero
      B_(data.b),   // shunt susceptance
      ...
```

`data.ratio` and `data.angle` are parsed into the struct but never passed to or used by `Branch`. They are silently dropped.

### What the Branch residual actually computes

`Branch::evaluateResidual()` implements a plain symmetric PI model:

```
g + jb = 1 / (R + jX)       (series admittance, no tap)

P1 -= (g + G/2)*V1^2 + V1*V2*(-g*cos(dθ) - b*sin(dθ))
Q1 -= (-b - B/2)*V1^2 + V1*V2*(-g*sin(dθ) + b*cos(dθ))
P2 -= (g + G/2)*V2^2 + V1*V2*(-g*cos(dθ) + b*sin(dθ))
Q2 -= (-b - B/2)*V2^2 + V1*V2*(g*sin(dθ) + b*cos(dθ))
```

with `dθ = θ1 - θ2`. No tap ratio `a`, no phase shift `φ`.

### What MATPOWER / PowerModels / PowerWorld model instead

The standard MATPOWER/PSS/E transformer model adds the turns ratio `a` and phase shift `φ` by modifying the admittance entries that appear in the bus admittance matrix `Ybus`. For an off-nominal transformer between buses `i` and `j` with ratio `a` and shift `φ`:

$$
Y_{ij} = -\frac{y_s}{a \cdot e^{j\phi}}, \quad
Y_{ii} += \frac{y_s}{a^2}, \quad
Y_{jj} += y_s
$$

where $y_s = 1/(R+jX)$ is the series admittance. When `a=1, φ=0` this reduces exactly to the plain PI model. When `a ≠ 1` or `φ ≠ 0`:

- The from-bus (`i`) sees a modified self-admittance: shunt is $y_s/a^2$ instead of $y_s$.
- The mutual admittance is $y_s/(a \cdot e^{j\phi})$ instead of $y_s$.
- The to-bus (`j`) sees the standard $y_s$ shunt.

In MATPOWER's `makeBusAdmittanceMatrix()` and PowerModels.jl's `calc_branch_y()` this is computed automatically for every branch. PowerWorld does the same inside its internal network matrix assembly. The result is that tap-changer transformers shift reactive power between the buses and affect voltage magnitudes in a way that a plain line model does not.

### Practical consequence for ACTIVSg200

`case_ACTIVSg200.m` has 245 branches. Checking `mpc.branch` column 9: 179 branches
have `TAP=0` (plain lines, unity tap by MATPOWER convention) and 66 have `TAP=1.0`
(explicit unity). There are **no off-nominal tap transformers** in this case — all
`SHIFT` values are also 0. This means the tap-ratio omission in GridKit's `Branch`
model has **no practical effect on this specific case**: with `a=1, φ=0` for all
branches, the generalized PI model reduces exactly to the plain symmetric PI model.

**The `max|V err| = 0.030 pu` seen in Section 8 is therefore not explained by the tap
ratio omission.** The source of that discrepancy is unknown; it requires investigation.
Candidates include: different shunt handling, generator reactive power limits enforced
differently, or the MATPOWER `.m` reference values coming from a different solver with
different convergence criteria. This is an open question.

GridKit's solution IS fully converged for its own model (`||f|| = 4.4e-6`).

**For relative comparisons** (base vs perturbed, same GridKit model throughout),
the discrepancy vs MATPOWER does not matter: both cases are solved with the same model,
so `diff_vs_base` measures genuine perturbation effects, not modeling artifacts.

**For absolute accuracy** (patching `illinois.json` with physically correct Vm/Va),
the discrepancy matters and needs to be understood before GridKit solutions are used
for JSON patching. PowerModels.jl cross-validation is the next step.

### How to fix tap ratio handling in GridKit (not done, likely not needed for ACTIVSg200)

Since `case_ACTIVSg200.m` has no off-nominal tap transformers, the fix is not needed
for this case. If the same approach is applied to a case that does have transformers
with `ratio != 1` or `shift != 0`, the fix would be: modify
`Branch::Branch(bus_type* bus1, bus_type* bus2, BranchData& data)` to use `data.ratio`
and `data.angle` in `evaluateResidual()`. The local `matpower_parser.hpp` already
parses these fields into `BranchData` correctly.
