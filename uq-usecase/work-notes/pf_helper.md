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
| 3 | Does `Branch` handle non-unity transformer ratio (`ratio != 0`)? | Check `Branch` source; ACTIVSg200 has transformers with `ratio=0.975` etc. |
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
