# Power Flow Solutions for ACTIVSg200: GridKit and PowerModels.jl

This document describes power flow (PF) solution experiments on the 200-bus Illinois
synthetic case (`case_ACTIVSg200.m`) using two solvers:

- **GridKit `solve_pf`** — a Newton/KINSOL-based solver built from GridKit's `PowerFlow`
  module (`uq-usecase/pf-solver/`). Results in `pf_helper.ipynb`.
- **PowerModels.jl (PM.jl)** — a Julia/Ipopt industry-standard MATPOWER-equivalent solver.
  Results in `pm_helper.ipynb`.

**Key findings:**
1. Both solvers converge on the base case and all 14 perturbed cases (±80% load, 80%
   wind curtailment, 10 generators offline).
2. **GridKit `solve_pf` solves the wrong physical problem.** It has no PV→PQ switching:
   every generator is assumed to have unlimited reactive capability, so every PV bus
   voltage is pinned at its setpoint unconditionally regardless of generator Q limits.
   PM.jl (and all industry solvers) check Qmax/Qmin after each Newton step and let
   voltages float when a generator saturates. Verified by source inspection (Section 1)
   and behavioral fingerprint (Section 11). GridKit's solution is numerically converged
   but physically wrong whenever any generator is Q-limited.
3. **"Safe envelope" = the operating conditions under which no generator hits a Q limit
   in the real (PM.jl) solution.** Empirically for ACTIVSg200: ±80% random per-bus
   load perturbations around the base case (0.5x–1.5x total system load, N≤2 largest
   gens off). Inside this envelope the same generators are Q-limited in PM.jl as at
   base case, the GridKit voltage error is a constant ~0.030 pu offset that cancels in
   ΔV comparisons, and the response is angle-dominated (max|dV| < 0.01 pu — Section 5).
   Outside it, additional generators saturate in PM.jl, voltages sag, and GridKit's
   error grows sharply: 0.089 pu at 2.0x load, 0.279 pu at N=25 largest gens off, with
   qualitatively different violation counts (Section 11).
4. **When to use which solver:**
   - **GridKit**: only inside the safe envelope above, and only for angle/flow QoIs
     (not absolute voltages). Fast angle-only sanity checks during development.
   - **PM.jl**: everything else — absolute voltages, any operating point outside the
     safe envelope, Vr/Vi patching into case JSON, voltage constraint checking.
   - Using GridKit outside the safe envelope produces a physically inconsistent initial
     condition that is silently wrong: dynamic simulations start from a state that
     violates generator Q limits, biasing all post-fault metrics without any error signal
     (see Section 8 for the downstream consequences).
5. Within ±80% per-bus random load perturbation, the PF response is angle-dominated
   (`max|dV| < 0.01 pu`). Property of the network's 24.5% regulated-bus fraction
   (Section 5, confirmed by both solvers).
6. The voltage-stiff regime breaks at ~1.5x uniform system load; PV-curve nose between
   2x–3x base load. Gen outages exit the voltage-stiff regime faster per MW than uniform
   load increases (Section 9–10).
7. Two solution branches (warm-start high-V physical; flat-start low-V non-physical),
   ~0.044 pu gap. Gap grows with gen outages as reactive support is removed (Section 10).
8. **Verdict for the aleatoric UQ pipeline:** use PM.jl for all PF solves that feed
   initial conditions into phasor-dynamic runs. GridKit `solve_pf` is safe only inside
   the envelope defined in finding 3 above, and only for relative angle/flow comparisons.

   **What goes wrong if GridKit PF is used outside the safe envelope:**
   The PF solution becomes the initial condition (Vm, Va, and generator state) patched
   into the case JSON before the phasor-dynamic simulation. If that initial condition
   is wrong, the dynamic simulation starts from a physically impossible steady state:

   - **Wrong bus voltages at t=0.** GridKit keeps every PV bus clamped at its voltage
     setpoint even when the generator would physically have saturated its AVR. The
     patched Vm values are then 0.05–0.28 pu too high at many buses (Section 11).
     The dynamic solver begins from this artificially inflated voltage, producing a
     transient artifact at t=0 as the dynamic model tries to relax to a real
     equilibrium. For fault studies this corrupts the pre-fault steady-state baseline
     and makes post-fault metrics (nadir, settling time, first-swing amplitude)
     unreliable.

   - **Wrong generator reactive output at t=0.** Because GridKit never computed a
     physically consistent Q_g for over- or under-excited generators, their AVR
     initial conditions are wrong. In phasor-dynamic simulations the AVR/exciter
     starts integrated from those initial values; an inconsistent IC causes a
     non-zero drift even before the fault is applied. This appears as a spurious
     pre-fault transient that is indistinguishable from the fault response.

   - **Silently wrong, not loudly wrong.** The dynamic solver will not crash or return
     an error. The output trajectories look plausible; there is no in-band signal that
     the initial condition was infeasible. The only way to detect the contamination is
     to compare against a run with a correct IC (PM.jl), which is precisely the check
     in Section 11. For UQ, this means sensitivity indices (Sobol, PCE coefficients)
     computed from a contaminated ensemble are quietly biased toward zero for any
     generator whose Q limit was active in the true operating point, making H
     variations appear less detectable than they actually are.

   **In short:** using GridKit PF outside the envelope does not cause a crash; it
   silently poisons the ensemble initial conditions, biasing sensitivity estimates
   and potentially producing a spurious pre-fault drift in every run.

### TL;DR: why prefer PM.jl over GridKit `solve_pf` (unmodified) as the production PF solver

1. **PM.jl reproduces the TAMU/PowerWorld reference to numerical noise** (max|dV| = 1.1e-5 pu,
   max|dTheta| = 9.3e-4 deg, Section 4). GridKit deviates from both by ~0.030 pu on the
   same base case; the deviation is entirely on GridKit's side.
2. **GridKit has no PV→PQ switching** (Section 1, verified across six independent code
   paths and behaviorally in Section 11). This produces qualitatively wrong voltage
   violation counts outside a narrow envelope: 0 vs 44-85 over-voltage buses at 0.1x-0.5x
   load, 0 vs 85-177 under-voltage buses at N=10-25 largest gens off.
3. **The ~0.030 pu offset is only "constant" in a narrow band** (0.5x-1.5x uniform load,
   N≤2 largest gens off). Outside that band max|dV| grows to 0.089 pu at 2.0x load and
   0.279 pu at N=25 gens off (Section 11). The ΔV-cancellation argument (Section 7) does
   not survive stress.
4. **PM.jl handles tap ratios, phase shifters, and other MATPOWER features natively**
   (Section 14). Irrelevant for ACTIVSg200 but blocks any future extension to ACTIVSg2000
   or real-utility cases without additional GridKit development.
5. **PM.jl per-call cost is manageable**: ~8.5 s (5.5 s Julia startup + 3 s Ipopt solve).
   For 8760 hourly scenarios a batched-manifest mode (one Julia process loops over all
   cases) amortizes the startup cost; Slurm sbatch parallelization gives near-linear
   speedup.

**Bottom line**: PM.jl is an actively maintained, industry-validated, MATPOWER-equivalent
solver that requires no GridKit development to be correct on this network. GridKit
`solve_pf` in its current form is a fast angle-and-flow sanity checker, not a production
voltage solver. Section 16 outlines a concrete development plan for closing the gap on
the GridKit side if that path is chosen later.

---

## 1. GridKit PF solver: architecture and implementation

GridKit implements a Newton-based PF solver using SUNDIALS KINSOL. The model is assembled
from component objects (Bus, Branch, Generator, Load), each contributing residual and
Jacobian entries. The system is solved as $\mathbf{f}(\mathbf{x}) = 0$ for bus voltages
$(V_i, \theta_i)$.

Key source files:

| File | Role |
|---|---|
| `GridKit/Model/PowerFlow/SystemModelPowerFlow.hpp` | Assembles bus/gen/branch/load into a solvable system; calls KINSOL |
| `GridKit/Model/PowerFlow/MatpowerParser.hpp` | Parses MATPOWER v2 `.m` files |
| `GridKit/Model/PowerFlow/PowerFlowData.hpp` | `BusData`, `GenData`, `BranchData`, `LoadData` structs |
| `GridKit/Model/PowerFlow/Bus/` | `BusFactory` — creates slack, PQ, PV bus objects |
| `GridKit/Model/PowerFlow/Generator/` | `GeneratorFactory` — creates slack and PV gen models |
| `GridKit/Model/PowerFlow/Branch/` | `Branch` — symmetric PI-model (no tap ratio) |
| `GridKit/Solver/SteadyState/Kinsol.hpp` | KINSOL wrapper |

### Known model limitations (relevant to ACTIVSg200)

1. **No Q-limit enforcement / no PV→PQ switching (verified across six code paths).**
   All PV buses hold their voltage setpoint unconditionally. This is the cause of the
   ~0.030 pu voltage bias vs industry solvers (Section 6) and of the qualitatively
   wrong violation counts under stress (Section 11). Direct source inspection:

   1. **Grep across the entire `Model/` and `Solver/` trees**:
      `grep -rniE "qmax|qmin|pv[_ ]?to[_ ]?pq|q[_ ]?limit|typeSwitch" Model/ Solver/`
      returns only five hits, all in `MatpowerParser.hpp` (parsing the column) and
      `PowerFlowData.hpp` (`GenData::Qmax`, `GenData::Qmin` struct fields, and their
      diagnostic pretty-printer). Zero hits in any `Bus/`, `Generator/`, `Branch/`,
      or solver source file.
   2. **`Model/PowerFlow/Bus/BusPV.cpp`** (145 lines): the PV-bus residual is
      `P() = 0.0; Q() = 0.0;` (reset only; injections are added by component models).
      The class stores only `V_` (setpoint) and `theta0_` (angle initial guess). No
      Qmax/Qmin, no switching state, no bus-type flip mechanism.
   3. **`Model/PowerFlow/Generator/GeneratorPV.cpp`** (109 lines): the reactive power
      member `Q_(data.Qg)` is commented out at construction and the residual
      contribution `bus_->Q() += Q_` is commented out. Only `bus_->P() += P_` remains.
      A PV generator never touches the bus's reactive residual, so there is nothing
      to compare against Qmax/Qmin.
   4. **`Model/PowerFlow/Bus/BusFactory.hpp`**: dispatches on `data.type` once at
      construction via `switch(data.type){case 1: new BusPQ; case 2: new BusPV;
      case 3: new BusSlack;}`, returning a raw pointer to a concrete subclass. No
      mechanism to swap the object type at runtime; bus type is fixed for the
      object's lifetime.
   5. **`Solver/SteadyState/Kinsol.cpp:121-127`** (`runSimulation()`):
      `retval = KINSol(solver_, yy_, KIN_LINESEARCH, scale_, scale_);` called
      exactly once. No wrapper loop.
   6. **`uq-usecase/pf-solver/solve_pf.cpp:97`** (the driver that produces the
      `solve_pf` binary): `ret = kinsol->runSimulation();` called exactly once.
      The local parser (`matpower_parser.hpp:284-285`) reads `Qmax`/`Qmin` and
      unit-scales them, but nothing downstream ever consults the values.

   All industry solvers (MATPOWER, PM.jl, PSS/E, PowerWorld) implement PV→PQ
   switching as an outer loop: after each Newton solve, check `Q_g` at each PV bus,
   flip to PQ (fix `Q_g = Qmax` or `Qmin`, let `|V|` float) if the limit is
   exceeded, and re-solve. GridKit has none of these ingredients: no per-Newton-step
   Q evaluation at PV buses, no bus-type flip mechanism, no outer switching loop
   around KINSOL. Behavioral confirmation (Section 11): PM.jl reports 44-85
   over-voltage violations at 0.1x-0.5x load (Qmin hits missed by GridKit) and
   85-177 under-voltage violations at N≥10 largest gens off (Qmax hits missed by
   GridKit); both directions match the "no switching" fingerprint.
   GridKit's own `Model/PowerFlow/Bus/README.md` acknowledges PV→PQ switching is
   what a real PF should do; the code does not implement it.
2. **Tap ratio/phase shift dropped.** `Branch` uses a plain symmetric PI model (`a=1,
   φ=0`). Fields are parsed but not used. Irrelevant for ACTIVSg200 (no off-nominal
   transformers), but would matter for cases with LTCs.

---

## 2. `solve_pf` binary: what was built

### Problem

GridKit's library parser could not read `case_ACTIVSg200.m` due to underscore-containing
field names (`mpc.bus_name`) and extra columns in `mpc.gen`/`mpc.branch`.

### Solution

A local parser (`pf-solver/matpower_parser.hpp`) populates the same `SystemModelData`
struct, handling unknown fields and extra columns gracefully. Two critical unit conversions
were added (degrees→radians for Va, MW→pu for Pd/Qd/Pg/Qg) that GridKit's internal
parser does not perform.

### `solve_pf` behavior

- Reads `.m` from `argv[1]`, assembles the system, runs KINSOL (warm-start from `.m` Vm/Va)
- Prints `bus <i> V=<pu> theta_deg=<deg> type=<1|2|3>` to stdout
- Accepts convergence if `||f|| < 1e-4` (handles stagnation at warm-start point)
- `--output-m <path>`: writes solved `.m` with updated Vm/Va columns
- `--flat-start`: resets all buses to Vm=1, Va=0 before solving
- Exit 0 = converged, 1 = not converged

Build: `bash ~/gridkit/uq-usecase/pf-solver/build.sh`

---

## 3. PowerModels.jl (PM.jl): setup

PM.jl is an industry-standard Julia package for AC/DC power flow and OPF. It reads
MATPOWER `.m` files natively, enforces generator Q limits via PV→PQ switching, handles
tap ratios and shunts, and uses Ipopt as the NLP solver.

`pm_solve.jl` (in `uq-usecase/pm-solver/`) is a thin Julia CLI driver:
```
julia112 --project=pm-solver/ pm_solve.jl <input.m> [--output-m ...] [--pf-type ac|dc] [--tol ...] [--max-iter ...] [--flat-start]
```

Packages pinned: PowerModels v0.21.6, Ipopt v1.15.0, JuMP v1.31.1.

Per-call overhead: ~8.5 s (5.5 s Julia startup + 3 s Ipopt solve). For 8760 scenarios
a batched manifest mode or Slurm sbatch parallelization is planned.

---

## 4. Base-case results

Both solvers converge on the base case with a warm start from the `.m` file's Vm/Va.

### GridKit

KINSOL return 0, `||f|| = 4.4e-6`, nni=4. V range: [1.009, 1.043] pu.

### PM.jl

Ipopt 4 iterations, `LOCALLY_SOLVED`. V range: [1.010, 1.056] pu.

### Three-way comparison (PM.jl, GridKit, TAMU/PowerWorld reference)

| pair | max \|dV\| (pu) | max \|dTheta\| (deg) |
|---|---|---|
| PM.jl vs TAMU/PowerWorld | 0.000011 | 0.000926 |
| PM.jl vs GridKit | 0.029728 | 0.190703 |
| GridKit vs TAMU/PowerWorld | 0.029730 | 0.190720 |

PM.jl reproduces the TAMU/PowerWorld solution to numerical noise. GridKit deviates from
both by ~0.030 pu, confirming the deviation is entirely in GridKit's model (Q-limit
non-enforcement), not between the two independent reference solvers.

---

## 5. Perturbation studies: angle-dominated response (GridKit)

### 5.1 Load perturbation (±5% to ±80%, random per-bus, seed=42)

| pct | max\|dV\| (pu) | max\|dTheta\| (deg) | \|\|f\|\| | V violations |
|-----|----------------|---------------------|-----------|--------------|
| ±5%  | 0.000640 | 0.139 | 4.44e-06 | 0 |
| ±10% | 0.001276 | 0.277 | 4.49e-06 | 0 |
| ±20% | 0.002534 | 0.553 | 4.72e-06 | 0 |
| ±40% | 0.004999 | 1.103 | 5.81e-06 | 0 |
| ±80% | 0.009726 | 2.193 | 8.13e-07 | 0 |

**Key finding:** perfectly linear response; doubling perturbation doubles both metrics.
`max|dV| < 0.01 pu` even at ±80%. Voltage is stiff; angles carry the response.

This is the **voltage-stiff (DC-approximation-valid) regime**, also described informally
as "angle-dominated." In power systems literature the same condition appears as:
- **"DC power flow validity"**: the DC PF assumption ($|V_i| \approx 1$ pu everywhere,
  only $\theta$ varies) holds when voltage magnitudes barely move under real-power
  perturbations.
- **"small-signal / linear response regime"**: the nonlinear AC PF equations behave
  approximately linearly; doubling the input perturbation doubles the output change.
- **"angle stability regime"** (IEEE/NERC classification): the perturbation stress falls
  in the angle-stability domain rather than the voltage-stability domain.

The voltage-stiff regime is a property of the network's regulated-bus fraction (Section 8):
with 24.5% of buses holding voltage setpoints, reactive redistribution suppresses |V|
deviations and the response is carried almost entirely by bus angles.

### 5.2 Wind curtailment (up to 10%, 20%, 40%, 80%)

All converge. `max|dV| < 0.001 pu` for all levels. Angle shifts grow proportionally
with curtailment (up to ~10 deg at 80% curtailment of 536 MW total wind). Slack bus
absorbs the entire wind reduction.

### 5.3 Generator offline (1, 2, 3, 5, 10 gens)

| Test | MW dropped | max\|dV\| (pu) | max\|dTheta\| (deg) |
|------|-----------|----------------|---------------------|
| bus 147 | 92.4 | 0.002758 | 10.842 |
| 2 rand (seed=99) | 70.4 | — | — |
| 3 rand | 72.0 | 0.000557 | 4.387 |
| 5 rand | 165.3 | 0.002703 | 7.315 |
| 10 rand | 306.3 | 0.004528 | 16.482 |

`max|dV| < 0.005 pu` even with 306 MW dropped (21% of total generation).

---

## 6. PM.jl cross-validation: voltage bias is constant across the mild-perturbation envelope

14 perturbed cases solved with PM.jl and compared to GridKit's solved files:

| metric | min across 14 cases | max across 14 cases | base case |
|---|---|---|---|
| max \|dV\| PM.jl vs GridKit (pu) | 0.027721 | 0.029884 | 0.029728 |
| max \|dTheta\| PM.jl vs GridKit (deg) | 0.194444 | 0.234170 | 0.190703 |
| Ipopt iterations | 4 | 4 | 4 |

**Within the mild-perturbation envelope, the 0.030 pu offset is constant.** It does not
grow with stress *for these 14 cases* (±80% random per-bus load, ±80% wind curtailment,
N≤10 random gens off). Since it cancels in `ΔV = V_perturbed - V_base`, GridKit's relative
voltage response matches PM.jl's inside the envelope. The voltage-stiff regime finding
(Section 5) is confirmed as a network property, not a GridKit artifact.

**⚠ Envelope boundary matters.** Section 11 revisits this comparison under uniform load
scaling (0.1x-2.0x) and largest-gen outages (N=1-25). The offset stays near 0.030 pu only
in the narrow band 0.5x-1.5x uniform load with N≤2 largest gens off. Outside that band the
offset grows sharply (0.089 pu at 2.0x, 0.279 pu at N=25) and violation counts diverge
qualitatively. The relative-ΔV verdict below applies only *inside the envelope*.

**Split verdict for UQ (within the Section 6 envelope):**
- **Angles**: GridKit agrees with PM.jl within 0.23 deg across all cases. Suitable for
  angle-based QoIs (line flows, rotor angle ICs).
- **Absolute voltages**: 0.030 pu offset (3x the 0.01 pu threshold). Use PM.jl for
  voltage-sensitive outputs (Vm constraint checking, Vr/Vi JSON patching).
- **Relative ΔV**: offset is constant *inside the envelope* (variation < 0.003 pu across
  the 14 cases). GridKit is suitable for voltage-sensitivity UQ *within this envelope*.
  Outside the envelope the offset ceases to be constant (Section 11) and this verdict
  does not hold.

---

## 7. Two solution branches (GridKit Section 14, confirmed by PM.jl)

Both solvers show two distinct converged solutions depending on initialization:

| start | V range (pu) | nni (GridKit) / iters (PM.jl) |
|-------|-------------|-------------------------------|
| warm (base Vm/Va) | [1.007, 1.043] | 4–5 / 4 |
| flat (Vm=1, Va=0) | [0.965, 1.000] | 7–8 / 4 |

Separation: ~0.044 pu, constant across base case and all perturbed cases tested. The
high-voltage branch matches the TAMU/PowerWorld reference (the physically operated
solution). The low-voltage branch is a non-operational equilibrium.

This is a structural property of the AC PF equations for ACTIVSg200, not a solver
artifact.

---

## 8. Why voltage is stiff: PV bus regulation (novice explanation)

### Bus types in MATPOWER

| Type | Code | Fixed | Solved | Physical meaning |
|------|:----:|-------|--------|------------------|
| Slack | 3 | \|V\|, θ | — | Reference generator; absorbs P/Q imbalance |
| PV | 2 | \|V\|, P | θ, Q_g | Generator with AVR; voltage held at setpoint |
| PQ | 1 | P, Q | \|V\|, θ | Load bus; voltage determined by network |

At a PV bus, |V| is a boundary condition: the solver never changes it (unless Q limits
are hit). The **regulated fraction** (PV + slack buses / total) determines how much of
the network is voltage-pinned.

### ACTIVSg200: 24.5% regulated

| Case | Buses | PV | Slack | Regulated fraction |
|---|---:|---:|---:|---:|
| ACTIVSg200 (Illinois) | 200 | 48 | 1 | **24.5%** |
| Hawaii40 | 37 | 8 | 1 | **24.3%** |

This is typical for a bulk transmission network (real systems: 15–30%). With ~24% of
buses voltage-pinned, real power perturbations propagate as angle shifts, with reactive
redistribution suppressing voltage deviations. The "angle carries the response, magnitude
stays put" behavior is exactly the assumption behind DC power flow:

$$
P_{ij} \approx \frac{\theta_i - \theta_j}{X_{ij}}
$$

### PV → PQ switching: how Q limits break the stiff-voltage regime

All industry solvers (MATPOWER, PM.jl, PSS/E, PowerWorld) enforce generator Q limits
via **PV → PQ switching** during the solve:

1. **Newton step**: solve for free variables (θ at PV buses, both |V| and θ at PQ buses).
2. **Q-limit check**: for each PV bus, compute the reactive power needed to hold the
   setpoint (`Q_needed`).
   - `Q_min ≤ Q_needed ≤ Q_max`: keep as PV (generator can hold voltage).
   - `Q_needed > Q_max`: exciter saturated. Switch to PQ, fix `Qg = Qmax`, let |V| float.
   - `Q_needed < Q_min`: under-excited limit. Switch to PQ, fix `Qg = Qmin`.

**Physically**: the AVR (automatic voltage regulator) on a synchronous generator
continuously adjusts field current to hold terminal voltage. When the rotor winding
hits its thermal limit (Qmax), the AVR saturates at full field current and terminal
voltage slips below the setpoint. The generator loses voltage control. If many
generators saturate simultaneously (e.g. during a heat wave with high system-wide
demand), cascading voltage sag leads to **voltage collapse**.

**GridKit does not implement this** (verified from source, Section 1). All PV buses stay
PV unconditionally, making the network appear stiffer than it physically is. The 0.030 pu
baseline offset observed at the base case indicates some generators are Q-limited in
PM.jl's base-case solution; identifying the specific generators is a nice-to-have audit
(Section 17, open item #3) but not needed to establish the mechanism. The fact that the
offset is constant across all 14 perturbed cases (Section 6) shows the set of Q-limited
generators does not change under those moderate perturbations; the extra reactive demand
or relief does not push additional generators over their limits. At heavier loading
(Section 9) and large gen outages (Section 10), more generators hit Q limits, PV→PQ
switching activates for additional buses in PM.jl, and voltages sag progressively.
The GridKit-vs-PM.jl gap in that regime is quantified directly in Section 11: at 2.0x
uniform load max|dV| = 0.089 pu (3x the base offset) and at N=10 largest gens off
max|dV| = 0.177 pu (6x the base offset), with the two solvers reporting qualitatively
different violation counts.

---

## 9. PV-curve nose and voltage collapse (PM.jl Section 8a)

### What is the PV-curve nose? (novice explanation)

The **P-V curve** plots bus voltage vs total load. As load increases, a new equilibrium
exists at each step until the **nose** (maximum loadability point):

```
|V|
1.0|****
   |    ***          <- upper branch: V slowly drops as load increases (physical)
   |       **
   |         * <- NOSE (max loadability)
   |       **
   |    ***          <- lower branch: same P, much lower V (unstable, not physical)
 0 |****
   +-------------------
   0                 P_max   P (total load)
```

- **Upper branch**: physically operated; V decreases slowly as load increases. Real
  grids operate here.
- **Lower branch**: mathematically valid solution at the same P values but much lower V.
  Physically unstable; no real grid operates here.
- **Nose**: the two branches meet at maximum loadability P_max. The PF Jacobian is
  singular at this point; no solution exists beyond it. = voltage collapse.

The distance from operating point to nose (in MW) is the **voltage stability margin**.

### ACTIVSg200 load scaling results

Uniform load scaling (all Pd/Qd × `scale`), PM.jl warm-start:

| scale | total Pd (MW) | converged | V_min (pu) | violations | max\|dV\| vs base (pu) | max\|dTheta\| vs base (deg) |
|------:|-------------:|:---------:|----------:|----------:|----------------------:|--------------------------:|
| 0.1x | 147.6 | yes | 1.036 | 85 (over-V) | 0.038 | 28.1 |
| 0.5x | 737.8 | yes | 1.033 | 44 (over-V) | 0.023 | 15.8 |
| 1.0x | 1475.7 | yes | 1.010 | 1 | 0 | 0 |
| 1.5x | 2213.5 | yes | 0.980 | 0 | 0.035 | 17.3 |
| 2.0x | 2951.4 | yes | 0.926 | 31 (under-V) | 0.096 | 38.4 |
| 3.0x+ | 4427+ | **no** | — | — | — | — |

**Nose location: between 2x and 3x base load (2951–4427 MW).** Both warm and flat start
fail at 3x (Ipopt reports infeasibility, not convergence failure). Voltage stability
margin ≥ 100% additional load.

### Voltage violations as a leading indicator

- **0.1x–0.5x**: over-voltage (generators over-excite with low reactive demand)
- **1.0x**: 1 marginal violation
- **1.5x**: **zero violations** — extra reactive demand pulls voltages into band
- **2.0x**: 31 buses below 0.95 pu — PV→PQ switching activated, collapse precursors
- **3.0x**: no solution (past the nose)

### Where the voltage-stiff regime ends

The voltage-stiff behavior from Section 5 holds only for moderate perturbations:

| scale | max\|dV\| (pu) | max\|dTheta\| (deg) | regime |
|------:|---------------:|--------------------|--------|
| 1.0x (base) | 0 | 0 | — |
| ±80% random per-bus | 0.010 | 2.2 | voltage-stiff / DC-valid (Section 5) |
| 1.5x uniform | 0.035 | 17.3 | transitional |
| 2.0x uniform | 0.096 | 38.4 | nonlinear (V + θ both large) |

At 1.5x uniform loading, V_min drops from 1.010 to 0.980 pu (a 0.030 pu drop); at 2.0x,
V_min = 0.926 pu. The PV buses that were unconditionally regulating voltage in Section 5
have begun hitting Qmax and switching to PQ, allowing bus voltages to sag. The DC PF
approximation becomes inaccurate.

**The voltage-stiff regime is bounded at roughly ±80% random per-bus load** (equivalent to
~0.2x–1.8x total system load). Beyond that, uniform load increases stress all generators
simultaneously, triggering Q limits earlier than random perturbations (which partially
cancel across buses).

---

## 10. Generator outage stress test (PM.jl Section 8b)

The 1, 2, 5, 10, 15, 20, 25, 30 largest non-slack generators (by base-case Pg) were
taken offline in sequence. `case_ACTIVSg200` has 49 generators total; base total Pd =
1475.7 MW. Pg of all remaining non-slack generators is fixed at base-case values; the
slack bus absorbs the entire real power shortfall with no upper limit in the PF
formulation.

| N off | MW dropped | warm conv | warm V_min (pu) | violations | dV vs base (pu) | dTheta vs base (deg) | flat conv | warm vs flat (pu) |
|------:|----------:|:---------:|----------------:|-----------:|----------------:|---------------------:|:---------:|------------------:|
| 1  | 154.8  | yes | 1.010 | 1   | 0.005 | 7.6  | yes | 0.045 |
| 2  | 288.7  | yes | 1.006 | 1   | 0.011 | 9.6  | yes | 0.045 |
| 5  | 648.5  | yes | 0.989 | 0   | 0.039 | 21.0 | yes | 0.048 |
| 10 | 919.7  | yes | 0.841 | **85**  | 0.177 | 29.0 | yes | 0.084 |
| 15 | 1020.8 | yes | 0.821 | **96**  | 0.197 | 31.3 | yes | 0.098 |
| 20 | 1062.7 | yes | 0.809 | **99**  | 0.208 | 32.4 | yes | 0.110 |
| 25 | 1086.9 | yes | 0.736 | **177** | 0.279 | 36.3 | **no** | — |
| 30 | 1095.6 | **no** | — | — | — | — | **no** | — |

### Why the PF converges with most generation offline

The slack bus has no Pg limit in an AC PF formulation. The solver finds (V, θ)
satisfying Kirchhoff at every bus; whatever real power the slack needs to supply
emerges as a consequence, not a constraint. So removing 25 generators (1087 MW) is
mathematically solvable because the slack can cover the shortfall. In a real dispatch
sense this is infeasible, but the PF equations do not represent dispatch.

**"Converges" does not mean "acceptable."** The voltage column is the honest indicator.
At N=10 (920 MW dropped), V_min = 0.841 pu with 85 buses below 0.95 pu: a voltage
emergency by any standard. At N=25, V_min = 0.736 pu with 177 violations. The system
is in deep voltage collapse territory even though a mathematical PF solution still exists.

### Asymmetric branch disappearance

In the load-scaling experiment (Section 9) both branches (warm and flat) fail
simultaneously at 3x load. Here they fail asymmetrically: flat-start (low-voltage
branch) fails at N=25 while warm-start (high-voltage branch) survives to N=30.
Generator outages remove reactive support (PV buses go offline), which stresses the
two equilibria differently from uniform load increases. The low-voltage branch is less
robust to reactive support loss and merges with the upper branch first.

### Two-branch gap grows with outage count

The warm-vs-flat gap is constant at 0.044 pu under moderate perturbations (Section 9
base case, Section 6 all 14 perturbed cases). Here it grows monotonically:
0.045 (N=1) → 0.048 (N=5) → 0.084 (N=10) → 0.110 (N=20). Each generator outage
removes a voltage-supporting PV bus, shrinking the distance between the two equilibria.

### Voltage-stiff regime breaks earlier per MW than load scaling

| stress | max dV (pu) | max dTheta (deg) |
|--------|------------|------------------|
| ±80% random per-bus load | 0.010 | 2.2 |
| N=2 gen off (289 MW) | 0.011 | 9.6 |
| N=5 gen off (649 MW) | 0.039 | 21.0 |
| 1.5x uniform load (Section 9) | 0.035 | 17.3 |

Gen outages exit the voltage-stiff regime much faster per MW than uniform load scaling.
A 289 MW outage already produces dTheta = 9.6 deg (vs ~3 deg for 289 MW of extra
random load). The slack concentrates reactive deficit at one location rather than
distributing it across all buses.

---

## 11. GridKit vs PM.jl under stress (Section 9 in `pm_helper.ipynb`)

Section 6 established that GridKit's ~0.030 pu voltage bias vs PM.jl is constant across
the 14 mild perturbations (±80% random per-bus load, ±80% wind curtailment, N≤10 random
gens off). Two questions were left open: (a) does the constant offset hold under uniform
stress? (b) what happens at large gen outages? The stress-test comparison in
`pm_helper.ipynb` Section 9 answers both.

Method: both solvers wrote solved `.m` files to `pm-solver/m-cases/stress-test/`
(GridKit via `run_solve_pf_out` in `pf_helper.ipynb` Section 8; PM.jl via
`run_pm_solve_out` in `pm_helper.ipynb` Section 8). Section 9 reads both sides and
computes bus-by-bus `|dV|`, `|dTheta|`. No re-solving.

### 11.1 Load-scale stress (warm-start, both solvers)

| scale | PM V_min | GK V_min | PM viols | GK viols | max\|dV\| (pu) | max\|dTheta\| (deg) |
|-----:|--------:|--------:|--------:|--------:|--------------:|--------------------:|
| 0.1x | 1.036 | 1.029 | **85** | 0 | 0.0373 | 0.36 |
| 0.5x | 1.033 | 1.029 | **44** | 0 | 0.0355 | 0.21 |
| 1.0x | 1.010 | 1.009 | 1 | 0 | 0.0297 | 0.19 |
| 1.5x | 0.980 | 0.981 | 0 | 0 | 0.0305 | 0.33 |
| 2.0x | 0.926 | 0.929 | **31** | **21** | 0.0890 | 0.45 |
| 3.0x-8.0x | PM.jl did not converge (past PV-curve nose, Section 9) |

The 0.030 pu constant offset holds only in a narrow band (0.5x-1.5x). Outside that
band, on both sides:
- **Under-loaded (0.1x, 0.5x)**: PM.jl reports 44-85 **over-voltage** violations
  (V > 1.05 pu); GridKit reports zero. Under light load, physical generators must
  absorb reactive power to keep buses at setpoint; when they hit `Qmin` (under-excited
  limit) they switch to PQ and voltage climbs. GridKit's non-enforcement of `Qmin`
  keeps its generators absorbing unlimited reactive, so voltages stay clamped at the
  setpoint.
- **2.0x**: max|dV| jumps from 0.030 to 0.089 pu (3x baseline), max|dTheta| from
  0.19 to 0.45 deg. GridKit still gets the qualitative picture (V_min agrees to
  0.003 pu, 21 vs 31 violations), but the constant-offset assumption has broken.

### 11.2 Gen-outage stress (warm-start, both solvers): sharp transition at N=10

| N off | PM V_min | GK V_min | PM viols | GK viols | max\|dV\| (pu) | max\|dTheta\| (deg) |
|-----:|--------:|--------:|--------:|--------:|--------------:|--------------------:|
| 1  | 1.010 | 1.009 | 1 | 0 | 0.0294 | 0.20 |
| 2  | 1.006 | 1.007 | 1 | 0 | 0.0288 | 0.20 |
| 5  | 0.989 | 1.002 | 0 | 0 | 0.0390 | 0.17 |
| 10 | **0.841** | **0.996** | **85** | **0** | **0.177** | **3.85** |
| 15 | 0.821 | 0.994 | 96 | 0 | 0.197 | 4.86 |
| 20 | 0.809 | 0.993 | 99 | 0 | 0.208 | 5.50 |
| 25 | 0.736 | 0.993 | **177** | **0** | **0.279** | 9.72 |
| 30 | PM.jl did not converge (network has no equilibrium, Section 10) |

The transition between N=5 and N=10 is the **operational boundary of GridKit's
usability on ACTIVSg200** for gen-outage stress:
- N=1, 2: baseline 0.029 pu offset (matches Section 6).
- N=5: 0.039 pu, first mild expansion.
- **N=10: 0.177 pu (6x baseline), max|dTheta| jumps 0.17 → 3.85 deg (22x jump).**
- N=25: 0.279 pu, max|dTheta| = 9.72 deg. GridKit reports V_min = 0.993 pu with **zero
  violations** while PM.jl reports V_min = 0.736 pu with **177 violations**. The two
  solvers report qualitatively different network states.

Mechanism: dropping generators strands reactive support at neighboring PV buses; those
buses hit `Qmax` in PM.jl and switch to PQ (voltages sag); in GridKit they keep holding
their setpoints, so the network appears artificially stiff. Same PV→PQ non-enforcement
mechanism as under-load, opposite Q-limit sign.

### 11.3 Behavioral confirmation of the Section 1 source-inspection finding

The two stress sweeps together provide behavioral evidence that agrees exactly with the
source-inspection result of Section 1 (`grep` found zero uses of `Qmax`/`Qmin` in
GridKit's PF solver path):

| stress direction | PM.jl behavior | GK behavior | Q-limit missed |
|------------------|---------------|-------------|----------------|
| Light load (0.1x-0.5x) | 44-85 buses > 1.05 pu | 0 violations | Qmin (absorbing limit) |
| Base case (1.0x) | 1 violation, V_min = 1.010 | 0 violations, V_min = 1.009 | 0.030 pu offset persists (a few gens at limit) |
| Heavy load (2.0x) | 31 buses < 0.95 pu | 21 buses < 0.95 pu | Qmax (producing limit), partial |
| N=10-25 gens off | 85-177 violations | 0 violations | Qmax, extreme |

Both Q-limit signs are missed; the offset is symmetric. This is exactly the fingerprint
of no PV→PQ switching, independently confirming the Section 1 source-inspection
finding.

### 11.4 Threshold-based verdict

| threshold | 11.1 (load-scale) violated at | 11.2 (gen-outage) violated at |
|-----------|------------------------------|------------------------------|
| max\|dV\| < 0.010 pu | never (baseline offset already fails) | never (same reason) |
| max\|dV\| < 0.030 pu (base-case bias, ΔV-cancellation limit) | breaks at 2.0x (0.089 pu); qualitatively at 0.1x-0.5x (viols mismatch) | breaks at N=10 (0.177 pu) |
| max\|dTheta\| < 1.0 deg (Section 6 threshold) | never in the converged range | breaks at N=10 (3.85 deg) |

The `DV_ABS_THRESH_PU` of 0.010 pu was already violated by the base offset, so
absolute-voltage QoIs require PM.jl regardless (as established in Section 6). What is
new: the `ΔV_offset_variation` threshold of 0.005 pu from Section 7 is exceeded at
2.0x uniform load and at N=10 gens off. GridKit's relative-ΔV verdict (Section 7)
holds only inside the Section 6 envelope; it does **not** extend to uniform stress or
large gen outages.

### 11.5 Boundary of GridKit usability (revised)

The Section 7 verdict ("angles PASS, ΔV PASS within envelope") is confirmed only for
the 14 mild perturbations (±80% random per-bus load, ±80% wind, N≤10 random gens off).
It does **not** extend to:
- Uniform load < 0.5x or > 1.5x (over/under-voltage violations mis-counted).
- Uniform load ≥ 2.0x (max|dV| grows past 0.030 pu, ~3x baseline).
- N ≥ 10 largest gens off (max|dV| and max|dTheta| both explode; violation counts
  qualitatively wrong).

---

## 12. Implication for UQ

The aleatoric UQ perturbations (±5–20% per-bus random load) sit deep in the
voltage-stiff (DC-approximation-valid, angle-dominated) regime:
- Far from the PV-curve nose (margin > 100% additional load)
- Far from Q-limit activation (no PV→PQ switching at these perturbation levels)
- Linear angle response (confirmed: doubling perturbation doubles metrics)

### Scope: this applies only to PF initial conditions

The findings above describe the **steady-state PF solution** that becomes the initial
condition (t=0) for the downstream phasor-dynamic simulation. They say that across
perturbed operating points, the `.m` file's Vm/Va values change very little (voltage)
and modestly (angle) compared to the base case.

The aleatoric UQ study will apply a **fault** at t > 0 and monitor dynamic response
(rotor angles, rotor speed, bus voltage magnitude and angle, frequency nadir, damping)
over ~10 s of simulated time. A steady-state voltage-stiff PF **does not** guarantee that
`V_mag(t)` will be stiff during the fault transient. Dynamic V_mag response depends on
generator inertia (H), AVR/exciter dynamics, fault location and clearing time, and
whether AVRs saturate during the transient. Analogous statements hold for rotor angle
swings, which depend on synchronizing torque coefficients, not just steady-state
angle differences.

### Implications for the dynamic UQ pipeline

- **Do not drop V_mag(t) as a dynamic-simulation QoI a priori.** That decision must be
  made from dynamic-simulation output, not from these PF-only findings.
- **PF-side QoIs** (base-case Vm, base-case Va, base-case branch MW flow) *can* be
  characterized directly from the PF sweep: their variation across 8760 hourly scenarios
  is bounded by the response linearity shown in Sections 5–6 and 9.
- **Dynamic QoIs** (max |dTheta| during swing, frequency nadir, first-swing V_mag dip,
  settling time) still need to be measured on actual phasor-dynamic runs. Surrogate
  modeling of these QoIs (polynomial chaos expansion (PCE), Gaussian processes (GPs))
  is a downstream step after enough dynamic runs are collected.
- The 8760-scenario aleatoric run is the natural next step: solve PF for each hourly
  operating point, patch into `illinois.json`, run the phasor-dynamic simulation with
  a fixed fault scenario, collect dynamic QoIs. See `uq_plan.md` Task 0b for the
  prototype-then-scale plan.

### Solver choice for PF stage

- PM.jl for any PF output requiring absolute voltage accuracy (`Vr`/`Vi` JSON patching
  into `illinois.json`, absolute voltage constraint checking).
- GridKit for rapid relative comparisons where only angle/flow changes matter.

### When the voltage-stiff PF assumption itself breaks

If the study is extended to correlated large-scale load increases (e.g. heat wave
+30–50% system load), or to hourly PCM cases with dispatch far from the base case,
the voltage-stiff assumption breaks down at the PF stage (Section 9 shows the boundary
is around 1.5x uniform load). At that point PM.jl must be used, and the linearity of
PF-side QoIs vs perturbation is no longer valid.

---

## 13. GridKit `solve_pf` implementation details

### Local parser (`matpower_parser.hpp`)

Replaces GridKit's `MatpowerParser.hpp` for this use case:

| Feature | GridKit parser | Local parser |
|---------|---------------|--------------|
| Unknown field names (e.g. `bus_name`) | regex fail | skip silently |
| Extra columns in gen/branch | throws | ignores extras |
| Cell arrays `{...}` | not handled | consumed and skipped |
| `BusData.Va` unit | degrees (raw) | converted to radians |
| `Pd`, `Qd`, `Pg`, `Qg` unit | MW/MVAr (raw) | divided by baseMVA → pu |

### `--output-m`: solved .m output

Rewrites the input `.m` with Vm/Va columns (8-9 of `mpc.bus`) updated to the PF
solution. All other content passed through unchanged.

Pipeline:
```
modified .m → solve_pf input.m --output-m solved.m → m_to_case.py solved.m → case.json
```

### Build

```bash
bash ~/gridkit/uq-usecase/pf-solver/build.sh
# output: uq-usecase/pf-solver/build/solve_pf
```

---

## 14. Transformer tap ratio: what GridKit omits (reference)

### What the MATPOWER model includes

For a transformer with turns ratio `a` and phase shift `φ`:

$$
Y_{ij} = -\frac{y_s}{a \cdot e^{j\phi}}, \quad
Y_{ii} += \frac{y_s}{a^2}, \quad
Y_{jj} += y_s
$$

When `a=1, φ=0` this reduces to the plain PI model that GridKit implements.

### ACTIVSg200: no off-nominal transformers

All 245 branches have `TAP=0` (179 plain lines) or `TAP=1.0` (66 explicit unity),
`SHIFT=0` for all. The tap-ratio omission has **no effect** on this case. Tap modeling
is not the source of the 0.030 pu bias; that is caused by Q-limit non-enforcement
(verified from source, Section 1).

### GridKit code path

`BranchData` parses `ratio` and `angle` from columns 9-10, but `Branch::Branch()`
only reads `data.r`, `data.x`, `data.b`. The transformer fields are silently dropped.
A fix would be needed for cases with LTC transformers (`ratio != 1`).

---

## 15. Stagnation investigation (GridKit Section 14)

**Concern**: warm-start perturbed solutions might be KINSOL stagnating at the initial
point rather than finding a genuine new equilibrium.

**Resolution**: `nni=4-5` for all warm-start perturbed cases (stagnation would show
`nni=1`). Flat-start consistently converges to the low-voltage branch (V~[0.965, 1.000])
with `nni=7-8`, confirming different Newton trajectories. Warm-start solutions are
genuine equilibria of GridKit's model.

---

## 16. Recommendation: use PM.jl as the production PF solver

### Why PM.jl over GridKit `solve_pf` for the aleatoric UQ pipeline

GridKit `solve_pf` is functional and useful for angle-based comparisons, but has two
verified structural gaps for the ACTIVSg200 UQ study:

1. **No PV→PQ switching** (Section 1, verified across six independent code paths;
   Section 11, confirmed behaviorally). Produces a ~0.030 pu voltage bias in the
   voltage-stiff regime that is *constant only* inside the Section 6 envelope
   (0.5x-1.5x uniform load, N≤2 largest gens off). Outside that envelope the
   divergence is measured (Section 11): max|dV| = 0.089 pu at 2.0x uniform load
   (3x baseline), 0.177 pu at N=10 largest gens off (6x baseline), 0.279 pu at N=25
   (with GridKit reporting 0 violations vs PM.jl's 177). Both Q-limit signs are
   missed: over-voltage violations at 0.1x-0.5x load (Qmin hits) and under-voltage
   violations at large gen outages (Qmax hits).
2. **No transformer tap ratio / phase shift** (Section 14). Irrelevant for ACTIVSg200
   (all TAP=1, SHIFT=0), but blocks any future extension to networks with LTC
   transformers (which the ACTIVSg2000 and most real industry cases contain).

PM.jl (via `pm_solve.jl`, Section 3) is a MATPOWER-equivalent solver with correct
PV→PQ switching, correct tap/shift handling, and cross-validation against the TAMU
reference to 1e-5 pu (Section 4). Its per-call overhead is ~8.5 s (5.5 s Julia startup +
3 s Ipopt); for the 8760-scenario aleatoric run, a batched-manifest mode (one Julia
process loops over all cases) or Slurm sbatch parallelization keeps total wall time
tractable.

### Future development plan: how to bring GridKit `solve_pf` to PM.jl parity

If a decision is made later to invest in GridKit rather than lock in PM.jl, the
following sequence of work items is required. Each item is scoped to a concrete
surface in the current source tree (Section 1 identifies file locations).

**Phase 1: PV→PQ switching (mandatory; unblocks Section 11 correctness).**
This is the single largest gap. Required pieces:

1. **Store Q at PV buses.** Uncomment the `Q_(data.Qg)` member and the
   `bus_->Q() += Q_` residual contribution in `Model/PowerFlow/Generator/GeneratorPV.cpp`
   (currently both are commented out; a PV generator injects only P into the residual
   today, so the reactive residual is structurally absent). This requires deciding
   whether `Q_g` is a solver unknown (as it is in MATPOWER's PV bus formulation) or
   a post-processed diagnostic.
2. **Expose Qmax/Qmin to the solver.** They already live in `GenData` (Section 1);
   plumb them into `GeneratorPV` at construction and into a per-Newton-step check
   inside `SystemModelPowerFlow::evaluateResidual()`.
3. **Add a bus-type flip mechanism.** `BusFactory` currently dispatches on `data.type`
   once at construction and returns a raw pointer to a concrete `BusPQ`/`BusPV`/`BusSlack`
   subclass (Section 1, `Bus/BusFactory.hpp`). To flip a bus mid-solve, either:
   (a) merge the three subclasses into a single `Bus` class carrying a runtime `type_`
   field and dispatching residual behavior inside `evaluateResidual()`, or
   (b) keep the polymorphic hierarchy but tear down/rebuild the object when a flip
   is triggered, then rewire the `Generator` pointer to the new object. (a) is a
   cleaner long-term choice but more invasive.
4. **Wrap KINSOL in an outer switching loop.**
   `Solver/SteadyState/Kinsol.cpp::runSimulation()` calls `KINSol(...)` exactly once
   today. Add an outer loop in `SystemModelPowerFlow` (or a new `KinsolPowerFlow`
   wrapper) that: (i) calls `runSimulation()`, (ii) checks Q at every PV bus against
   its limits, (iii) flips the bus type and re-initializes the residual for any bus
   past its limit, (iv) re-runs KINSOL, (v) terminates when no bus flips in an iteration
   or a max-outer-iter cap is hit. Match MATPOWER's default of ≤ 10 outer iterations.
5. **Cross-validate against MATPOWER.** Ship a small `pytest`-style suite that runs
   the base case + all 17 Section 11 stress cases through both GridKit and PM.jl
   (or a MATPOWER reference), and asserts per-bus max|dV| < 1e-4 pu and
   max|dTheta| < 1e-3 deg. This is what "parity" means operationally.

**Phase 2: Transformer tap ratio and phase shift (needed for ACTIVSg2000 and real cases).**
Currently `Branch::Branch()` in `Model/PowerFlow/Branch/Branch.cpp` reads only
`data.r`, `data.x`, `data.b`; the `ratio` and `angle` fields are silently dropped
(Section 14). Extend `Branch` with the standard MATPOWER off-nominal transformer
formulation ($Y_{ij} = -y_s/(a e^{j\phi})$, $Y_{ii} \mathrel{+}= y_s/a^2$,
$Y_{jj} \mathrel{+}= y_s$). Small change; no solver-loop implications.

**Phase 3: Bus shunts (Gs, Bs) — smaller, but touches nodal reactive balance.**
ACTIVSg200 has ~24 buses with nonzero Bs. The contribution has not been audited
(Section 17, open item #2). Add `data.Gs`, `data.Bs` reads and inject
$P_{shunt} = |V|^2 G_s$, $Q_{shunt} = -|V|^2 B_s$ into the bus residual.

**Phase 4: Optional MATPOWER features that PM.jl already supports.**
DC lines, HVDC, area/zone dispatch coupling, unit commitment linkage. Deferred
unless the study scope changes.

**Effort estimate.** Phase 1 is the significant one (multiple weeks of dev + validation).
Phases 2-3 are days each. Phase 4 is open-ended. Even after Phase 1 completes, GridKit
would only reach MATPOWER-equivalent PF accuracy on the network types PM.jl already
handles for free.

**Why PM.jl is still the current recommendation.** All of Phases 1-3 are already
solved and actively maintained inside PowerModels.jl. Investing GridKit dev-time only
makes sense if the downstream project needs tight coupling between PF and GridKit's
transient/dynamic solvers (`SystemModel`, `ModelEvaluatorImpl`) that would be awkward
to express through a subprocess call to a Julia binary. If the PF stage stays
decoupled from the dynamic stage (which is the current architecture: PF output is
patched into `illinois.json`), PM.jl is the better default indefinitely.

### Recommended split for production

| Task | Solver |
|------|--------|
| PF solves feeding `Vr`/`Vi` into `illinois.json` (all 8760 scenarios) | **PM.jl** |
| Base-case constraint checking, voltage stability screening | **PM.jl** |
| Rapid what-if angle/flow comparisons during development | GridKit (still useful for relative deltas within the voltage-stiff regime) |
| Any case with tap ratios != 1, phase shifters, LTC transformers | **PM.jl** |
| Any case near or past the PV-curve nose (e.g. contingency screening past 1.5x load) | **PM.jl** |

Both solvers stay in the repo. PM.jl becomes the primary; GridKit becomes a fast
second opinion for angle-only checks.

---

## 17. Remaining open items

1. **Regulated-fraction check on ACTIVSg2000**: verify similar regulated fraction for
   the 2000-bus Texas case.
2. **Bus shunts (Gs, Bs)**: GridKit's `Branch` hardcodes shunt conductance to zero.
   ACTIVSg200 has ~24 buses with nonzero Bs (capacitor banks). Their contribution to
   the nodal reactive balance has not been audited.
3. **Reactive dispatch audit** (low priority, nice-to-have): tabulate `Qg vs
   [Qmin, Qmax]` from PM.jl base-case solution to identify *which* generators are
   Q-limited. Not needed to establish the mechanism (verified from GridKit source in
   Section 1 and behaviorally in Section 11), only informational.

**Recently closed:**
- **GridKit vs PM.jl at 1.5x+ uniform load and large gen outages** (was highest
  priority): closed by Section 11. GridKit's error grows sharply outside the Section 6
  envelope; boundary quantified for both stress directions.