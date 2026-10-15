# m_to_case helper: MATPOWER `.m` → GridKit `case.json`

Companion doc for [`notebooks/m_to_case_helper.ipynb`](../notebooks/m_to_case_helper.ipynb).

Purpose: document and rigorously test the `.m` → `case.json` mapping, so the
aleatoric UQ workflow can rest on it.

## Scope

- **In scope**: the production `.m` + `.dyr` → `case.json` contract; the
  one-time `.aux` bootstrap that validates it; scenario perturbations (load,
  wind, gen-off); tolerance philosophy.
- **Out of scope**: PF solving itself (see [`pm_helper.md`](pm_helper.md) /
  PM.jl), running `case.json` through `DynamicSimulation` (`aleatoric_helper.md`,
  to come), unit-commitment decisions.

## Two phases: bootstrap once, then produce forever

The distinction that matters most in this document. PowerWorld appears **only**
in the bootstrap, and never at runtime.

```
BOOTSTRAP  (once per case, to establish trust)

   ACTIVSg200.aux + ACTIVSg200_dynamics.dyr
        │
        ├─ build_case_from_aux_dyr()  ─→ rebuild of illinois.case.json
        │       └─ diff vs the shipped file → 0 diffs (with creator overlay)
        │          proves the whole field mapping is right
        │
        └─ write_matpower_m()         ─→ case_ACTIVSg200_pw.m
                └─ a seed .m whose operating point matches the shipped case
   ── PowerWorld exits here, permanently ──


PRODUCTION  (every scenario)

   seed .m ──(turn gens off: thin inertia)──→ base .m          [once]
        │
        └─ hourly load / wind  ──(PM.jl)──→ solved .m          [×8760]
                │
                └─ build_case_from_m(base_case, solved_m, dyr)
                        └─→ case.json ──→ DynamicSimulation
```

Both production inputs are open standards: MATPOWER `.m` (Sienna, PowerModels.jl
and every PCM tool speak it) and PSS/E `.dyr`. The third input is the base
GridKit `case.json`, which supplies topology, device structure, the ZIP mix and
the fault devices — none of which a PCM run produces.

Because generator commitment is fixed during bootstrap, the 8760-scenario loop
never adds or removes devices. It only patches bus `init.Vr/Vi`, Genrou
`p0/q0`, and `LoadZIP.{Pnom,Qnom}` — the narrowest and best-tested code path.

## Entry points

| Function | Phase | Use when |
|---|---|---|
| `build_case_from_aux_dyr()` | bootstrap | reproduce a shipped case from TAMU source data; generate the creator overlay |
| `write_matpower_m()` | bootstrap | emit a MATPOWER `.m` from `.aux`-derived frames; no template `.m` is read |
| `build_case_from_m()` | production | map a solved `.m` (+ `.dyr`) onto the base case |
| `patch_case_from_m()` | legacy | epistemic track (`uq_setup.ipynb`); bus init + Genrou `p0/q0` only. See the import bug in [Known limitations](#known-limitations) |


## How the pipeline was validated

Nothing here rests on "the mapping looks right". Every stage was checked against
an independent oracle, and the checks are cheap enough to rerun as regressions.

### The `.m` builder is derived only from the `.aux`

`write_matpower_m()` reads **no template `.m`**. The MATPOWER table layout is a
documented constant in [`py-utils/aux_to_m.py`](../py-utils/aux_to_m.py)
(`BUS_COLS`, `GEN_COLS`, `BRANCH_COLS`, `GENCOST_COLS`), and every value comes
from the `.aux`.

The original `case_ACTIVSg200.m` is deliberately **not** used, for three reasons
that were measured rather than assumed:

1. It is a different operating point (see
   [Where `illinois.case.json` actually comes from](#where-illinoiscasejson-actually-comes-from-verified)).
2. It is missing branch **(64, 82)**, which exists in the `.aux` and in the
   shipped case as `BR_82_64_1` (R=0.023, X=0.158, B=0.026). The `.m` has 245
   branches; the `.aux` and `case.json` have 246.
3. It is a rounded export. Over the 245 shared branches,
   **max |dR| = 5.0e-07, max |dX| = 5.0e-07, 245/245 within 1e-6** — the `.aux`
   carries ~15 significant figures where the `.m` stores ~7.

Two `.aux` details that are easy to get wrong:

- The `Branch` object appears in **two `DATA` blocks with different schemas**:
  180 **Lines** carry impedance in `LineR/LineX/LineC`, while 66
  **Transformers** use `LineR:1/LineX:1/LineC:1`. Reading `LineR` for a
  transformer silently yields `NaN`.
- All 66 transformers are `TAP = 1.0`, `SHIFT = 0.0` (`LineTap:1`, `LinePhase`).
  ACTIVSg200 has no off-nominal taps and no phase shifters.

`gencost` is fully derivable — PowerWorld stores a cubic input/output curve
scaled by fuel price:

```
c0 = GenFixedCost             (49/49 exact vs the original .m)
c1 = GenIOB * GenFuelCost     (8.796296 * 2.16 = 19.0)
c2 = GenIOC * GenFuelCost     (0.000926 * 2.16 = 0.002)
c3 = GenIOD * GenFuelCost     (zero throughout -> NCOST = 3, quadratic)
```

`genfuel` comes from `GenFuelType`, `bus_name` from `BusName`.

**Not present in a `.aux`**, so MATPOWER defaults are used and reported by the
function through `AUX_ABSENT`: gen `PC1/PC2/QC1MIN/QC1MAX/QC2MIN/QC2MAX/RAMP_*/APF`,
branch `RATE_B/RATE_C/ANGMIN/ANGMAX`, gencost `STARTUP/SHUTDOWN`. `baseMVA` is a
parameter (100). Per-bus voltage-limit fields are empty in the ACTIVSg exports,
so `VMAX/VMIN` default to 1.1/0.9 — a `write_matpower_m` argument, not a
hardcode. **This matters if you run OPF with voltage constraints.**

### Oracle 1 — AC power-flow mismatch

Compute `S = V · conj(Y_bus V)` at the stored voltages and compare against the
scheduled injections. If any impedance, tap, shunt or sign convention were
wrong, the mismatch would be large at the affected buses.

| | template-based build | **pure `.aux` build** |
|---|---|---|
| max \|ΔP\| over PQ+PV buses | 0.0078 MW | **0.0063 MW** |
| max \|ΔQ\| over PQ buses | 0.0045 MVAr | **0.0018 MVAr** |
| network losses vs `PG − PD` | 24.749 / 24.75 MW | **24.749 / 24.749 MW** |

The residual is the ~8-decimal write precision, not a modelling error. This is
also what confirms branch (64, 82) belongs: omit it and buses 64 and 82 show
large mismatch.

### Oracle 2 — independent confirmation by PM.jl

```
parsed case_ACTIVSg200_pw.m: 200 buses, 49 gens (9 offline), baseMVA=100
Number of Iterations....: 4
Constraint violation....: 1.6e-13
EXIT: Optimal Solution Found.   termination_status=LOCALLY_SOLVED
```

Converging in 4 iterations at a 1e-13 violation is what an already-solved input
looks like.

### Oracle 3 — round trip back to the shipped case

| Stage | Result |
|---|---|
| `.aux` + `.dyr` → `case.json`, **with** creator overlay | **0 diffs** vs shipped |
| `.aux` + `.dyr` → `case.json`, **without** overlay | only `197_1` differs above 1e-4 (14 params) |
| `case_ACTIVSg200_pw.m` → `case.json` | **0 diffs @ 1e-6**, 1 @ 1e-7, 13 @ 1e-8 |

Below 1e-4 the extra rows are `.dyr` print rounding: `H` on all 40 Genrou,
`Tgov1.T2/T3` on 24 devices, and `76_1.Xdpp`. The `.dyr` and `.dyd` print 4
decimals while the JSON keeps full precision (`49_1` `H`: `.dyr` 5.2115, JSON
5.21145296).

### Operating-point facts that fall out of the build

| Quantity | Value |
|---|---|
| total PD / QD | 2178.00 MW / 620.73 MVAr |
| total PG (online) | 2202.75 MW |
| losses | 24.749 MW |
| gens online / total | 40 / 49 |
| bus types | 151 PQ, 48 PV, 1 slack (bus 189) |
| branches (transformers) | 246 (66) |

**The slack is at its limit**: bus 189 has PG = 569.14 MW against PMAX = 569.1.
Any generator switched off in a plain PF dumps its output entirely on the slack,
so switching off anything larger than a few MW needs an explicit redistribution
rule rather than letting the slack absorb it.

There are **0** `PG > PMAX` violations against the `.aux` limits, which are wider
than the original `.m`'s. Five `PG < PMIN` rows remain, all wind (65, 104, 105,
114, 147) curtailed below economic minimum — fine for PF, relevant only if OPF
is added, where non-dispatchables need `PMIN` relaxed.

## PF solver choice and the solved `.m` contract

The pipeline's PF-solve stage (stage 3) can use either solver. Both read a
MATPOWER `.m` and produce a solved `.m` with updated bus voltages:

| Solver | Wrapper | Updates in solved `.m` | PV→PQ switching |
|---|---|---|---|
| PM.jl | `pm_utils.run_pm_solve_out` | `mpc.bus` VM, VA; **`mpc.gen` PG, QG** | Yes (correct Q limits) |
| GridKit `solve_pf` | `pf_utils.run_solve_pf_out` | `mpc.bus` VM, VA only | No (see below) |

Both return `(bus_df, stderr, return_code)` with columns
`[bus_i, V_pu, theta_deg, type]`.

### What `build_case_from_solved_m` reads from the solved `.m`

| Field | Source column | Used for | PF-solve role |
|---|---|---|---|
| VM, VA | `mpc.bus[:, 8:9]` | bus `init.Vr/Vi` | Solved (all buses). Both solvers update. |
| PD, QD | `mpc.bus[:, 3:4]` | LoadZIP `Pnom/Qnom` | Input (unchanged by PF solve). |
| PG (non-slack) | `mpc.gen[:, 2]` | Genrou `p0` | Input (dispatch from PCM). Unchanged by both solvers. No gap. |
| PG (slack gen) | `mpc.gen[:, 2]` | Genrou `p0` | Solved (absorbs losses + imbalance). PM.jl updates; GridKit does not. |
| QG (slack + PV gens) | `mpc.gen[:, 3]` | Genrou `q0` | Solved (reactive needed to hold voltage). PM.jl updates; GridKit does not. |
| GEN_STATUS | `mpc.gen[:, 8]` | offline-gen removal | Input (unchanged by PF solve). |

**Summary of the gap**: for ACTIVSg200, GridKit's output `.m` has stale values
in 1 slack-gen PG + 49 gen QG (1 slack + 48 PV-bus). All other PG values are
dispatch inputs and are correct in both solvers.

### The PG/QG gap in GridKit's `solve_pf --output-m`

GridKit's KINSOL solver correctly solves the power balance equations for all
bus voltages and angles. The solved values are internally consistent: the
slack bus P injection and PV-bus Q injections are implicitly determined by
the converged (V, θ) solution. However, `solve_pf.cpp`'s output writer only
extracts V and θ from the bus objects and writes them to `mpc.bus`. It never
performs the post-solve step of computing generator injections from the
converged network state.

**This is a missing output feature, not a missing solve capability.** The
solved (V, θ) contain all the information needed to recover PG and QG. The
post-solve computation is:

For any generator at bus $i$, the net complex power injection at bus $i$ is:

$$S_i = V_i \sum_{j} Y_{ij}^* V_j^*$$

where $Y_{ij}$ are elements of the bus admittance matrix and $V_k = |V_k| e^{j\theta_k}$.
Expanding into real and reactive components:

$$P_i = |V_i| \sum_{j} |V_j| \bigl( G_{ij} \cos\theta_{ij} + B_{ij} \sin\theta_{ij} \bigr)$$

$$Q_i = |V_i| \sum_{j} |V_j| \bigl( G_{ij} \sin\theta_{ij} - B_{ij} \cos\theta_{ij} \bigr)$$

where $\theta_{ij} = \theta_i - \theta_j$. Then:

- **Slack gen PG**: $PG_\text{slack} = (P_i + PD_i) \times \text{baseMVA}$, where $P_i$
  is the net injection and $PD_i$ is the load at the slack bus. This reflects loss
  pickup: the slack absorbs the system real-power mismatch.
- **PV-bus gen QG**: $QG_k = (Q_i + QD_i) \times \text{baseMVA}$, where $Q_i$ is the
  net reactive injection and $QD_i$ is the reactive load at bus $i$. At multi-gen
  buses, QG is split among generators proportional to their Qmax (MATPOWER convention).

Since we wrote `solve_pf.cpp` (see [`pf_helper.md`](pf_helper.md) Section 1 for
architecture), adding this post-solve PG/QG writeback is a straightforward extension.
The bus admittance matrix $Y_\text{bus}$ is already assembled internally by
`SystemModelPowerFlow`; the post-solve step computes bus injections from the
converged state and writes updated PG/QG columns to `mpc.gen` in the output `.m`.

### PV→PQ switching: a separate (and larger) issue

The PG/QG output gap above is purely a writeback issue. A more fundamental
limitation is that GridKit's PF solver does not implement PV→PQ switching
(Q-limit enforcement). This means the solved V and θ themselves may be wrong
when generators hit reactive limits. See [`pf_helper.md`](pf_helper.md)
Section 1 (verified from six independent code paths) and Section 11
(behavioral confirmation under stress).

For the aleatoric UQ pipeline, the perturbation levels (±5-20% per-bus load)
sit inside the "safe envelope" where no additional generators hit Q limits
beyond those already limited at base case (see `pf_helper.md` Section 6: the
~0.030 pu offset is constant across all 14 mild perturbations). At these
levels, the V/θ solution from GridKit is usable for relative comparisons. For
absolute-voltage accuracy, PM.jl remains required.

### Current default and swap path

**Current default**: PM.jl (correct PV→PQ switching, correct PG/QG in output).

**To enable GridKit as an alternative**:
1. Add post-solve PG/QG computation to `solve_pf.cpp` (formulas above).
2. Add a `run_pf_solve(m_path, out_m_path, solver="pm", **kwargs)` dispatcher
   in `pf_utils.py` that routes to either `pm_utils.run_pm_solve_out` or
   `pf_utils.run_solve_pf_out`.
3. Notebook config sets `PF_SOLVER = "pm"` (or `"gridkit"`); all PF calls go
   through the dispatcher.

The swap is meaningful only inside the safe envelope (mild perturbations).
Outside it, PM.jl is required regardless of the PG/QG fix, because the V/θ
solution itself is wrong without PV→PQ switching.

## Field-by-field mapping

Reference: MATPOWER column definitions,
[`cases/illinois.md`](../cases/illinois.md) sections "Bus initial conditions",
"Generator dispatch", "Load demand".

All `case.json` per-unit values are on the system MVA base
(`va_base` in the JSON header, `100.0` for ACTIVSg200).

### 1. Bus initial conditions

For every JSON `bus` with `number == BUS_I` in `mpc.bus`:

```
Vr = VM * cos(VA * pi / 180)
Vi = VM * sin(VA * pi / 180)
```

Written into `bus.init.Vr` and `bus.init.Vi`. `VM` is already per-unit; `VA`
is degrees in MATPOWER and must be converted to radians before `cos`/`sin`.
Any JSON bus with no matching `BUS_I` is left unchanged (this should be zero
for ACTIVSg200 → `illinois.json`, which is 1:1 by construction).

Actual `Bus` device from `illinois.case.json` (bus 49) — only `init.Vr/Vi`
is touched, `name`/`params.kv`/`mon` are never modified:

```json
{
    "number": 49,
    "class": "Bus",
    "name": "RANTOUL 2 1",
    "init": {"Vr": 0.16975389826059847, "Vi": -1.023389141624529},
    "params": {"kv": 13.800000190734863},
    "mon": ["Vr", "Vi"]
}
```

### 2. Generator dispatch (online)

For every JSON device with `class == "Genrou"`, decode
`id = "{bus}_{rank}"` (see [`cases/illinois.md`](../cases/illinois.md) —
"Generator dispatch"). Find the row in `mpc.gen` at that bus with that rank
(ranks are 1-based among generators at the same bus, in the order they appear
in `mpc.gen`). If `GEN_STATUS == 1`:

```
p0 = PG / baseMVA
q0 = QG / baseMVA
```

Written into `Genrou.params.p0` and `Genrou.params.q0`.

See [`cases/illinois.md` — example JSON pieces (bus 49, id `49_1`)](../cases/illinois.md#example-json-pieces-bus-49-id-49_1)
for the full `Genrou` + `SexsPti` device pair, including the dynamics
parameters (`H`, `Xd`, ...) that this mapping stage never touches (see
[§5 below](#5-what-is-not-touched)).

### 3. Offline generator handling (new; not in `patch_case_from_m`)

If `GEN_STATUS == 0` in the `.m` **and** a matching Genrou exists in the base
JSON, remove the `{Genrou, Tgov1, companion-exciter}` device triplet from
`case["devices"]`.

**Companion classes recognised**: `Tgov1`, `SexsPti` (Illinois), `Ieeet1`
(Hawaii), plus `Gast`, `Hygov` for future cases.

**Companion identification**: a companion `dev` is matched to a Genrou with
id `gid` iff:
- `dev.class` is in the companion-class list, AND
- `dev.id == gid` OR `dev.id.startswith(gid + "_")`.

The `dev.id.startswith(gid + "_")` prefix rule handles the observed
`"49_1_sexs_pti"` naming pattern in Illinois. The `dev.id == gid` case is
present because some GridKit cases give companions the same id as the Genrou
(bus-1 convention for Tgov1 in some Hawaii versions).

**Documented limitation**: if the `.m` marks a generator online
(`GEN_STATUS == 1`) but no matching Genrou exists in the base JSON, this
requires a full device-list rebuild (create the Genrou, Tgov1, and exciter
from scratch with reasonable defaults). That is deferred to
`aleatoric_setup.ipynb`. In `build_case_from_solved_m`, this condition
raises `ValueError` when `strict=True` (default) so the caller notices
immediately, or emits a warning when `strict=False`.

For the prototype (`aleatoric_helper.ipynb`), gen-off perturbations are
constrained to gens already online in `illinois.json`. That sidesteps the
rebuild until production.

### 4. LoadZIP demand (new; not in `patch_case_from_m`)

For every JSON device with `class == "LoadZIP"`, group by `ports.bus`. For
each bus with load data in the `.m` (indexed by `BUS_I`):

```
Pnom_pu = PD / baseMVA
Qnom_pu = QD / baseMVA
```

> **Do not write `Vnom`.** `Vnom` is *not* a `LoadZIP` JSON parameter. GridKit
> derives it internally from the connected bus's initial voltage magnitude.
> Writing it into `case.json` has no effect on the simulation and emits one
> `[ERROR]` log line per device. See
> [The `Vnom` non-parameter](#the-vnom-non-parameter-verified) for the full
> investigation. `m_to_case.py` currently still writes it; that is a bug to
> remove.

Because GridKit sets `Vnom = |V|` at the bus and `G = Pnom / Vnom²`, the load's
initial consumed power is `G·Vnom² = Pnom` for any ZIP mix. So `Pnom` **is** the
initial load power by construction, and patching `bus.init.Vr/Vi` and `Pnom`
independently stays self-consistent with no extra coupling.

Actual single-device `LoadZIP` from `illinois.case.json` (bus 2 — note there is
no `Vnom`, and there never should be):

```json
{
    "class": "LoadZIP",
    "ports": {"bus": 2},
    "id": "load_2_1",
    "params": {"Pnom": 0.10819626, "Qnom": 0.03083585, "alphaI": 0.0, "alphaP": 0.0}
}
```

Then split the LoadZIP devices at each bus into two groups:

- **Real loads**: `Pnom > 0` OR (`Pnom == 0` AND `Qnom >= 0`), AND `id` does
  not start with `"shunt_"`.
- **Shunts / capacitors**: `Pnom == 0` AND `Qnom < 0`, OR `id.startswith("shunt_")`.

Shunt-like devices are excluded from PD/QD patching (they represent reactive
compensation, not real load) — their `Pnom`/`Qnom` are left exactly as found.
Nothing else about them is touched; their ZIP reference voltage tracks the new
bus voltage automatically because GridKit re-derives `Vnom` from `bus.init` at
initialization.

There are **4** shunt-like `LoadZIP` devices in `illinois.json`, not just the
one at bus 15 used as the worked example in [§4f](#4f-shunt-exclusion) below:
`shunt_15_1` (Qnom=-0.3233), `shunt_95_1` (Qnom=-0.3143), `shunt_100_1`
(Qnom=-0.8443), `shunt_194_1` (Qnom=-0.5307). Example from the base JSON:

```json
{
    "class": "LoadZIP",
    "ports": {"bus": 15},
    "id": "shunt_15_1",
    "params": {"Pnom": 0.0, "Qnom": -0.32326737, "alphaI": 0.0, "alphaP": 0.0}
}
```

For the real loads at each bus:

- **N=1**: straightforward — set `Pnom=Pnom_pu`, `Qnom=Qnom_pu`.
- **N≥2 (multi-ZIP)**: proportional scaling to preserve each device's share of
  the bus load:

  ```
  S_P = sum(Pnom_old for devs)      # existing sum in base JSON
  S_Q = sum(Qnom_old for devs)
  for each dev:
      Pnom_new = Pnom_old * (Pnom_pu / S_P)
      Qnom_new = Qnom_old * (Qnom_pu / S_Q)
  ```

  This preserves the ratio between multiple loads at the same bus (Illinois has
  a number of multi-ZIP buses each with 2 LoadZIP devices; `illinois.md` says
  38, but that count has not been independently re-verified — see the caveat
  in [the identity round-trip test](#the-identity-round-trip-test) above).

- **Zero-load bus (PD=QD=0)**: set all real LoadZIP `Pnom=0, Qnom=0`.
- **Degenerate case (S_P=0 but PD>0)**: distribute PD equally across the N
  devices; emit a warning. Should not occur in the Illinois base case.

`alphaI` and `alphaP` are **never** changed — they encode the load's ZIP mix,
which is a modeling choice, not a per-scenario quantity.

### 5. What is NOT touched

- Dynamic parameters on Genrou (`H`, `D`, `Xd`, `Xdp`, `Ra`, ...).
- Companion device parameters (Tgov1 gains, SexsPti/Ieeet1 gains).
- Branch data (`R`, `X`, `G`, `B` on Branch devices).
- Fault specification (`BusFault` devices, `solver.json`).
- Monitor lists.
- Any base-case network topology.

## The identity round-trip test

**Setup**: solve `case_ACTIVSg200.m` with PM.jl → `base_pm_solved.m`. Feed that
into `build_case_from_solved_m(illinois.json, base_pm_solved.m, out.json)`.

**Expected**:

- Structural: `out.json` and `illinois.json` have the same set of devices
  (same classes, same ids, same `ports`).
- Numerical: `bus.init.Vr/Vi` and `LoadZIP.{Pnom,Qnom}` and
  `Genrou.{p0,q0}` may differ, but only within tolerance.

**Tolerance philosophy**: the base `illinois.json` was built from a *different*
PF solver (per illinois.md, "The JSON p0/q0 were derived from a different PF
solution (likely PowerWorld)"). So the identity round-trip does not compare
against a mathematical truth; it compares against a *different-solver
equilibrium*. Two sources of residual:

1. **Solver-equilibrium difference**: PM.jl (Ipopt-based) may converge to a
   slightly different (still-valid) equilibrium than PowerWorld. Expect
   `|dVm|` on the order of `1e-3` to `1e-2` pu depending on network
   ill-conditioning; `|dVa|` on the order of `1e-2` to `1e-1` deg. See
   `pm_helper.md` cross-comparison for the observed magnitudes.
2. **Load-model residual**: illinois.md shows JSON `Pnom` is systematically
   ~1.49x higher than `PD/100` from the base `.m`. After rebuilding via
   `build_case_from_solved_m`, the sum-Pnom-per-bus will match `PD/baseMVA`
   exactly. So the identity round-trip is expected to *change* LoadZIP values;
   the invariant to check is that the multi-ZIP proportion between devices at
   the same bus is preserved.

Concrete tolerances used in the notebook:

- Bus count updated: exactly `n_buses` in JSON.
- Genrous updated: **38**, not all 40 present in the base JSON. Buses **161**
  and **197** are `GEN_STATUS=0` in the base `.m` (verified in
  [`cases/illinois.md` — Discrepancy vs. .m file](../cases/illinois.md#discrepancy-vs-m-file))
  but are present in `illinois.json` with non-zero `p0/q0` — since the identity
  round-trip PM.jl-solves the *same* base `.m` (PM.jl does not change
  `GEN_STATUS`), both are caught by the offline-removal rule below rather than
  patched.
- Genrous removed: **2** (`161_1`, `197_1`) — 6 devices total once their
  `Tgov1`/`SexsPti` companions are included. This is the concrete case that
  motivates `remove_offline_gens` existing at all; see
  [What we already know about the target cases](#what-we-already-know-about-the-target-cases)
  above.
- LoadZIP updated: 164 total `LoadZIP` devices in the base JSON = **160 real
  loads + 4 shunts**. The real loads sit on **108 buses**, of which **36** carry
  more than one device (22 buses with 2, 12 with 3, 2 with 4) — verified against
  the `.aux` `Load` table, which has the identical multiplicity distribution.
  `cases/illinois.md`'s "109 buses / 38 multi-ZIP" is wrong and also wrongly
  assumes multi-ZIP means exactly 2.
- Shunts skipped: **4** (`shunt_15_1`, `shunt_95_1`, `shunt_100_1`,
  `shunt_194_1`) — not just bus 15. Bus 15 is used as the single worked
  example in §4f below, but it is one of four, not the only one.
- `Vnom` is **not** a `LoadZIP` parameter and must never appear in `out.json`.
  The base file has zero `"Vnom"` hits and so should every file this pipeline
  produces. Assert `"Vnom" not in dev["params"]` for all 164 LoadZIP devices.
- Any structural change in `out.json` vs `illinois.json` **other than** the
  four field families above → hard fail.

**Caveat**: `cases/illinois.md` predates these checks and is stale in several
places — it says "exactly 1 shunt (bus 15)", gives "109 buses / 38 multi-ZIP",
and lists `Vnom` as a `LoadZIP` parameter. Treat this document as the source of
truth until that one is corrected.

## Perturbation-driven tests (non-identity)

Beyond the identity round-trip, the notebook tests each perturbation family
individually to confirm the code path is exercised and the right fields
change.

### 4a. Bus init only

Manually edit one bus's `VM`, `VA` in a copy of the base `.m` (without solving),
run `build_case_from_solved_m` with `patch_gen_dispatch=False`,
`patch_load=False`, `remove_offline_gens=False`. Verify:

- Only the target bus's `init.Vr/Vi` changed in the output JSON.
- The magnitude of the change matches the trig-derived expected value to
  double-precision.

### 4b. Genrou dispatch only

Manually edit one gen's `PG` in a copy of the base `.m` (no PF solve), run
with only `patch_gen_dispatch=True`. Verify only the target Genrou's `p0`
changed by the expected `dPG / baseMVA`.

### 4c. Offline-gen removal

Set `GEN_STATUS=0` on one gen that IS in the base JSON. Run with
`remove_offline_gens=True`. Verify:

- The `{Genrou, Tgov1, SexsPti}` triplet with that id prefix is removed.
- No other device is affected.
- Total device count decreased by exactly 3 (Illinois).

Then repeat with `strict=True` on a gen that's offline in the base `.m` but
whose companion is **online** in the base JSON — should raise `ValueError`
(the online→JSON-absent case, which is the rebuild-required direction).

### 4d. LoadZIP single-bus patch

Manually edit one single-LoadZIP bus's `PD` and `QD` in the `.m`. Verify:

- Only that bus's LoadZIP `Pnom` and `Qnom` changed.
- `alphaI`, `alphaP` unchanged.
- No `Vnom` key is written (see
  [The `Vnom` non-parameter](#the-vnom-non-parameter-verified)).

### 4e. LoadZIP multi-ZIP proportional scaling

Pick a bus with 2 LoadZIP devices (e.g. bus 39, per illinois.md). Change its
`PD` by a factor of 2 in the `.m`. Verify:

- `sum(Pnom_new) == PD_new / baseMVA` exactly (to floating-point tolerance).
- The ratio `Pnom[0] / Pnom[1]` is preserved to within `1e-12`.

### 4f. Shunt exclusion

Bus 15 in illinois.json has one LoadZIP with `Pnom=0, Qnom<0` (capacitor).
Run with `patch_load=True`. Verify:

- Bus 15's LoadZIP `Pnom` still `0`, `Qnom` still the original negative
  value (unchanged).
- The device is byte-identical to the base entry — nothing is written to it
  at all. Repeat for the other three shunts (`shunt_95_1`, `shunt_100_1`,
  `shunt_194_1`); bus 15 is one of four, not the only one.

### 4g. Zero-load bus

Force `PD=QD=0` on a bus that has one or more LoadZIPs. Verify:

- All real LoadZIPs at that bus have `Pnom=0, Qnom=0`.

### 4h. Full end-to-end with a real perturbation

Take `make_load_scale_m(base_m, 0.8)` → PM.jl solve → run
`build_case_from_solved_m`. Verify:

- Sum of `LoadZIP.Pnom` in the output ≈ `0.8 * sum(base LoadZIP Pnom)` within
  PM.jl solve tolerance.
- Sum of `Genrou.p0` in the output ≈ base sum (slack absorbs the ~20% cut in
  ACTIVSg200; the actual distribution depends on PM.jl's equilibrium choice).

## Reference: MATPOWER columns used

| Field | Column index | Units | Used for |
|---|---|---|---|
| `mpc.bus[:, 1]` (`BUS_I`) | 1 | int | JSON bus `number` |
| `mpc.bus[:, 3]` (`PD`) | 3 | MW | LoadZIP `Pnom` |
| `mpc.bus[:, 4]` (`QD`) | 4 | MVAr | LoadZIP `Qnom` |
| `mpc.bus[:, 8]` (`VM`) | 8 | pu | Bus init `Vr, Vi` |
| `mpc.bus[:, 9]` (`VA`) | 9 | deg | Bus init `Vr, Vi` |
| `mpc.gen[:, 1]` (`GEN_BUS`) | 1 | int | Genrou id decode |
| `mpc.gen[:, 2]` (`PG`) | 2 | MW | Genrou `p0` |
| `mpc.gen[:, 3]` (`QG`) | 3 | MVAr | Genrou `q0` |
| `mpc.gen[:, 8]` (`GEN_STATUS`) | 8 | 0/1 | offline removal |
| `mpc.baseMVA` | header | MVA | per-unit conversion |

## Where `illinois.case.json` actually comes from (verified)

This matters because it explains why nothing in the base JSON lines up with
`case_ACTIVSg200.m`, and it is the context for the `197_1` anomaly below.

**`illinois.case.json` was built from the PowerWorld `.aux` plus the PSS/E
`.dyr`, not from `case_ACTIVSg200.m`.** The `.m` is a parallel, independent
export of a *different* operating point. Source files, all in
`/kfs2/projects/scidac/scidac-data/ACTIVSg200/raw-tamu-data/`:

| File | Supplies |
|---|---|
| `ACTIVSg200.aux` | `Bus` (PUVolt/Angle), `Gen` (MW/Mvr setpoints, status, MVA base), `Load` (S/I/Z split), `Shunt` |
| `ACTIVSg200_dynamics.dyr` | `GENROU` / `SEXS` / `TGOV1` records for all 49 generators |
| `ACTIVSg200_dynamics.dyd` | same dynamics data, PowerWorld format |
| `case_ACTIVSg200.m` | a different operating point; **not** the JSON's source |

Verified correspondences (all computed, not assumed):

| JSON field | `.aux` source | Result |
|---|---|---|
| `Bus.init.Vr/Vi` | `Bus.BusPUVolt`, `Bus.BusAngle` | 300/400 values bit-exact double, 100/400 within `1e-7`, zero mismatches |
| `Genrou.p0/q0/mva` | `Gen.GenMWSetPoint`/`GenMvrSetPoint` ÷ 100, `Gen.GenMVABase` | **40/40 exact** |
| which 40 Genrou exist | `Gen.GenStatus` | 40 `"Closed"`, 9 `"Open"` — exactly the 40 present / 9 absent |
| `LoadZIP` structure | `Load` (160 records, 108 buses) | per-bus multiplicity `{1:72, 2:22, 3:12, 4:2}` — **identical** in both. 164 LoadZIP = 160 loads + 4 shunts |
| shunt `Qnom` | `Shunt.SSNMVR` | `Qnom = -(SSNMVR/100)·|V|²` holds for all 4 (residuals `7e-9`…`3e-8`) |

`Pnom`/`Qnom` match `Load.LoadSMW`/`LoadSMVR` ÷ 100 to a worst residual of
`5.0e-09` across all 320 values — that is the JSON writer's ~8-decimal output
precision, not a data discrepancy.

Bus counts corrected: the `.aux` has **108** load buses and **36** multi-ZIP
buses (22 with 2 devices, 12 with 3, 2 with 4). `cases/illinois.md`'s
"109 buses / 38 multi-ZIP" is wrong, and it assumed multi-ZIP meant exactly 2.

### Device parameter mappings (proven from the `.dyr`)

```
GENROU  = [Tdop, Tdopp, Tqop, Tqopp, H, D, Xd, Xq, Xdp, Xqp, Xdpp, Xl, S10, S12]
          -> JSON Genrou, same names.  Xqpp = Xdpp  (PSS/E GENROU has one X").
          -> mva = MBASE = aux GenMVABase.
SEXS    = [TA/TB, TB, K, TE, Emin, Emax]
          -> JSON SexsPti: Ta = dyr[0]*dyr[1], Tb, K, Te, Efdmin, Efdmax.
TGOV1   = [R, T1, Vmax, Vmin, T2, T3, Dt]
          -> JSON Tgov1: R, T1, Pvmax=Vmax, Pvmin=Vmin, T2, T3, Dt, Trate=MBASE.
```

39 of 40 Genrou devices reproduce the `.dyr` exactly when compared at the
`.dyr`'s own 4-decimal print precision. The JSON carries **more** precision
than the `.dyr`/`.dyd` text (e.g. `49_1` `H`: `.dyr` prints `5.2115`, JSON
stores `5.21145296`, which rounds to `5.2115`), so the full-precision values
came from the PowerWorld side, not from the rounded text exports.

Because the `.dyr` covers **all 49** generators — including all 9 absent from
the JSON (78, 79, 92, 164, 165, 166, 168, 169, 196) — adding a generator for a
unit-commitment scenario needs no invented defaults for this case. A `.dyr`
parser plus the three mappings above is sufficient.

### The `197_1` anomaly (the one departure from the `.dyr`)

`197_1` is the **one** Genrou whose parameters do not come from the `.dyr`, and
it differs in **all 14** of them:

| param | `.dyr` (bus 197) | `illinois.case.json` |
|---|---|---|
| `H` | 2.0687 | **3.0** |
| `Tdop` / `Tdopp` | 9.5968 / 0.0306 | **7.0** / **0.04** |
| `Tqop` / `Tqopp` | 1.3646 / 0.0422 | **0.75** / **0.05** |
| `Xd` / `Xq` | 1.5534 / 1.4877 | **2.1** / **0.5** |
| `Xdp` / `Xqp` | 0.2436 / 0.2218 | **0.2** / **0.5** |
| `Xdpp` (= `Xqpp`) | 0.1773 | **0.18** |
| `Xl` | 0.132 | **0.15** |
| `S10` / `S12` | 0.1888 / 0.6127 | **0.0** / **0.0** |

The JSON values are round numbers and the saturation terms are zeroed — it is
the only machine in the fleet with `S10 = S12 = 0`. These look like generic
GENROU defaults rather than TAMU data, but **why the case creator departed from
the `.dyr` for this one unit is not known** and should be asked directly.

Scoping facts:

- Its **companions are correct**. `197_1_sexs_pti` matches `.dyr` `SEXS`
  `[0.1, 10, 250, 0.1, -4, 5]` (`Ta = 0.1 × 10 = 1.0`), and `197_1_tgov1`
  matches `.dyr` `TGOV1` `[0.04, 0.5, 1, 0, 0.918892, 2.7567, 0]` with
  `Trate = 81 = MBASE`. Only the `Genrou` device is affected.
- `161_1` — the other generator that is `GEN_STATUS=0` in the `.m` but present
  in the JSON — matches the `.dyr` perfectly. So this is not a
  "both anomalous gens are hand-filled" pattern; it is specific to `197_1`.
- The discrepancy is not cosmetic for UQ purposes: `H` is overstated by ~45%
  (3.0 vs 2.0687) on an 81 MVA machine, and zeroed saturation changes the
  field-voltage response.

**Pipeline rule**: the `.dyr` is the default. `197_1` is rebuilt from the `.dyr`
like every other machine, so the pipeline has one self-consistent provenance
chain with no per-device exceptions baked into the code. The case creator's
values are recoverable through the **creator overlay** below.

### The creator overlay

`197_1` is unlikely to be the only place a case author departed from the source
data — that judgement is part of the craft of building a dynamics case, and
Hawaii and the larger ACTIVSg cases may well have their own. Rather than
special-casing each one, the deviations are captured as data:

```
creator_overlay.json
{
  "197_1": {
    "H": 3.0, "Tdop": 7.0, "Tdopp": 0.04, "Tqop": 0.75, "Tqopp": 0.05,
    "Xd": 2.0999999, "Xdp": 0.2, "Xdpp": 0.18000001,
    "Xq": 0.5, "Xqp": 0.5, "Xqpp": 0.18000001,
    "Xl": 0.15000001, "S10": 0.0, "S12": 0.0
  }
}
```

The overlay is **generated, never hand-written**: rebuild the case from
`.aux` + `.dyr`, diff against the shipped `illinois.case.json`, and record every
parameter that differs by more than the JSON write precision. That keeps it
automatically in sync — if GridKit ships an updated case, re-running the diff
surfaces any new or removed deviations instead of silently diverging.

Builder behaviour:

| `creator_overlay` | Result |
|---|---|
| `None` (default) | pure `.aux` + `.dyr`; `197_1` gets its real TAMU parameters |
| path / `"auto"` | overlay applied; output reproduces the shipped `illinois.case.json` |

Both modes must be recorded in the run metadata, because they produce
physically different machines at bus 197 (`H` 2.0687 vs 3.0 on an 81 MVA unit).

**Consequence for the round-trip test**: with the default (`None`), the rebuilt
case will *not* match the shipped file — it will differ in exactly the 14
parameters above and nowhere else. The identity test therefore has two arms:
run with the overlay to assert byte-equivalence, and run without it to assert
that the diff is *exactly* the known overlay and contains no surprises.

## The `Vnom` non-parameter (verified)

Short version: `Vnom` is not a `LoadZIP` input. It never was. It entered this
project's own docs and code by mistake and should be removed from both.

### Where it came from

Not from any case file — `illinois.case.json` and `hawaii.case.json` both have
**zero** `"Vnom"` occurrences. It originates in
[`uq_plan.md`](uq_plan.md#L285), which lists the GridKit `LoadZIP` parameters as
`Pnom, Qnom, Vnom, alphaI, alphaP`. That reads `Vnom` out of the model
*equations* (`G = Pnom/Vnom²`) and mistakes it for a model *input*.
`m_to_case.py` then implemented that reading, and earlier revisions of this
document and `cases/illinois.md` propagated it.

### What GridKit actually does

From `GridKit/Model/PhasorDynamics/Load/LoadZIP/`:

```
LoadZIPData.hpp   enum LoadZIPParameters { Pnom, Qnom, alphaI, alphaP };   // no Vnom
LoadZIPImpl.hpp   Vnom_ = vm0;                       // derived from bus init voltage
LoadZIPImpl.hpp   G_ = Pnom_/Vnom2;  B_ = Qnom_/Vnom2;
LoadZIPImpl.hpp   alphaZ_ = 1 - alphaI_ - alphaP_;
```

`ComponentDataJSONParser.hpp` whitelists every `params` key through
`magic_enum::enum_cast<Parameters>`. An unrecognised key takes the `else`
branch and emits `Log::error() << "Initial parameter \"...\" has no value."`.
`Logger.cpp` shows `error()` only writes to a stream (default `std::cout`) — it
does not throw or change the exit code.

### Empirical confirmation

Injected `"Vnom": 0.9998` into `load_2_1` and ran `DynamicSimulation`
(`tmax = 0.05`) against the unmodified case:

| run | exit | `[ERROR]` lines | `mon.csv` |
|---|---|---|---|
| base | 0 | 0 | 27,443,596 bytes |
| with `Vnom` | 0 | **1** | 27,443,596 bytes — **byte-identical** |

Message emitted:
`Initial parameter "Vnom" has no value. See the "LoadZIP" device with "id": "load_2_1" ...`

So the value is ignored completely, and the only effect is log noise — one line
per device, which would be 164 lines for Illinois if every LoadZIP were patched.
Because it goes to **stdout**, it would also be captured into `stdout.txt` by
`run_sample` and could false-positive any QC step that greps logs for `ERROR`.

Hawaii cannot be used to cross-check this: it uses `LoadZ`, whose parameter enum
is `{R, X}` only.

### Consequences

1. `m_to_case.py` must stop writing `Vnom` (currently written in five places:
   the shunt branch, the zero-load branch, and the three scaling branches).
2. `alphaI`/`alphaP` with both at `0.0` mean `alphaZ = 1` — **100% constant
   impedance**, not constant power. The ZIP mix is a case-level modelling
   choice and the pipeline preserves it verbatim.
3. *(unverified)* The `.dyd` `lodrep` section holds one system-wide
   `wlwscc  ... 0 1 1 0 0 0 0 0 0 0 0` record. Fields 2–3 being `1 1` would be
   consistent with 100% constant-Z for P and Q, which would justify
   `alphaI = alphaP = 0`. PowerWorld's `wlwscc` field order is not known to the
   author, so this is a working hypothesis, not a confirmed provenance.

## Known limitations

Ordered by when they will bite. The first two are guarded so they cannot fail
silently; the rest are documented behaviour.

### 1. Multiple generators per bus — guarded, not yet supported

`parse_dyr` keys on `(bus, model)` and **drops the PSS/E unit id**. A bus with
two machines of the same model would keep only the last record. Likewise
`gen_ids()` builds `"{bus}_{rank}"` where rank is the **positional** order of
rows in `mpc.gen`; if an upstream tool (Sienna, a PCM export) reorders the rows
of a multi-unit bus, dynamics would attach to the wrong machine.

Neither bites on ACTIVSg200 — verified: 49 generators on 49 distinct buses,
every PSS/E unit id is `'1'`, zero `(bus, model)` collisions. **Both will bite on
ACTIVSg2000 and ACTIVSg10k.**

Guards are in place so this cannot happen unnoticed:

- `parse_dyr` **raises** on any `(bus, model)` collision, naming the buses.
  `allow_multi_unit=True` is a deliberate last-record-wins override.
- `gen_ids` **warns** when any bus carries more than one generator.

The fix is to key both on the unit id, `(bus, model, unit)` and
`"{bus}_{unit}"`. Deferred until a case actually needs it, and it should be done
together with the Sienna handoff, since that is where row ordering stops being
under our control.

### 2. Switching off a generator overloads the slack

Bus 189 sits at PG = 569.14 MW against PMAX = 569.1. In a plain PF the slack
absorbs the entire output of anything switched off. Removing a 150–200 MW unit
drives the slack 25–35% past its rating, and its Genrou `p0` with it.

Until `perturb.py` grows a redistribution rule (share the lost MW across online
units by headroom), keep test gen-off cases small. `65_1` (PAXTON 1 1) is the
useful probe: H×mva = 1135 MWs — 7.2% of the 15,694 MWs fleet total — but only
4 MW of dispatch, so it changes inertia substantially while barely moving the
power balance.

### 3. Adding a generator is untested

`_add_gens_from_dyr()` exists and builds `{Genrou, exciter, governor}` triplets
plus their `signals` entries from `.dyr` records. Every round trip so far has
reported `gen_added: 0`, so **this path has never executed against real data**.
Commitment is fixed during bootstrap in the current workflow, so it is not on
the critical path; treat it as unvalidated until a UC study needs it.

### 4. `LoadZ` (Hawaii) not implemented

Illinois uses `LoadZIP`; **Hawaii uses `LoadZ`** with params `R`, `X` — not a
class named `Load`, as `cases/hawaii.md` currently states. The patch algorithm
differs: `G = PD_pu / VM²`, `B = -QD_pu / VM²`, `R = G/(G²+B²)`, `X = -B/(G²+B²)`.

### 5. Multi-ZIP with more than one shunt at a bus

Shunts are identified and excluded from PD/QD scaling, but the code assumes at
most one shunt-like device per bus. Several shunts on one bus would need
reactive splitting. Not observed in Illinois (4 shunts, 4 distinct buses).

### 6. ZIP mix is preserved verbatim

`alphaI`/`alphaP` are never modified. A scenario needing a different ZIP mix
belongs in a separate device-model editor, not here.

### 7. Bugs in the superseded `m_to_case.py`

Both are live in the old module, which `m_to_case_utils.py` replaces:

- `from matpowercaseframes import Case` — the installed v2.1.0 exports
  `CaseFrames` only, so `_load_m` raises. The `try/except` rewrites it into a
  misleading "pip install matpowercaseframes" message for a package that is
  already installed. This means `patch_case_from_m` and
  `build_case_from_solved_m` have never been runnable in `h-py312-basic`.
- It writes `Vnom` at five call sites — see
  [The `Vnom` non-parameter](#the-vnom-non-parameter-verified).

Only `m_to_case_helper.ipynb` imports it; `uq_setup.ipynb` does not, so the
completed epistemic runs are unaffected.


## Findings

### Verified (2026-09-29, independent of the notebook)

| Finding | Result |
|---|---|
| Provenance of `illinois.case.json` | PowerWorld `.aux` + PSS/E `.dyr`, **not** `case_ACTIVSg200.m` |
| `Genrou.p0/q0/mva` vs `.aux` `Gen` | 40/40 exact |
| Bus `init.Vr/Vi` vs `.aux` `Bus` | 300/400 bit-exact double, 100/400 within `1e-7`, 0 mismatches |
| `LoadZIP.Pnom/Qnom` vs `.aux` `Load` | 320/320 agree, worst residual `5.0e-09` (JSON write precision) |
| Shunt `Qnom = -(SSNMVR/100)·\|V\|²` | 4/4, residuals `6.8e-9` … `2.7e-8` |
| `Genrou` dynamics vs `.dyr` | 39/40 exact at the `.dyr`'s 4-dp print precision; `197_1` is the sole exception |
| `Vnom` in `case.json` | **ignored by the simulator** — `mon.csv` byte-identical, 1 spurious `[ERROR]` line |
| Load-bus counts | 108 load buses, 36 multi-ZIP (`{1:72, 2:22, 3:12, 4:2}`) — corrects `cases/illinois.md` |

### Pending the notebook run

- **Identity round-trip Vr/Vi residual**: [TODO]
- **Identity round-trip LoadZIP change vs base**: [TODO]
- **Genrou removal test on `illinois.case.json`**: [TODO]
- **Multi-ZIP scaling exact-match**: [TODO]
- **End-to-end 0.8x load perturbation Pnom sum**: [TODO]
- **Creator-overlay diff contains exactly `197_1` and nothing else**: [TODO]

## Forward pointers

- `aleatoric_helper.md` (to come): consumes `build_case_from_solved_m` end-to-end
  (perturb → PM.jl → build case.json → DynamicSimulation → QoI compare).
- `aleatoric_setup.md` / `aleatoric_setup.ipynb` (to come, production):
  scales the workflow to N hours of user-provided load/wind timeseries; adds
  device-list rebuild and kill-list builder.
