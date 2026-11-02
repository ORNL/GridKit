# m_to_case helper: source data → GridKit `case.json`

Design notes for the case-construction pipeline behind
[`notebooks/m_to_case_helper.ipynb`](../notebooks/m_to_case_helper.ipynb) and
[`notebooks/8760_aleatoric.ipynb`](../notebooks/8760_aleatoric.ipynb).

Last substantive revision: 2026-10-19, after the upstream `Cleanup Examples`
(#577) restructure replaced the Illinois example with a differently-shaped
ACTIVSg200 case.

---

## Scope

Build GridKit `case.json` files whose **operating point** we control, so a
dynamic simulation can be run for an arbitrary load / generation scenario.

Two notebooks, deliberately split:

| Notebook | Owns |
|---|---|
| `m_to_case_helper.ipynb` | conversion, validation, **base case construction** |
| `8760_aleatoric.ipynb` | consuming those base cases for a large campaign |

The split matters: there is one definition of what a base case is, it is
validated where it is built, and the campaign notebook cannot silently drift
from it.

---

## What each source file contributes

The two TAMU files are complements, not alternatives. Neither alone is enough.

| File | Supplies | Analogous to |
|---|---|---|
| `.aux` | topology, bus voltages, load, dispatch, **which units are online** | a MATPOWER `.m` |
| `.dyr` | machine dynamic models: `GENROU`, `SEXS`, `TGOV1`, **147 records = 49 × 3** | nothing in `.m` |

A MATPOWER `.m` has nowhere to put an inertia constant, so the `.m` is
power-flow only. Dynamics always come from the `.dyr`, directly or indirectly.

---

## The chain

```
   BOOTSTRAP  (once per case)
   ───────────────────────────────────────────────────────────────────
     .aux ─┐
           ├──> BASE_CASE (case.json) ──────────────────────┐
     .dyr ─┘         validated: 0 diffs vs shipped          │
                                                            │ dynamics
     .aux ──> seed .m                                       │ inherited
                 │                                          │ from here
   ──────────────┼────────────────────────────────────────  │
   CASE MOD      │  structural, a handful, made in advance   │
   (helper §8)   │  gen on/off  ->  changes the device list  │
                 ▼                                          │
            base_<c>.m ──[OPF]──> base_<c>_solved.m         │
                 │                                          │
   ──────────────┼────────────────────────────────────────  │
   PERTURBATION  │  parametric, 8760, one per hour           │
   (8760 nb)     │  load / wind  ->  changes values only     │
                 ▼                                          │
            hour_<h>.m ──[OPF or PF]──> hour_<h>_solved.m    │
                                              │             │
                                              ▼             ▼
                                   build_case_from_m ───────┘
                                              │
                                        case.json
                                              │
                                     DynamicSimulation ──> mon.csv
```

`[OPF]` **is implemented and works**: `pm_solve.jl` exposes
`--pf-type ac|dc|ac-opf|dc-opf`. Measured on ACTIVSg200: DC OPF solves in
0.026 s, AC OPF in 3.02 s, AC PF in 2.94 s, all `LOCALLY_SOLVED`.

After bootstrap nothing reads PowerWorld. Everything downstream speaks MATPOWER,
which is also what Sienna and PowerModels.jl speak.

### What the production path overwrites

`build_case_from_m(SPEC, BASE_CASE, solved_m)` deep-copies `BASE_CASE` and
replaces only the operating point:

| Overwritten from the `.m` | Inherited from `BASE_CASE` |
|---|---|
| bus `init.Vr`/`Vi` (from `VM`/`VA`) | Genrou `H`, `Xd`, `Tdop`, ... |
| Genrou `p0`/`q0`/`mva` | SexsPti and Tgov1 parameters |
| `LoadZIP` `Pnom`/`Qnom` | branches, `BusFault`, signals, monitors |

So **the `.dyr` is normally not read again**. Dynamics reach a scenario case by
inheritance, not re-derivation, which guarantees every scenario in a campaign
shares byte-identical machine parameters. Any difference in results is then
attributable to the operating point.

---

## Case mods vs perturbations

Two kinds of change. The test is **does it change the device list in
`case.json`?**

| | **Case mod** | **Perturbation** |
|---|---|---|
| changes | the device list | parameter values only |
| example | generator on/off | bus load, wind output |
| how many | a handful, pre-made | 8760, one per hour |
| function | `perturb.set_gen_status` | `perturb.set_bus_load`, `set_gen_output` |
| needs the `.dyr`? | only to **recommit** an `.aux`-offline unit | never |
| changes system inertia? | **yes** | no |
| changes `n_genrou`? | **yes** | no |
| moves `element_id`? | no (one `BusFault`, always 0) | no |
| lives in | helper notebook §8 | `8760_aleatoric.ipynb` |

**Inertia only moves with a case mod.** A perturbation changes injections at
fixed inertia. The campaign is therefore `case mod × perturbation`: mods set the
stiffness of the system, perturbations sample aleatoric variability within it.

### Commitment is asymmetric

`BASE_CASE` materialises only the **40 units online in the `.aux`**, not all 49.

* **Off** needs nothing extra: the whole triplet is dropped from the copy.
* **On**, for one of the 9 units the `.aux` had offline (`78_1`, `79_1`, `92_1`,
  `164_1`, `165_1`, `166_1`, `168_1`, `169_1`, `196_1`), requires `dyr=`,
  because there is no device to inherit from. Without it `build_case_from_m`
  raises. With it, `_add_gens_from_dyr` appends the full triplet. Verified:
  `gen_added: 1`, ids `78_1_genrou` / `78_1_tgov1` / `78_1_sexs_pti`, 0 machines
  missing companions. Structurally correct; not yet exercised in a dynamic run.

`set_gen_status(..., 1)` sets the status but leaves `PG` at 0, so a recommitted
unit also needs an output assigned.

---

## Where OPF belongs

Both legs need a re-dispatch step, for different reasons. The generated `.m`
already carries everything an OPF needs: `VMIN`/`VMAX` = 0.9/1.1, `RATE_A`
nonzero on all 246 branches (max 740 MVA), `PMIN < PMAX` on all 49 units,
`QMIN < QMAX` on all 49, and quadratic `gencost` (`MODEL 2`, `NCOST 3`,
`c2 = 0.002`, `c1 = 19.0`).

**After a case mod.** Decommitting 450 MW leaves a hole.
`redistribute="headroom"` fills it proportionally to spare capacity: a heuristic
with no economics and no network awareness. An OPF gives a defensible reference
dispatch instead. Once per mod, so cost is irrelevant.

**After a perturbation.** More urgent. The hourly loop currently changes `PD`
and lets the PF slack absorb the whole mismatch, and the slack is nuclear
`189_1` with **0.008 MW** of headroom. A ±15 % load profile pushes it hundreds
of MW past `PMAX`, and PF converges anyway because it does not enforce limits.
An OPF over the *fixed* committed set (economic dispatch, not unit commitment)
is the correct fix.

**An AC OPF does not strictly need a following PF, but one is worth running.**
`solve_ac_opf` uses the full nonconvex ACP formulation, so its solution is an
AC-feasible point in principle, and it returns both the dispatch and the bus
voltages that `build_case_from_m` consumes. A follow-up PF is still a cheap
independent check: Ipopt converges only to a tolerance, the operating point
round-trips through a text `.m`, and the `.m` carries neither bus type nor
which reactive limits were binding. If the PF moves the point materially,
something upstream is wrong.

For a DC OPF the following AC PF is **mandatory**, not optional: DC ignores
losses entirely, so the dispatch under-generates and the PF slack absorbs the
difference. Measured on ACTIVSg200 at peak, that put the slack 27.05 MW over
`PMAX`; one loss-aware iteration reduced it to 0.34 MW.

Keep PF available for three things: the identity checks in helper sections 2–3,
where the point is to solve a *given* dispatch rather than improve it; a
fallback when an OPF is infeasible; and backend comparison.

One caveat: OPF enforces constraints PF ignores, so it will *fail* on operating
points PF would happily return. That is the desired behaviour, but it means the
campaign needs an infeasibility path, not just a convergence check.

---

## Entry points

| Module | Role |
|---|---|
| `m_to_case_utils.py` | `.aux`/`.dyr` parsing, case assembly, diffing, overlays |
| `aux_to_m.py` | MATPOWER frames in/out (`m_to_frames`, `write_matpower_m`) |
| `perturb.py` | scenario edits on frames (load, gen output, commitment) |
| `solver_utils.py` | `solver.json`, fault placement, validation |
| `pm_utils.py` | PowerModels.jl power-flow driver |

Key functions:

* `build_case_from_aux_dyr(spec, base_case, aux, dyr, creator_overlay=None)`
  — bootstrap path. Rebuilds a case from source data.
* `build_case_from_m(spec, base_case, m_path)`
  — production path. Maps a solved `.m` onto `base_case`'s device list.
* `make_creator_overlay(shipped, rebuilt, tol)` / `apply_overlay`
  — capture and re-apply the case author's departures from the `.dyr`.
* `diff_case_json(ref, new, tol)` — device-by-device comparison; an empty frame
  means agreement.

---

## The #577 migration (2026-09-21)

Upstream deleted `examples/PhasorDynamics/Large/Illinois` and introduced
`cases/PhasorDynamics/ACTIVSg200`. This is **not** a rename. The underlying
network is the same ACTIVSg200 system, but the `case.json` is shaped
differently:

| | old `illinois.case.json` | new `ACTIVSg200.case.json` |
|---|---|---|
| buses / Branch / Genrou | 200 / 246 / 40 | 200 / 246 / 40 |
| **LoadZIP** | 164 (160 real + 4 shunt), up to 4/bus | **109, exactly 1/bus** |
| **BusFault** | **200** (one per bus) | **1**, bus 149 `SPRINGFIELD 5 0` |
| device ids | `49_1`, `load_12_1` | `49_1_genrou`, `load_2` |
| Genrou params | no `Ra` | includes `Ra` |
| solver | plain | carries validation keys |

The case now lives in `cases/`, the solver stays beside the example, so
`CaseSpec` grew a `solver_dir` field. `CaseSpec.ILLINOIS` remains as an alias
for `ACTIVSG200`.

### Three porting bugs this caused, all found by diffing

1. **Suffixed device ids.** `49_1_genrou` no longer matched the bare `49_1` from
   `mpc.gen`, so every machine was silently dropped (`genrou: 0`). Fixed with
   `CaseSpec.machine_suffix` and `gen_id_of()`, which maps any of
   `49_1_genrou` / `49_1_tgov1` / `49_1_sexs_pti` / `49_1` to `49_1`.
2. **One load per bus.** The `.aux` has 160 load records across 108 buses. With
   one `LoadZIP` per bus the records must be **summed**, not rank-matched.
3. **Shunt folded into load.** Where a bus has both a load and a shunt, the new
   case merges them into the load's `Qnom`. Verified exactly at buses 95, 100
   and 194:

   ```
   bus  95: load Q=+0.011663  shunt Q=-0.314276  sum=-0.302614  shipped=-0.302614
   bus 100: load Q=+0.168582  shunt Q=-0.844345  sum=-0.675763  shipped=-0.675763
   bus 194: load Q=+0.161215  shunt Q=-0.530665  sum=-0.369450  shipped=-0.369450
   ```

   Bus 15 has a shunt but no load, so it stays a standalone `Pnom = 0` device.

Both builders now support either convention: one device per bus takes the
combined total, several devices per bus keep the old split behaviour.

---

## Validation

The health check in section 2 of the helper notebook is the trust anchor.
Everything downstream inherits from it, and it asserts rather than reports.

### Step 1: the mapping is complete

Rebuild `case.json` from `.aux` + `.dyr` alone, diff against the shipped file,
capture the difference as `OVERLAY`, re-apply, and require **zero** remaining
differences.

```
rebuild counts : {'bus': 200, 'genrou': 40, 'sexs': 40, 'tgov1': 40,
                  'load': 108, 'shunt': 1, 'shunt_folded': 3, 'missing_dyr': []}
overlay        : {'197_1_genrou': 14}
diffs no ovl   : 14
diffs w/ ovl   : 0
```

Zero is the point. It means every difference between source data and shipped
case is *accounted for*, not merely small. The overlay is discovered by diffing,
never hand-written, so it self-corrects if GridKit reships the case.

### Step 2: the seed `.m` carries the operating point

```
seed .m -> case.json, diffs vs BASE_CASE @1e-6: 0
```

### The `197_1` anomaly

Of 49 machines, `197_1` is the only one whose Genrou parameters in the shipped
case do not come from the `.dyr`. All 14 differ:

| param | shipped | `.dyr` |
|---|---|---|
| H | 3.00 | 2.0687 |
| Tdop | 7.00 | 9.5968 |
| Xd | 2.10 | 1.5534 |
| S10 / S12 | 0.00 / 0.00 | 0.1888 / 0.6127 |

Its **companions are consistent** with the `.dyr` (`SexsPti Ta = 0.1 × 10`,
`Tgov1 Trate = MBASE`), and `161_1` matches the `.dyr` perfectly. So this looks
like a deliberate edit to one machine, not a systematic conversion difference.
We have not established *why*, which is exactly why it is switchable rather than
baked in.

The anomaly survived the #577 restructure unchanged, which is reassuring about
both the finding and the rebuild.

### `USE_CREATOR_OVERLAY`

| Setting | `197_1` params | Use when |
|---|---|---|
| `True` | the author's values | reproducing or validating against GridKit |
| `False` | straight from the `.dyr` | the `.dyr` is the authority, e.g. calibration |

Only machine dynamics change. The operating point always comes from the `.aux`,
and the seed `.m` is identical either way, because MATPOWER has nowhere to store
a Genrou parameter.

**`SHIPPED` vs `BASE_CASE`.** `SHIPPED` is the reference file; `BASE_CASE` is
what every scenario inherits dynamics from. Downstream code uses `BASE_CASE`
exclusively. Passing `creator_overlay` to `build_case_from_m` alongside
`SHIPPED` was a no-op, because `SHIPPED` already contains the author's values;
the real lever is which base case is handed downstream.

---

## Identifier spaces

A recurring source of confusion. Five distinct namespaces:

| Id | Example | Lives in |
|---|---|---|
| bus number | `149` | `.aux`, `.m` (`BUS_I`), `case.json` (`buses[].number`) |
| bus name | `"SPRINGFIELD 5 0"` | `.aux`, `case.json` (`buses[].name`) |
| device id | `49_1_genrou` | `case.json` only |
| `element_id` | `0` | `solver.json` only |
| gen id | `65_1` | `<bus>_<unit>`, `perturb` and the `.dyr` |

Verified: bus numbers are **identical sets** across the seed `.m` and
`ACTIVSg200.case.json` (200 each). Bus names round-trip too, and the trailing
digit is part of the name, not the number: bus 1 is `"CREVE COEUR 0"` and bus 2
is `"CREVE COEUR 1"`.

### `element_id` indexes faults only

From `SystemModelImpl.hpp:773`:

```cpp
IdxT gridkit_fault_id = static_cast<IdxT>(gridkit_fault_indices_.size());
```

It counts **`BusFault` devices only**, assigned in registration order. Two
consequences:

* it is not a bus number, and does not track bus order
* **adding or removing generators does not move it**

With ACTIVSg200 shipping exactly one fault, `element_id` is always 0. Verified:
decommitting `126_1` removed 3 devices and left the fault at `element_id` 0,
bus 149.

This is a change from the old case. With 200 faults, commitment had to be frozen
for a campaign or indices would shift. That constraint no longer exists.

Always resolve with `su.fault_element_id(case, bus)`, which accepts a bus number
or name and raises on an ambiguous match. Never type an `element_id`.

---

## Field-by-field mapping

### Bus initial conditions

`.aux` `BusPUVolt` / `BusAngle` → `buses[].init.Vr` / `Vi` as
`V·cos(θ)`, `V·sin(θ)`. From a `.m`, `VM`/`VA` the same way.

### Generator dispatch

`.aux` `GenMWSetPoint` / `GenMvrSetPoint` / `GenMVABase` → `Genrou.params`
`p0` / `q0` / `mva`, on system base. From a `.m`, `PG`/`QG`/`MBASE`.

Offline generators (`GenStatus != Closed`, or `GEN_STATUS == 0`) have their
**whole device triplet removed**: `Genrou`, `SexsPti`, `Tgov1`.

### Machine dynamics, from the `.dyr`

`GENROU` record order is
`[Tdop, Tdopp, Tqop, Tqopp, H, D, Xd, Xq, Xdp, Xqp, Xdpp, Xl, S10, S12]`.
`Xqpp` is set equal to `Xdpp`; `Ra` is emitted as 0.0.

`SEXS`: `Ta = (TA/TB) × TB`. `TGOV1`: `Trate = MBASE`.

Residual differences from the `.dyr`'s 4-decimal rounding show at the 1e-6 level
and are expected.

### Loads and shunts

Per bus. One `LoadZIP` takes the summed load plus the shunt folded into `Qnom`;
several keep the old split with `shunt_` devices separate. Shunt contribution is
`Qnom = -(SSNMVR / baseMVA) · |V|²`, or from a `.m`, `-(BS / baseMVA) · Vm²`.

Shunts are rebuilt from `mpc.bus` `BS`, **never** from `QD`.

### Not touched

Branches, `BusFault` devices, signals, monitors, and `params` are inherited from
the base case unchanged.

---

## The `Vnom` non-parameter

`Vnom` must never be written to a `LoadZIP`. It is not a GridKit parameter.

* `LoadZIPData.hpp` enumerates exactly `{Pnom, Qnom, alphaI, alphaP}`
* `LoadZIPImpl.hpp:144` sets `Vnom_ = vm0` internally, from the initial voltage
* `ComponentDataJSONParser.hpp` whitelists keys via `magic_enum::enum_cast`; an
  unknown key logs `Log::error()` and is **ignored**, without failing

Empirically: injecting `"Vnom": 0.9998` gave exit 0 both ways, ERROR lines 0 vs
1, and a **byte-identical** `mon.csv` (27,443,596 bytes).

So the failure mode is silent. The audit asserts the count is zero, and both
builders pop the key defensively. The idea originated in `uq_plan.md` line 285,
which lists `Vnom` as a LoadZIP parameter; that line is wrong.

---

## Wind is modelled as synchronous machines

ACTIVSg200 has **no inverter-based generation**. This is TAMU's source data, not
a GridKit conversion choice: `ACTIVSg200_dynamics.dyr` contains only `GENROU`,
`SEXS` and `TGOV1`, 49 of each.

All six wind-labelled units carry a full synchronous machine:

| gen | bus | H | mva | H·MVA |
|---|---|---|---|---|
| 65_1 | 65 | 6.29 | 180.5 | 1135.2 |
| 104_1 | 104 | 2.40 | 118.8 | 285.1 |
| 105_1 | 105 | 1.50 | 237.6 | 356.7 |
| 114_1 | 114 | 8.92 | 2.0 | 18.2 |
| 115_1 | 115 | 2.73 | 180.0 | 490.8 |
| 147_1 | 147 | 3.98 | 120.6 | 479.9 |

The `genfuel` label is the only thing that says "wind", and nothing dynamic
reads it. Consequences:

* wind contributes synchronous **inertia** (2,766 MW·s, ~18 % of the fleet)
* wind provides **governor droop** through `Tgov1`
* wind contributes **synchronous short-circuit current** during a fault
* varying wind `PG` changes injection only, never inertia

None of which a Type 3/4 converter would do.

**What this bounds.** An 8760 campaign here measures injection variability *at
fixed inertia*. The dominant real-world pathway, where high wind displaces
synchronous units through unit commitment and inertia falls, is structurally
absent. That pathway is reachable only by varying commitment, which is why the
base case library exists.

**Forward path.** GridKit implements `Regca` (full C++ with Enzyme derivatives)
and `SystemModelDataJSONParser.hpp:138` dispatches `"Regca"`, so a converter can
be instantiated from `case.json` today. But `REECA` and `REGCB` are README-only
stubs; `Regca`'s `ipcmd`/`iqcmd` inputs are optional, so it would run open-loop
on fixed current commands. Enough to remove inertia and synchronous fault
current, not enough for a credible ride-through study.

---

## Commitment and the base case library

Built in section 8 of the helper notebook, consumed by `8760_aleatoric`.

### Thermal only

Wind is never a commitment candidate: it is non-dispatchable, it will change
representation when IBR models land (thermal units will not), and it holds most
of the system's spare capacity.

`perturb.NONDISPATCHABLE_FUELS = {wind, solar, pv, photovoltaic}`.
`perturb.dispatchable_gens()` drops those fuels and the slack, leaving **33
candidates** (25 coal, 8 gas). `perturb.set_gen_status` excludes them from
headroom redistribution by default.

### The feasibility constraint

Eligible headroom is only **641.5 MW**, and it is concentrated:

| fuel | n | PG (MW) | headroom (MW) |
|---|---|---|---|
| coal | 25 | 1444.8 | 473.8 |
| natural gas | 8 | 71.9 | 167.7 |
| nuclear (slack) | 1 | 569.1 | 0.0 |
| wind (excluded) | 6 | 117.0 | 582.6 |

`135_1` (250.0) and `136_1` (223.8) hold 74 % of eligible headroom; everything
else is at `PMAX`. So **decommit units at `PMAX`** and let those two absorb.
Decommitting *them* collapses absorbing headroom to 167.7 MW and dumps the
shortfall on a slack with no room.

Starting ladder against the 15,158 MW·s fleet:

| Base | Decommitted | Lost MW | Inertia |
|---|---|---|---|
| `c0_full` | none | 0 | 0 |
| `c1_m150` | `126_1` | 150.0 | −1,139 (7.5 %) |
| `c2_m300` | `126_1`, `127_1` | 300.0 | −1,667 (11.0 %) |
| `c3_m450` | `126_1`, `127_1`, `125_1` | 450.0 | −2,009 (13.3 %) |
| `c4_m527` | + `154_1` | 527.2 | −2,791 (18.4 %) |

Adding `153_1`, `152_1`, `155_1` reaches 758.9 MW, exceeds 641.5, and is **not**
feasible. `set_gen_status` warns when redistribution is capped.

Inertia, not headroom, moves the frequency nadir. The ladder is the knob for
making the case sensitive enough that load and wind variation produce a visible
stability response.

Each base directory holds `base.m`, `base_solved.m`, `pf.log` and `meta.json`,
so a campaign can be traced back to exactly which units were offline.

---

## The solver

Start from the shipped `solver.json` and patch it with a dict, same style as
`uq_setup.ipynb`.

### Validation keys must be stripped

The shipped solver is a regression test. It carries `reference_file`,
`error_type`, `error_tolerance` and `abs_err_threshold`, which make
`DynamicSimulation` grade the run against a stored reference trace. **Any
perturbed scenario fails that comparison by construction**, exiting `rc=1`.

`apply_solver_overrides(..., strip_test_keys=True)` removes them, and this is
the default. The failure is loud (`rc=1`) but the cause is not obvious, and it
would otherwise poison an entire campaign.

`output_file` is deliberately **not** stripped: without it `DynamicSimulation`
writes no monitor csv at all.

### Fault impedance

`R`/`X` live on the `BusFault` device in **`case.json`**, not in `solver.json`.
Use `su.set_fault_impedance`. ACTIVSg200 ships `R = 0.0`, `X = 0.01`.

### Validation

`su.validate_solver` catches clearing before onset, `tmax` before the last
event, and `fault_on`/`fault_off` pointing at different elements (a fault that
never clears). `su.describe_solver` prints the bus number and name behind each
`element_id`, so a misplaced fault is visible rather than hidden behind an int.

---

## Quantities of interest

### Observable vs internal

A digital twin is only useful if its outputs can be compared against what the
field reports. A PMU measures the **voltage phasor at a bus**, nothing inside a
machine.

| Quantity | Measurable? |
|---|---|
| bus `\|V\|`, bus angle `θ` | **yes**, PMU |
| bus frequency from `dθ/dt` | **yes**, this is what a PMU calls frequency |
| rotor angle `δ`, rotor speed `ω` | **no**, internal states |

Both groups are computed. The internal ones are the ground truth the terminals
do not show, and the gap between them is itself informative.

### `omega` is a deviation

Verified in `GenrouImpl.hpp`: the state is initialised to `omega(0.0)` (line
465) and enters the algebraic equations as `(ONE + omega)` (line 481). So
synchronous speed is `omega == 0`, absolute pu speed is `1 + omega`, and
frequency is `freq_base · (1 + omega)`. `freq_base` is read from
`case.json` `params`, not assumed; ACTIVSg200 sets 60.0 and
`Component.hpp:370` defaults to the same.

### ROCOF needs a blanking window

Bus frequency is differentiated from the unwrapped angle. The angle **steps** at
fault clearing, and differentiating a step gives a spike that is numerical, not
physical. Without blanking, an early run reported nadir 56.6 Hz and ROCOF
414 Hz/s; with `F_BLANK_S = 0.05` the same run gives 59.91 Hz and 2.1 Hz/s.
These numbers are sensitive to `F_BLANK_S` and `dt_monitor`; treat them as
indicative.

### Loads are 100 % constant impedance

Every `LoadZIP` has `alphaI = alphaP = 0`, hence `alphaZ = 1`. Load current
falls with voltage instead of holding, so **fault-induced delayed voltage
recovery cannot occur in this case by construction**. The `v_recovery_s` metric
reads 0.0 for every scenario run so far, which is structural, not luck.
Modelling FIDVR would need a motor load model this case does not have.

---

## Known limitations

Ordered by when they are likely to bite.

1. **One fault location.** ACTIVSg200 ships a single `BusFault` (bus 149).
   Timing is variable, location is not. Upstream is actively developing fault
   placement; revisit on the dev branch.
2. **No IBR.** See above. Bounds every wind-related conclusion.
3. **Headroom is scarce and concentrated.** 641.5 MW eligible, 74 % of it in two
   units. Commitments must be checked for feasibility before building.
4. **Multi-unit buses.** `parse_dyr` raises and `gen_index` warns if a bus
   carries more than one unit; ids are positional. ACTIVSg200 has 49 generators
   on 49 distinct buses so it never fires, but ACTIVSg2000 and larger will need
   `(bus, model, unit)` keying.
5. **Adding a generator is untested.** `_add_gens_from_dyr` exists but every
   round trip so far reports `gen_added: 0`. Removal is well tested; addition is
   not.
6. **`LoadZ` (Hawaii) not implemented.** Only `LoadZIP` is mapped.
7. **ZIP mix is preserved verbatim.** `alphaI`/`alphaP` are never rewritten, so
   a case with voltage-dependent load keeps it; we only change `Pnom`/`Qnom`.
8. **Synthetic profiles.** `8760_aleatoric` falls back to an invented AR(1) wind
   series and diurnal load shape when `PROFILE_CSV` is absent. It is labelled at
   every output. Nothing computed from it is a result.

---

## Tooling notes

* **GridKit has no `solve_pf` binary** in this build. PowerModels.jl is the only
  power-flow option, via `pm_utils.run_pm_solve_out`. Typical solve: 4 Ipopt
  iterations, ~9 s wall, almost all of it Julia startup.
* **Cost at scale.** 8760 serial runs is roughly a day. Two independent fixes:
  a SLURM array (each hour is self-contained), and amortising Julia startup by
  looping `pm_solve.jl` over a batch in one process.
* **Disk.** `mon.csv` is ~27 MB per run at the shipped `dt_monitor`. 8760 hours
  is ~230 GB. Use `keep_mon=False` to reduce to QoI scalars, or raise
  `dt_monitor`.
* **Python.** `matpowercaseframes` 2.1.0 exports `CaseFrames`, not `Case`, and
  does not parse `mpc.genfuel` (commented out in its `constants.py`), so
  `aux_to_m` parses the fuel cell arrays itself.

---

## Status

### Verified

* `.aux` + `.dyr` → `case.json` reproduces the shipped ACTIVSg200 **exactly**
  (0 diffs) with a 14-parameter overlay on `197_1`
* seed `.m` → `case.json` matches `BASE_CASE` to 1e-6
* commitment changes device counts without moving the fault index
* decommitted base (`c3_m450`, 450 MW) runs end to end: PF `LOCALLY_SOLVED`,
  Vmin 0.9697, 37 Genrou, dynamics `rc=0` with 0 ERROR lines, nadir 59.754 Hz
* wind is excluded from headroom redistribution (wind PG held at 117.00 MW)
* `Vnom` is inert in GridKit, confirmed by source inspection and experiment

### Open

* real hourly profiles (`PROFILE_CSV` does not exist yet)
* the SLURM worker `TODO`: factor the shared setup into `py-utils/campaign.py`
  so array tasks do not duplicate notebook code
* fault placement, pending upstream
* `cases/illinois.md` and `cases/hawaii.md` still describe the pre-#577 shapes
