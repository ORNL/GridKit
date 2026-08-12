# m_to_case helper: MATPOWER `.m` → GridKit `case.json`

Companion doc for [`notebooks/m_to_case_helper.ipynb`](../notebooks/m_to_case_helper.ipynb).

Purpose: document and rigorously test the `.m` → `case.json` mapping in isolation,
before any downstream workflow (aleatoric UQ, dynamic simulation) relies on it.

## Scope

- **In scope**: field-by-field mapping from a MATPOWER PF solution to a GridKit
  `case.json`; supported perturbation types (load, wind, gen-off); tolerance
  philosophy for the identity round-trip.
- **Out of scope**: producing the solved `.m` in the first place (that belongs
  to [`pm_helper.md`](pm_helper.md) / PM.jl), running the resulting `case.json`
  through `DynamicSimulation` (that belongs to `aleatoric_helper.md`, to come),
  full device-list rebuild for adding NEW generators (deferred to production
  `aleatoric_setup.ipynb`).

## Two entry points (they coexist)

| Function | Use when | Modifies |
|---|---|---|
| `patch_case_from_m()` | Only `PG`, `QG`, `Vm`, `Va` changed (same commitment, same load structure) | bus `init.Vr/Vi`, Genrou `p0/q0` |
| `build_case_from_solved_m()` | Load, wind, or generator commitment changed | Same as above **plus** `LoadZIP.{Pnom,Qnom,Vnom}`, plus removal of `{Genrou, Tgov1, SexsPti/Ieeet1}` triplets for offline gens |

`patch_case_from_m` is unchanged and remains the default for the epistemic track
(`uq_setup.ipynb`, illinois-v2), which never varies load or commitment. The new
`build_case_from_solved_m` targets the aleatoric track.

## The pipeline (context)

```
base .m                                           ← case_ACTIVSg200.m
   │
   ├─ perturb  (make_load_scale_m / make_wind_dispatch_m / make_gen_off_m)
   │
   └─→ perturbed .m
       │
       ├─ PM.jl solve PF  (run_pm_solve_out)     ← pm_helper.md
       │
       └─→ solved .m                              ← same MATPOWER format, updated
           │                                        VM/VA/PG/QG (and possibly
           │                                        GEN_STATUS)
           │
           ├─ build_case_from_solved_m(base_json, solved_m, out.json)  ← THIS DOC
           │
           └─→ case.json
               │
               └─ DynamicSimulation  (aleatoric_helper.md, to come)
```

Everything before the fourth stage is already covered in `pf_helper.md`,
`pm_helper.md`, and `pf_utils.py`. This doc concentrates on the fourth stage
in isolation.

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
Vnom_pu = VM
```

Then split the LoadZIP devices at each bus into two groups:

- **Real loads**: `Pnom > 0` OR (`Pnom == 0` AND `Qnom >= 0`), AND `id` does
  not start with `"shunt_"`.
- **Shunts / capacitors**: `Pnom == 0` AND `Qnom < 0`, OR `id.startswith("shunt_")`.

Shunt-like devices are excluded from PD/QD patching (they represent reactive
compensation, not real load). Their `Vnom` **is** refreshed to `VM` so their
per-unit reference tracks the new bus voltage.

For the real loads at each bus:

- **N=1**: straightforward — set `Pnom=Pnom_pu`, `Qnom=Qnom_pu`, `Vnom=VM`.
- **N≥2 (multi-ZIP)**: proportional scaling to preserve each device's share of
  the bus load:

  ```
  S_P = sum(Pnom_old for devs)      # existing sum in base JSON
  S_Q = sum(Qnom_old for devs)
  for each dev:
      Pnom_new = Pnom_old * (Pnom_pu / S_P)
      Qnom_new = Qnom_old * (Qnom_pu / S_Q)
      Vnom_new = VM
  ```

  This preserves the ratio between multiple loads at the same bus (Illinois has
  38 buses with 2 LoadZIP devices each; see illinois.md).

- **Zero-load bus (PD=QD=0)**: set all real LoadZIP `Pnom=0, Qnom=0`, `Vnom=VM`.
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
- Numerical: `bus.init.Vr/Vi` and `LoadZIP.{Pnom,Qnom,Vnom}` and
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
- Genrous updated: exactly `n_online_gen_in_JSON` (40 for illinois-v2 base).
- Genrous removed: 0 (nothing is offline in the base `.m` that's online in
  the JSON).
- LoadZIP updated: exactly the 108-load-bus × device count from illinois.md.
- Shunts skipped: exactly 1 (bus 15).
- Any structural change in `out.json` vs `illinois.json` **other than** the
  four field families above → hard fail.

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
- `Vnom` updated to `VM`.
- `alphaI`, `alphaP` unchanged.

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
- Its `Vnom` **is** updated to bus 15's `VM` in the `.m`.

### 4g. Zero-load bus

Force `PD=QD=0` on a bus that has one or more LoadZIPs. Verify:

- All real LoadZIPs at that bus have `Pnom=0, Qnom=0`, `Vnom=VM`.

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
| `mpc.bus[:, 8]` (`VM`) | 8 | pu | Bus init `Vr, Vi`; LoadZIP `Vnom` |
| `mpc.bus[:, 9]` (`VA`) | 9 | deg | Bus init `Vr, Vi` |
| `mpc.gen[:, 1]` (`GEN_BUS`) | 1 | int | Genrou id decode |
| `mpc.gen[:, 2]` (`PG`) | 2 | MW | Genrou `p0` |
| `mpc.gen[:, 3]` (`QG`) | 3 | MVAr | Genrou `q0` |
| `mpc.gen[:, 8]` (`GEN_STATUS`) | 8 | 0/1 | offline removal |
| `mpc.baseMVA` | header | MVA | per-unit conversion |

## Known limitations (to lift in production)

1. **Device-list rebuild for new-online gens** (Section 3): if a gen is online
   in the `.m` but absent from the base JSON, `build_case_from_solved_m` raises
   in `strict` mode. Production `aleatoric_setup.ipynb` will need a companion
   `add_genrou_from_m()` helper that synthesises reasonable Genrou / Tgov1 /
   SexsPti defaults from the `.m` gen row (`Pmax`, `Qmax`, `Vg`, etc.), and
   from a fuel-type-based lookup table for the dynamic parameters (`H`, `D`,
   `Xd`, ...).
2. **`LoadZ` support** (constant-impedance load model). Illinois uses only
   `LoadZIP`; Hawaii uses a simpler `Load`. If a case uses `LoadZ`, the
   patch algorithm changes: `G = PD_pu / VM^2`, `B = -QD_pu / VM^2`, `R =
   G / (G^2 + B^2)`, `X = -B / (G^2 + B^2)`. Not implemented until a case
   requires it.
3. **Multi-ZIP with mixed real+shunt at the same bus**: current code
   correctly identifies shunts and skips them from PD/QD scaling, but
   assumes at most one shunt per bus. Any bus with multiple shunt-like
   LoadZIPs would need reactive splitting. Not observed in Illinois.
4. **`alphaI`/`alphaP` sensitivity**: the ZIP mix is preserved verbatim. If
   a scenario needs a different ZIP mix (e.g. varying constant-power
   fraction with load level), that is out of scope for the mapping helper
   and belongs in a separate device-model editor.

## Findings

Filled in after running [`m_to_case_helper.ipynb`](../notebooks/m_to_case_helper.ipynb).

- **Identity round-trip Vr/Vi residual**: [TODO]
- **Identity round-trip LoadZIP change vs base**: [TODO]
- **Genrou removal test on illinois.json**: [TODO]
- **Multi-ZIP scaling exact-match**: [TODO]
- **End-to-end 0.8x load perturbation Pnom sum**: [TODO]

## Forward pointers

- `aleatoric_helper.md` (to come): consumes `build_case_from_solved_m` end-to-end
  (perturb → PM.jl → build case.json → DynamicSimulation → QoI compare).
- `aleatoric_setup.md` / `aleatoric_setup.ipynb` (to come, production):
  scales the workflow to N hours of user-provided load/wind timeseries; adds
  device-list rebuild and kill-list builder.
