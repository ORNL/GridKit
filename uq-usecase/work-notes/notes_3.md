# Work Notes — Hawaii v6–v10 sensitivity runs (2026-08-08)

## The H parameter

**H** (inertia constant) is the UQ input parameter varied in all v6–v10 runs.

$$H = \frac{E_k \text{ (stored kinetic energy at rated speed, MJ)}}{\text{machine rating (MVA)}}$$

Units: MJ/MVA = **seconds**. Physical interpretation: if the machine were supplying rated
power with zero mechanical input, its stored kinetic energy would be exhausted in H seconds.
Typical range for synchronous generators: 2–10 s (large steam turbines toward the high end,
small gas turbines toward the low end). KAHE138 units have H = 2.48 s (low end);
WAIPAHU69 units have H = 6.15 s (high end).

H is a per-unit quantity on each machine's own MVA base, not on the system base. The
**physically meaningful measure of a machine's inertia contribution to the grid is
H × mva (MWs)** — see the H×mva section below.

## Key quantities

**Machine MVA (`mva`)**: the generator nameplate rating in MVA; serves as the per-unit
base for H and machine-level power quantities. Each machine has its own `MBASE` in the
case JSON. Values in the Hawaii fleet span a wide range (see fleet table below). Not to be
confused with the system base.

**Grid base MVA (`S_base`)**: the system-wide per-unit base used for all network quantities
(bus voltages, line flows, admittances). The Hawaii case uses S_base = 100 MVA (standard
convention). Voltages in pu and network-level quantities are on this base; H and per-machine
p0 are on each machine's own `mva` base, not this one.

**Generator dispatch (`p0`)**: the active power output at the pre-fault operating point,
in MW (= p0_pu × mva). Determines how much mechanical power the prime mover is delivering
before the fault; a lightly loaded machine (small p0) contributes little to the post-fault
swing even if H×mva is moderate. p0 values for all fleet machines are listed in the fleet
table below.

## Motivation

Previous runs (v1–v5) used a bus-diverse 4-generator selection (one per cluster: ALOHA69,
WAIPAHU69, SCHOFIELD69, KALAELOA138). Analysis downstream raised the question: **why is
H variation for some generators not detectable in post-fault signals?** Three candidate
explanations were identified:

1. **Small H** — the inertia constant itself is low, so the machine stores little energy
   and its speed response is weak relative to noise.
2. **Small machine size** (mva) — a tiny machine contributes negligible total stored energy
   (H×mva) even with a moderate H value; its electrical coupling to the network is also weak.
3. **Electrical distance** — generators far from the fault bus see an attenuated voltage
   perturbation and their response is decoupled from near-fault observables.
4. **Small dispatch** (p0) — a lightly loaded machine contributes little mechanical inertia
   to the post-fault swing; its angle/speed deviation is small relative to others.

Note: in the Hawaii fleet these criteria are **heavily correlated**. Small mva machines
(buses 2 and 36, mva=2.8–3.5 MVA) are simultaneously the smallest H×mva and smallest
dispatch machines. Only electrical distance is independent.

## Sampling range rationale

All v6–v9 use **uniform ±20%** of nominal H. The choice of ±20% reflects:

| range | interpretation |
|-------|---------------|
| ±2–5% | Manufacturer nameplate: rotor mass/geometry are well-measured |
| ±10–15% | NERC/IEEE standard epistemic assumption for utility database values |
| ±20% | "Fairly uncertain database value" — defensible calibration prior |
| ±50%+ | Stress-test / identifiability bound; not a realistic calibration prior |

v5 used normal std=12% (≈ one-sigma covers ±12%), which sits in the "database" tier.
±20% uniform is wider and flat — it makes no assumption about where in the range the
true value sits, consistent with treating H as unknown-but-bounded.

## H×mva: why it matters

`H` is in per-unit on each machine's own MVA base (`MBASE`), not on the system base
(100 MVA). Each machine has its own `mva` that varies from 2.8 MVA (ALOHA69 units) to
124.3 MVA (KALAELOA138 `35_8`). Total stored kinetic energy = H × mva (MWs). Two
machines with the same H contribute very different inertia to the grid:

| id | bus | H | mva (MVA) | H×mva (MWs) |
|----|-----|---|-----------|-------------|
| `2_1` | ALOHA69 | 3.69 | 2.8 | 10.3 |
| `23_1` | WAIPAHU69 | 6.15 | 85.6 | 526.4 |

Smallest H ≠ smallest inertia. The KAHE138 machines (`37_3`, `37_5`) have the lowest H
in the fleet (2.48 s) but are large (96–86 MVA) so H×mva ≈ 238–213 MWs — much more
inertia than the tiny ALOHA69 machines.

## Fleet overview (all 39 Genrou, key quantities)

One representative per bus (smallest H×mva within bus), sorted by H×mva:

| id | bus | bus_name | H | mva (MVA) | H×mva (MWs) | p0 (MW) | hop from fault |
|----|-----|----------|---|-----------|-------------|---------|----------------|
| `2_1` | 2 | ALOHA69 | 3.69 | 2.8 | 10.3 | 0.07 | 1 (transformer) |
| `36_1` | 36 | COGEN69 | 4.42 | 3.5 | 15.5 | 0.07 | 3+ |
| `26_2` | 26 | EWA BEACH69 | 3.00 | 11.2 | 33.6 | 1.14 | 2 |
| `34_1` | 34 | SCHOFIELD69 | 4.35 | 9.2 | 40.0 | 0.66 | 3+ |
| `23_9` | 23 | WAIPAHU69 | 3.00 | 16.2 | 48.6 | 2.03 | 2 |
| `26_1` | 26 | EWA BEACH69 | 3.00 | 22.0 | 66.0 | 4.40 | 2 |
| `27_2` | 27 | KAHUKU69 | 3.00 | 30.4 | 91.2 | 8.39 | far north |
| `33_1` | 33 | WAIANAE69 | 3.00 | 30.4 | 91.2 | 8.39 | 3+ west |
| `27_1` | 27 | KAHUKU69 | 3.00 | 33.0 | 99.0 | 9.90 | far north |
| `35_4` | 35 | KALAELOA138 | 5.22 | 22.0 | 114.8 | 4.40 | 1 (direct) |
| `37_5` | 37 | KAHE138 | 2.48 | 85.9 | 213.0 | 58.84 | 2 |
| `37_3` | 37 | KAHE138 | 2.48 | 95.9 | 237.8 | 55.94 | 2 |
| `28_2` | 28 | HALEIWA69 | 3.00 | 75.9 | 227.7 | 52.40 | far north |
| `28_1` | 28 | HALEIWA69 | 3.00 | 53.9 | 161.7 | 24.98 | far north |
| `23_1` | 23 | WAIPAHU69 | 6.15 | 85.6 | 526.4 | 59.30 | 2 |
| `35_8` | 35 | KALAELOA138 | 5.22 | 124.3 | 648.9 | 69.90 | 1 (direct) |

Fault is at bus 1 (ALOHA138). p0 is computed as p0_pu × mva (MW dispatch).

## Experiment designs

All runs (v6–v10): N=4000, LHS, uniform ±20%, `MONITORS_BY_CLASS = {bus: [Vm, Va], genrou: [delta, omega]}`.
Selection principle: one generator per bus for diverse dynamics parameters.

---

### v6 — smallest H (inertia constant)

**Hypothesis**: generators with the lowest H have the weakest speed/angle response to a
fault; their contribution to observable signals is intrinsically small regardless of size.

Selection: lowest H value per bus, one per bus, requiring distinct dynamics parameter sets.

| id | bus | bus_name | H (nominal) | H range [lo, hi] | mva (MVA) | H×mva (MWs) | p0 (MW) | hop |
|----|-----|----------|-------------|------------------|-----------|-------------|---------|-----|
| `37_3` | 37 | KAHE138 | 2.48 | [1.984, 2.976] | 95.9 | 237.8 | 55.94 | 2 |
| `26_2` | 26 | EWA BEACH69 | 3.00 | [2.400, 3.600] | 11.2 | 33.6 | 1.14 | 2 |
| `27_1` | 27 | KAHUKU69 | 3.00 | [2.400, 3.600] | 33.0 | 99.0 | 9.90 | far north |
| `33_1` | 33 | WAIANAE69 | 3.00 | [2.400, 3.600] | 30.4 | 91.2 | 8.39 | 3+ west |

Note: `37_3` has the lowest H in the fleet (2.48 s) but is a large machine (H×mva=238 MWs).
The three H=3.0 machines are chosen for geographic diversity (south/north/west periphery)
rather than being the globally smallest H×mva.

---

### v7 — farthest from fault (electrical periphery)

**Hypothesis**: generators electrically remote from the fault bus see a weak voltage
perturbation; their post-fault dynamics are largely decoupled from near-fault signals.

Selection: buses on the geographic and electrical periphery of the Oahu grid.

| id | bus | bus_name | H (nominal) | H range [lo, hi] | mva (MVA) | H×mva (MWs) | p0 (MW) | hop |
|----|-----|----------|-------------|------------------|-----------|-------------|---------|-----|
| `27_1` | 27 | KAHUKU69 | 3.00 | [2.400, 3.600] | 33.0 | 99.0 | 9.90 | far north shore |
| `28_1` | 28 | HALEIWA69 | 3.00 | [2.400, 3.600] | 53.9 | 161.7 | 24.98 | far north shore |
| `33_1` | 33 | WAIANAE69 | 3.00 | [2.400, 3.600] | 30.4 | 91.2 | 8.39 | west, hop 3+ |
| `34_1` | 34 | SCHOFIELD69 | 4.35 | [3.480, 5.220] | 9.2 | 40.0 | 0.66 | central, hop 3+ |

Note: v6 and v7 share `27_1` and `33_1` — these machines are both far from the fault
and have low H. This overlap is intentional: if both runs show low sensitivity for these
two machines, that is consistent with both hypotheses; the non-overlapping machines
help disambiguate (e.g. `37_3` in v6 is close to the fault, so if it also shows low
sensitivity, distance is not the explanation).

---

### v8 — smallest mva (machine rating)

**Hypothesis**: machines with the smallest MVA rating are electrically weakly coupled
to the network; their speed/angle perturbations contribute negligibly to bus voltages.
Also tests the confounded "small dispatch" hypothesis since small mva ↔ small p0.

Selection: smallest mva across all buses, one per bus, requiring distinct dynamics.

| id | bus | bus_name | H (nominal) | H range [lo, hi] | mva (MVA) | H×mva (MWs) | p0 (MW) | hop |
|----|-----|----------|-------------|------------------|-----------|-------------|---------|-----|
| `2_1` | 2 | ALOHA69 | 3.69 | [2.952, 4.428] | 2.8 | 10.3 | 0.07 | 1 (transformer) |
| `36_1` | 36 | COGEN69 | 4.42 | [3.536, 5.304] | 3.5 | 15.5 | 0.07 | 3+ |
| `34_1` | 34 | SCHOFIELD69 | 4.35 | [3.480, 5.220] | 9.2 | 40.0 | 0.66 | 3+ |
| `26_2` | 26 | EWA BEACH69 | 3.00 | [2.400, 3.600] | 11.2 | 33.6 | 1.14 | 2 |

Note: `2_1` is at bus 2 (ALOHA69) which is electrically very close to the fault
(hop 1 through the transformer), but extremely small (mva=2.8, p0=0.07 MW). If even
this machine shows no sensitivity despite proximity, size/dispatch is likely the explanation.

---

### v9 — largest H×mva (positive control)

**Hypothesis (null)**: machines with the largest total stored inertia should show the
strongest sensitivity in post-fault signals. If H variation is undetectable even here,
the issue lies elsewhere — in the choice of observable, metric, or analysis method.

Selection: largest H×mva across the fleet, one per bus, geographically diverse.

| id | bus | bus_name | H (nominal) | H range [lo, hi] | mva (MVA) | H×mva (MWs) | p0 (MW) | hop |
|----|-----|----------|-------------|------------------|-----------|-------------|---------|-----|
| `35_8` | 35 | KALAELOA138 | 5.22 | [4.176, 6.264] | 124.3 | 648.9 | 69.90 | 1 (direct branch) |
| `23_1` | 23 | WAIPAHU69 | 6.15 | [4.920, 7.380] | 85.6 | 526.4 | 59.30 | 2 |
| `37_3` | 37 | KAHE138 | 2.48 | [1.984, 2.976] | 95.9 | 237.8 | 55.94 | 2 |
| `28_2` | 28 | HALEIWA69 | 3.00 | [2.400, 3.600] | 75.9 | 227.7 | 52.40 | far north |

Note: `37_3` appears in both v6 and v9. In v6 it is the lowest-H machine; in v9 it is
the largest-H×mva machine at its bus (KAHE138). This overlap is useful: if this machine
shows sensitivity in v9 but not in v6, the H×mva magnitude matters more than H alone.

---

### v10 — small p0 (dispatch) with moderate H×mva

**Hypothesis**: even with moderate H×mva, a lightly dispatched machine contributes
negligible mechanical energy to the post-fault swing; its H variation is undetectable
because p0 ≈ 0 means the rotor is effectively coasting rather than driving the network.

**Caveat**: in this operating snapshot, small p0 and small H×mva are heavily correlated
across the fleet. A clean p0 isolation is not possible within a single snapshot. The best
available approach is to select machines that are lightly dispatched but at the *same
buses* as v9 machines, so electrical location is held constant and only p0 and mva differ.
A true p0 isolation requires different dispatch snapshots (future work).

**Key same-bus pairs vs v9** (same bus, different p0):

| v10 machine | v9 machine | bus | v10 H×mva (MWs) | v9 H×mva (MWs) | v10 p0 (MW) | v9 p0 (MW) |
|-------------|------------|-----|-----------------|----------------|-------------|------------|
| `35_4` | `35_8` | KALAELOA138 | 114.8 | 648.9 | 4.40 | 69.90 |
| `23_9` | `23_1` | WAIPAHU69 | 48.6 | 526.4 | 2.03 | 59.30 |

For the bus-35 pair, H is identical (5.22 s) in both machines — only mva and p0 differ.
This is the cleanest within-snapshot comparison the fleet allows.

Selection for v10: lightly dispatched machines at electrically close buses, one per bus.

| id | bus | bus_name | H (nominal) | H range [lo, hi] | mva (MVA) | H×mva (MWs) | p0 (MW) | hop |
|----|-----|----------|-------------|------------------|-----------|-------------|---------|-----|
| `35_4` | 35 | KALAELOA138 | 5.22 | [4.176, 6.264] | 22.0 | 114.8 | 4.40 | 1 (direct) |
| `23_9` | 23 | WAIPAHU69 | 3.00 | [2.400, 3.600] | 16.2 | 48.6 | 2.03 | 2 |
| `26_1` | 26 | EWA BEACH69 | 3.00 | [2.400, 3.600] | 22.0 | 66.0 | 4.40 | 2 |
| `27_2` | 27 | KAHUKU69 | 3.00 | [2.400, 3.600] | 30.4 | 91.2 | 8.39 | far north |

**Interpretation logic** (contingent on v9 results):

- v9 sensitive, v10 not: consistent with p0/dispatch being the limiting factor.
- v9 not sensitive: skip v10 or defer; the issue lies elsewhere.
- v10 also sensitive: p0 is not a barrier at these levels; H×mva may be the key quantity.

---

## Summary comparison table

| run | hypothesis | gens | N | pct | common theme |
|-----|-----------|------|---|-----|--------------|
| v6 | smallest H | `37_3`, `26_2`, `27_1`, `33_1` | 4K | ±20% | low H constant |
| v7 | farthest from fault | `27_1`, `28_1`, `33_1`, `34_1` | 4K | ±20% | electrical periphery |
| v8 | smallest mva + dispatch | `2_1`, `36_1`, `34_1`, `26_2` | 4K | ±20% | tiny machines |
| v9 | largest H×mva (positive control) | `35_8`, `23_1`, `37_3`, `28_2` | 4K | ±20% | high inertia contribution |
| v10 | small p0 / lightly dispatched | `35_4`, `23_9`, `26_1`, `27_2` | 4K | ±20% | low dispatch, moderate H×mva |

Overlapping generators across runs (useful for cross-run comparison).
Note: v10 shares no exact generator IDs with v6–v9, but uses machines at the same buses
as v9 (`35_4` vs `35_8`, `23_9` vs `23_1`) — the key same-bus comparison.

| gen | v6 | v7 | v8 | v9 | v10 |
|-----|----|----|----|----|-----|
| `27_1` | ✓ | ✓ | | | |
| `33_1` | ✓ | ✓ | | | |
| `34_1` | | ✓ | ✓ | | |
| `26_2` | ✓ | | ✓ | | |
| `37_3` | ✓ | | | ✓ | |

## Run sequence

One run at a time: select the run root (Step 1) and uncomment the matching PARAM_SPECS
block (Step 2) in the hawaii setup cell of `uq_setup.ipynb`, then run meta → samples →
run dirs → SLURM submit. Wait for completion before switching to the next version.

v10 should be run after reviewing v9 results. If v9 shows no H sensitivity at all, v10
can be skipped or deferred — the low-dispatch hypothesis is moot if the positive control fails.

Data lands in `/kfs2/projects/scidac/scidac-data/gridkit-runs/hawaii-v{6,7,8,9,10}/`.

### hawaii-ian-csv (Morris sensitivity design)

An externally provided Morris sensitivity design (`hawaii_morris_sensitivity_analysis.csv`)
covering all 39 Genrou devices in the Hawaii fleet was run using the **hawaii setup from CSV**
cell in `uq_setup.ipynb`. N=1000 samples (from the CSV), all 39 H parameters varied
simultaneously. Data lands in `/kfs2/projects/scidac/scidac-data/gridkit-runs/hawaii-ian-csv/`.
This run is independent of v6–v10 and was not designed by this session.
