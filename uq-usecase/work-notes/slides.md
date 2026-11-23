# Datacenter Load and Transient Security: Pilot Study Slides

---

## Highlight Slide: Overview

**Title:** Datacenter Load and Transient Security: An Aleatoric UQ Pilot on GridKit

**Central Question:**
Where, and how much, large datacenter load can be added to the Illinois grid before dynamic security degrades, given that IBR output makes system inertia vary hour to hour?

**Study Scope:**
- **Grid:** ACTIVSg200 (Illinois synthetic case, 200 buses)
- **Method:** Sampling-based aleatoric uncertainty quantification
- **Key Metric:** Critical Clearing Time (CCT) — longest fault system survives before generator loss of synchronism
- **Finding:** Constant-power datacenter load shifts CCT distribution significantly; 100 MW placement raises failure risk from ~50% to ~73%

**Outcomes:**
- Demonstrated GridKit as end-user tool on new application
- Identified cost drivers across sampling dimensions (operating conditions, faults, placements)
- Profiled scalability bottlenecks: sampling-wise analysis, fault screening, outcome categorization for 2000-bus deployment

---

## Slide 1: The Research Question

**Title:** Datacenter Load Impact on Stability: A UQ Sampling Study

**Central Question:**
How much large constant-power datacenter load can be safely added to a meshed grid with high IBR penetration, and how does placement location affect critical clearing time?

**Why It Matters:**
- High IBR penetration reduces system inertia variability across time and operating conditions
- Constant-power load representation (datacenter) differs fundamentally from traditional voltage-dependent loads
- Placement location, load size, and operating condition interact; quantifying dependencies requires sampling

**Study Approach:**
- Pilot on 200-bus Illinois synthetic grid with real hourly load/wind data
- Apply sampling-based aleatoric UQ to quantify stability margins
- Establish and validate sampling methodology for 2000-bus deployment

**[Figure: Illustration of grid topology with IBR + datacenter siting locations]**

---

## Slide 2: Glossary — Key Terms

| Term | Meaning |
|---|---|
| **Operating condition** | One (hour, commitment) pair. Hour fixes load and wind; commitment fixes which thermal units are online. |
| **Baseline** | Same operating condition with no datacenter. Change in CCT is always measured on the same condition. |
| **Net load** | Demand minus wind, system-wide. What the thermal fleet must serve. |
| **Critical Clearing Time (CCT)** | Longest fault duration system survives before generator loses synchronism. Larger is safer. |
| **Censored** | System survived the longest fault tested (0.5 s cap). We know CCT ≥ 0.5 s but not its value. |
| **Infeasible dispatch** | Committed generators cannot serve the load — finding about the system, not an error. |
| **ZIP load** | Load model splitting demand into constant-impedance (Z), constant-current (I), and constant-power (P) parts. Datacenter uses 80/20/0. |
| **Placement** | Where the datacenter sits: none (baseline), or one of buses 6, 108, 120, 148. |

---

## Slide 3: The Illinois Synthetic Grid (ACTIVSg200)

**Grid Characteristics:**
- **200 buses, 246 branches, 49 generators** (43 thermal: 25 coal, 17 gas, 1 nuclear; 6 wind IBR)
- **40 synchronous machines** (34 thermal with dynamic models online, 6 wind units, 9 offline thermal units can be committed)
- **109 loads** across the network
- Hourly load and wind are synthetic one-year series (TAMU repository)

**One-Time Seed Building:**

```
TAMU .aux (topology, dispatch) ──┐
                                  ├──> rebuild case.json ──> BASE_CASE
TAMU .dyr (dynamics) ────────────┘    (14 diffs, overlay applied)
                                  ↓
                           seed .m (dispatch)
                                  ↓
                        per condition:
                   build_case_from_m(BASE_CASE + seed)
                                  ↓
                           case.json (per condition)
```

**Base case equals shipped case after validation; dispatch overwritten per condition; six wind units become negative load stand-ins (no inertia).**

**[Figure 3.1: Seed building diagram showing blue one-time build, grey per-condition build]**

---

## Slide 4: Building One Operating Condition

**Two Decisions: Commitment and Dispatch**

### 4.1 Commitment: A Feasible Band

**Not an economic optimum.** Thermal units ranked by cost (nuclear always on). Committing *n* units switches on the first *n*.

**Feasible band computed per hour:**
- **Lower bound:** enough capacity to cover net load + 10% reserve, minimum 8 units (fewer → voltage instability)
- **Upper bound:** committed units' minimum output ≤ net load (else infeasible)
- **Example** (762 MW net load): band is 18 to 38 units → **1.56× inertia difference** between ends

**Sampling strategy:** Per hour, commit at two fixed positions: 25% and 75% from band lower to upper end.

**Why commitment matters:** At fixed load and wind, commitment sets inertia and headroom. Results reported **conditional on commitment**.

### 4.2 Dispatch: DCOPF → ACPF

```
DCOPF (lossless) ──> ACPF (solve V, Q, losses) ──> loss correction
                           │
                           ├──> PV→PQ switching (enforce Q limits)
                           │
                           └──> next DCOPF iteration
```

**Key features:**
- **One loss correction** (not two): 27.05 MW overshoot → 0.34 MW residual on worst hour
- **PV→PQ switching:** generators hitting Q limits fixed at limit, bus switched to PQ, re-solved
- **Q limits enforcement:** default ACPF left 2–11 generators outside limits; switching brings all within

---

## Slide 5: Operating Condition Dynamics

### 4.3 Wind Becomes Negative Load

**Wind units removed as synchronous machines.** IBR (inverter-based resources) provide no synchronous inertia—real wind is a **controllable current source** (WECC/Sandia model).

**Implementation:** Hourly wind MW pinned in DCOPF; ACPF solves reactive output Q. Write (P, Q) as negative-load ZIP with **constant-current model.**

**Result:** Wind injection scales with |V|, correct for constant-current control. At initial voltage reproduces power-flow solution exactly.

### 4.4 The Datacenter

**Representation:** Additional electrical load at one bus, **added before dispatch** (OPF sees it; commitment and dispatch produce and deliver it).

**ZIP Load Model (80/20/0):**
- **80% constant-power (P):** UPS + PFC in IT racks maintain constant draw across voltage sag
- **20% constant-current (I):** Motor-based cooling (chillers, pumps, fans) with voltage-dependent draw
- **0% constant-impedance (Z):** Pure datacenter IT portion

**Why this matters:** Constant-impedance load sheds itself during faults and helps recovery. Constant-power load keeps drawing, depressing voltage further → destabilizing.

**Placement:** Campaign faultscreen-2 at bus 108 (100 MW). Campaigns aleatoric-1 and aleatoric-2 at buses 6, 108, 120, 148 (50–400 MW).

---

## Slide 6: Critical Clearing Time by Bisection

**Stability Criterion:**

Pole slip relative to **Centre of Inertia (COI)** — the inertia-weighted average rotor angle:

$$\delta_{COI}(t) = \frac{\sum_i H_i S_i \delta_i(t)}{\sum_i H_i S_i}$$

**Run is unstable when:** $\max_{i,t} \bigl| \delta_i(t) - \delta_{COI}(t) \bigr| > \pi$

(Any machine leaves the collective inertial frame → loss of synchronism)

**CCT Search (Bisection):**
- Fault applied at $t_{on} = 1.0$ s, cleared at $t_{off} = 1.0 + \text{duration}$
- Search interval: [0.05 s, 0.5 s] (protection clears faults at ~0.05–0.5 s)
- Tolerance: 0.01 s
- **Censored run:** stable at top (0.5 s) → reported as CCT ≥ 0.5 s
- **Typical run:** boundary inside [0.05, 0.5] → 8 simulations via bisection

**Efficiency:** Recording only 21 of 441 output columns (rotor angles) → **13.8× storage reduction**

**[Figure 5.1: Rotor-angle deviation from COI showing unstable (red trace crosses π) vs censored (peaks at 0.47π) samples]**

---

## Slide 7: Sampling Space and Three Campaigns

**Three independent campaigns, all built the same way (§4), varying different dimensions:**

| Dimension | **aleatoric-1** | **aleatoric-2** | **faultscreen-2** |
|---|---|---|---|
| **Purpose** | How much operating conditions move CCT | Replication & convergence | Which faults separate conditions |
| **Hours** | 20 (4 per quintile) | 60 (12 per quintile) | 50 (10 per quintile) |
| **Conditions** | 40 | 180 | 100 |
| **Faults** | 1 (bus 66) | 1 (bus 66) | 17 locations |
| **Datacenter buses** | 6, 108, 120, 148 | 6, 108, 120, 148 | 108 only |
| **Datacenter sizes** | 50–400 MW | 50–400 MW | 100 MW |
| **Total runs** | 840 | 3,780 | 3,400 |
| **Runs with CCT** | 439 (52%) | 2,155 (57%) | 3,170 (93%) |

**Stratification:** Hours drawn randomly within net-load quintiles (8,760 hours split into 5 equal-count bins). Each quintile represents 20% of year's hours.

---

## Slide 8: Reading the Curves — Kaplan–Meier and Censoring

**Kaplan–Meier Estimator:**

Correctly handles **censored** runs (stable at 0.5 s cap, true CCT unknown):

$$\hat S(t)=\prod_{t_i\le t}\left(1-\frac{d_i}{n_i}\right), \quad P(\text{CCT}<t)=1-\hat S(t)$$

Where:
- $t_i$ = distinct measured CCT values
- $d_i$ = number unstable at $t_i$
- $n_i$ = at-risk count (measured at or after $t_i$, plus censored)

**Censored runs kept in denominator** (not dropped) → valid inference

**Greenwood Variance:** Estimates SE of KM estimate:

$$\widehat{\text{Var}}\big(\hat S(t)\big)=\hat S(t)^2\sum_{t_i\le t}\frac{d_i}{n_i(n_i-d_i)}$$

95% CI: $P(\text{CCT}<t) \pm 1.96 \cdot \text{SE}$

### 7.2 Simple Example: 200 MW at Bus 108, Fault at Bus 66

**Operating conditions:** Matched pair of 37 conditions, aleatoric-1

| | No datacenter | 200 MW at bus 108 |
|---|---|---|
| Conditions | 37 | 37 |
| Censored at 0.5 s | 46% | 22% |
| **P(CCT < 0.25 s)** | **0.49** | **0.73** |
| 95% Greenwood CI | 0.325–0.648 | 0.587–0.873 |
| Median CCT | 0.27 s | 0.18 s |

**Takeaway:** Datacenter raises failure risk by ~0.24 (95% CI overlap at n=37; paired test at n=171 shows significant shift).

**[Figure 7.2: KM P(CCT > t) baseline (blue) vs 200 MW at bus 108 (orange), 37 conditions, fault at bus 66]**

---

## Slide 9: Risk by Net Load

**Data:** aleatoric-1, all datacenter runs (all buses, all sizes), 1 fault location (bus 66)

| Net-load quintile | Mean load (MW) | Runs solved | Censored | **P(CCT < 0.25 s)** |
|---|---|---|---|---|
| q0 (lowest) | 809 | 75 of 160 | 83% | **0.07** |
| q1 | 1,001 | 82 of 160 | 72% | **0.16** |
| q2 | 1,221 | 82 of 160 | 0% | **0.73** |
| q3 | 1,477 | 80 of 160 | 0% | **1.00** |
| q4 (highest) | 1,660 | 80 of 160 | 0% | **1.00** |

**Key observation:** Sharpest step between q1 and q2 (1,001 → 1,221 MW net load).

**Why:** At low net load, committed inertia is high (many units at part load) → censored at 0.5 s. At high net load, fewer units carry more load → lower inertia, faster CCT.

**[Figure 7.5: KM P(CCT > t) by net-load quintile, all datacenter configurations]**

**Takeaway:** Risk concentrates in high-net-load quintiles. This is where to focus. Results conditional on runs that solved.

---

## Slide 10: Fault Screening — Which Faults Separate Operating Conditions

**Question:** Before running aleatoric campaign, which 17 candidate fault locations can tell operating conditions apart?

**Sampling:** faultscreen-2 uses 50 hours × 2 commitments = 100 conditions, each run at all 17 faults, with and without 100 MW datacenter at bus 108.

**Method:** Plot share of runs losing synchronism (before 0.5 s cap) vs net-load quintile. Separating faults show rising curve; flat line = no separation.

**Results:**
- **4 faults separate conditions:** buses 66, 68, 94, 181 (outcome changes with net load)
- **13 faults are flat:** censored in every quintile, cannot separate

**Implication for 2000-bus grid:** With 10× buses, contingency list grows ~115×. Screening cost grows ~19–25× per run. Need cheap pre-screen to identify discriminating faults before full campaign.

**[Figure 7.6: Fraction of runs losing synchronism vs quintile, 17 fault buses. Four non-flat lines (66, 68, 94, 181), 13 flat lines]**

**Takeaway:** Fault location is a major source of variance. Graph-based surrogate on these runs could identify active faults adaptively, reducing 2000-bus sampling cost.

---

## Slide 11: Scaling to Texas — ACTIVSg2000

**Small scoping test on ACTIVSg2000** (2,000 buses, TAMU repository). Question: how does cost grow?

| Quantity | ACTIVSg200 | ACTIVSg2000 | Ratio | Scaling |
|---|---|---|---|---|
| Buses | 200 | 2,000 | 10× | — |
| Branches | 246 | 3,206 | 13× | — |
| **DCOPF (warm)** | 0.011 s | 0.212 s | 19.1× | $N^{1.27}$ |
| **ACPF (warm)** | 0.018 s | 0.331 s | 18.5× | $N^{1.27}$ |
| **One dynamic run** (full monitors) | 0.28 s | 7.08 s | 25.3× | $N^{1.40}$ |
| Monitor output (raw CSV) | 27.3 MB | 262.5 MB | 9.6× | $N^{0.98}$ |

**Projected aleatoric-1 campaign on Texas (840 runs, ~7,560 simulations):**

| Resource | ACTIVSg200 | ACTIVSg2000 |
|---|---|---|
| Serial dynamic time | ~35 min | ~14.9 h |
| Wall time, 16 workers, 1 node | ~2.2 min (ideal) | ~56 min (ideal) |
| Actual (measured 2.7× slower) | 6.0 min | **~2.5 h** |
| Monitor storage (raw) | 22 GB | **~220 GB** |

**Takeaway:** Brute force survives one order of magnitude of grid size, not two. At 2000-bus scale:
- Design space grows 115× (siting 11.5× × faults 10×)
- Cost per run grows 19–25×
- Serial time grows 2,150–2,910×
- Storage becomes a bottleneck

**Solution:** Screening + surrogate models (e.g., graph NN trained on existing pilot runs) + adaptive sampling for next phase.

---

## Slide 12: Next Steps — FY27 Tracks

**Track 1: Spend fewer, better-chosen dynamic runs**

1. **Screen contingencies cheaply.** faultscreen-2 showed 13 of 17 faults are censored everywhere (do not separate). Texas grid: contingency list grows with grid size; screening cuts cost.

2. **Reduce CCT search cost per run.** Censored runs need 2 simulations. Measured runs need bisection. Options:
   - Start bisection from neighbouring dispatch (strongly correlated CCT)
   - Early termination detection (pole slip detectable before simulated window end)

3. **Stratified/importance sampling.** Variance decomposition shows fault location explains >50% variance. Stratify on top-ranked factors (fault location, net load, placement) rather than uniform sampling.

**Track 2: Graph surrogate as supplement to sampling**

Grid is intrinsically a graph: admittance matrix (weighted Laplacian), fault propagates as message-passing, siting is a graph question.

**Graph neural networks** carry inductive bias for power systems:
- Trained on existing pilot runs: (graph, operating point, parameters, fault) → trajectory map
- Emulate CCT or stability label **cheaply**
- Full runs confirm results and cover surrogate's uncertain regions
- Active learning chooses next runs

**Graph surrogate is supplemental to sampling-based UQ**, not a replacement:
- Sampling defines what estimate means (aleatoric vs epistemic, design space coverage)
- Surrogate lowers per-sample cost and enables adaptive refinement

---

## Slide 13: GridKit End-User Findings

**Findings from using GridKit in anger:**

| Issue | Impact |
|---|---|
| **REPCA cannot do frequency droop** | No primary frequency response from IBR; ROCOF studies incomplete |
| **BusFault closing at event-on causes IDA fail** | At alphaP 0.8, 97 runs fail initialization (should pass at 0.2). Mechanism: pre-fault state as IDA initial guess. At constant-power load, Newton method conditioning degrades. |
| **LoadZIP has no current limit** | Constant-power load cannot self-limit; physical gap for constant-current stand-in. Unknown if this drives the IDA failure. |
| **Fault at constant-power load always fails** | 40 of 40 runs fail (every alphaP, including 0.2). Most interesting contingency for datacenter study is unreachable. |
| **Events are FAULT_ON/FAULT_OFF only** | No generator trip or load step; true power-imbalance event not modeled. |
| **Monitors default to empty** | Post-#577 cases carry `"mon": []` on every bus. Nothing PMU-observable recorded unless explicitly enabled. |

**Portability:** GridKit built on RHEL 8 runs unchanged on RHEL 9 (GLIBC ≤ 2.14, GLIBCXX ≤ 3.4.31 required; WECC240_validation passes in 0.64 s).

