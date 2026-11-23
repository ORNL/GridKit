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
| **ZIP load** | Load model splitting demand into constant-impedance (Z), constant-current (I), and constant-power (P) parts. Datacenter uses 80/20/0; wind stand-ins use 0/100/0 (pure constant-current, models controllable current source). |
| **No outcome** | A run with no CCT value: DCOPF dispatch, ACPF power flow, or dynamic simulation failed. Reported separately, never as stable or unstable. |
| **Placement** | Where the datacenter sits: none (baseline), or one of buses 6, 108, 120, 148. |

---

## Slide 3: The Illinois Synthetic Grid (ACTIVSg200)

**Grid Characteristics:**
- **200 buses, 246 branches, 49 generators** (43 thermal: 25 coal, 17 gas, 1 nuclear; 6 wind IBR)
- **40 synchronous machines** (34 thermal with dynamic models online, 6 wind units, 9 offline thermal units can be committed)
- **109 loads** across the network
- Hourly load and wind are synthetic one-year series (TAMU repository)

**Base case:** One-time build from TAMU files (topology, dispatch, dynamics); dispatch overwritten per condition; six wind units become negative load stand-ins (no inertia).

![ACTIVSg200: 200-bus Illinois synthetic grid topology showing branches (gray), loads (dark grey, sized by MW), IBR/wind generators (blue), and thermal generators (brown).](uq-subsections/figs/uq-pilot/ACTIVSg200_slides_topology_loads_ibr.png)

*Figure 3.1. ACTIVSg200 (200-bus Illinois synthetic grid): topology, loads, and generator mix.*

---

## Slide 4: Building One Operating Condition

**• Unit Commitment: Not Real UC**
- We do not solve a real UC problem. Results are **conditional on this sampling strategy**; applying real UC distributions is future work.
- We rank units by gen costs (nuclear always on); sample *n* at the quarter and three-quarter positions within the feasible band per hour.
- Feasible band per hour: lower bound = capacity to cover net load + 10% reserve; upper bound = committed minimums ≤ net load. At identical load & wind, inertia varies 1.56× across the band.

**• Dispatch: DCOPF → ACPF**
- We use PowerModels.jl for both DCOPF and ACPF solves. As OPF functionality becomes available in GridKit (under development), we will substitute it in. PowerModels.jl serves as a validation tool and enables exploration of different ACOPF relaxation formulations.
- DCOPF (lossless) sets economic dispatch active power; ACPF (with losses) solves voltages, reactive power, and enforces Q limits.

---

## Slide 5: Operating Condition Dynamics

**• Wind as Negative Load**
- IBR provide no synchronous inertia; model as controllable current source (WECC/Sandia standard).
- Hourly wind MW pinned in DCOPF; ACPF solves reactive output Q. Wind injection scales with |V|, correct for constant-current control.

**• The Datacenter**
- Added as load at one bus before dispatch; OPF sees and must supply it. ZIP model 80/20/0: 80% constant-power, 20% constant-current, 0% impedance.
- Constant-power load keeps drawing during faults, depressing voltage further → destabilizing. Placement: faultscreen-2 at bus 108 (100 MW); aleatoric-1/2 at buses 6, 108, 120, 148 (50–400 MW).

---

## Slide 6: Critical Clearing Time by Bisection

**• Stability Criterion: Centre of Inertia (COI) Pole Slip**
- Pole slip measured relative to inertia-weighted average rotor angle: $\delta_{COI}(t) = \frac{\sum_i H_i S_i \delta_i(t)}{\sum_i H_i S_i}$
- Run unstable when: $\max_{i,t} |\delta_i(t) - \delta_{COI}(t)| > \pi$ (any machine leaves collective inertial frame).

**• CCT Search: Bisection [0.05, 0.5] s**
- Fault applied at $t_{on} = 1.0$ s, cleared at $t_{off} = 1.0 + \text{duration}$. Tolerance 0.01 s, typically 8 simulations per run.
- Censored run: stable at 0.5 s cap → reported as CCT ≥ 0.5 s. Efficiency: record only 21 of 441 rotor angles → 13.8× storage reduction.

---

## Slide 7: Sampling Space and Three Campaigns

**• Campaign Design: Stratified by Net-Load Quintile**
- Three independent campaigns built identically, varying fault (1 vs 17), datacenter buses (1 vs 4), sizes (single vs range). Each hour commits at two positions: quarter and three-quarter of feasible band.
- Stratification: 8,760 hours split into 5 equal-count bins by net load; hours drawn randomly within quintiles. Each quintile represents 20% of year's hours.

**• Campaign Scale**

| Dimension | **aleatoric-1** | **aleatoric-2** | **faultscreen-2** |
|---|---|---|---|
| **Hours** | 20 (4 per quintile) | 60 (12 per quintile) | 50 (10 per quintile) |
| **Conditions** | 40 | 180 | 100 |
| **Faults** | 1 (bus 66) | 1 (bus 66) | 17 locations |
| **Total runs** | 840 | 3,780 | 3,400 |
| **Runs with CCT** | 439 (52%) | 2,155 (57%) | 3,170 (93%) |

---

## Slide 8: Reading the Curves — Kaplan–Meier and Censoring

**• Kaplan–Meier Estimator with Censoring**
- Handles censored runs (stable at 0.5 s cap, true CCT unknown): $\hat S(t)=\prod_{t_i\le t}\left(1-\frac{d_i}{n_i}\right)$, where $n_i$ = at-risk count (measured at or after $t_i$, plus censored).
- 95% CI from Greenwood variance: $P(\text{CCT}<t) \pm 1.96 \cdot \text{SE}$. Censored runs kept in denominator → valid inference.

**• Example: 200 MW at Bus 108, Fault at Bus 66 (37 Conditions)**
- Datacenter raises failure risk: baseline P(CCT < 0.25 s) = 0.49 → 0.73 with 200 MW. Median CCT drops from 0.27 s to 0.18 s.
- 95% Greenwood CI overlap at n=37; paired test at n=171 shows significant shift.

---

## Slide 9: Risk by Net Load

**• Risk Concentrates in High Net-Load Quintiles**
- aleatoric-1 data: all datacenter runs, 1 fault (bus 66). Risk sharp step between q1 (1,001 MW net load, P=0.16) and q2 (1,221 MW, P=0.73).
- At low net load, inertia high (many units at part load) → censored at 0.5 s. At high net load, fewer units carry more load → lower inertia, faster CCT.

**• Takeaway:** Risk is conditional on runs that solved; focus on high-net-load scenarios.

| Net-load quintile | Mean load (MW) | Runs solved | Censored | **P(CCT < 0.25 s)** |
|---|---|---|---|---|
| q0 | 809 | 75 of 160 | 83% | 0.07 |
| q1 | 1,001 | 82 of 160 | 72% | 0.16 |
| q2 | 1,221 | 82 of 160 | 0% | 0.73 |
| q3–q4 | 1,477–1,660 | 160 of 160 | 0% | 1.00 |

---

## Slide 10: Fault Screening — Which Faults Separate Operating Conditions

**• Fault Screening Results (faultscreen-2: 50 hours × 2 commitments = 100 conditions, 17 faults)**
- 4 separating faults (buses 66, 68, 94, 181): outcome changes with net load. 13 flat faults: censored everywhere, cannot separate conditions.
- Separating faults show rising curve in failure share vs quintile; flat faults show no discrimination.

**• Implications for 2000-Bus Scale**
- Contingency list grows ~115× (siting 11.5× × faults 10×). Screening cost grows 19–25× per run. Fault location explains >50% of variance.
- Need adaptive pre-screen to identify discriminating faults before full campaign; enables graph-based surrogate to reduce 2000-bus sampling cost.

---

## Slide 11: Scaling to Texas — ACTIVSg2000

**• Per-Run Scaling (10× buses, TAMU repository)**
- DCOPF, ACPF: 19.1×, 18.5× cost growth ($N^{1.27}$). Dynamic simulation: 25.3× growth ($N^{1.40}$). Storage: 9.6× for monitors (raw CSV, $N^{0.98}$).
- Brute force survives one order of magnitude of grid size, not two.

**• Projected aleatoric-1 on Texas (840 runs, ~7,560 simulations)**
- Serial: ~35 min (ACTIVSg200) → ~14.9 h (ACTIVSg2000). Parallel (16 workers, 1 node): ~2.2 min ideal → ~56 min ideal; actual (2.7× overhead): ~2.5 h.
- Monitor storage: 22 GB → ~220 GB. Solution: Screening + surrogate models + adaptive sampling for next phase.

---

## Slide 12: Next Steps — FY27 Tracks

**• Track 1: Spend Fewer, Better-Chosen Dynamic Runs**
- Screen contingencies cheaply (faultscreen-2 showed 13 of 17 censored everywhere). Reduce CCT search cost: start bisection from neighbouring dispatch (CCT correlated) or early termination detection (pole slip detectable before simulated window end).
- Stratified/importance sampling on top-ranked factors (fault location >50% variance, net load, placement) rather than uniform sampling.

**• Track 2: Graph Neural Network Surrogate (Supplemental to Sampling)**
- Trained on existing pilot runs: (graph, operating point, parameters, fault) → trajectory map. Emulates CCT/stability label cheaply.
- Full runs confirm results and cover surrogate's uncertain regions; active learning chooses next runs. Surrogate lowers per-sample cost and enables adaptive refinement while sampling defines design space coverage (aleatoric vs epistemic).

---

## Slide 13: GridKit End-User Findings

**• Critical Modeling Gaps**
- REPCA cannot do frequency droop → no primary frequency response from IBR. BusFault closing at event-on causes IDA fail (97 runs at alphaP 0.8; Newton method conditioning degrades with constant-power load).
- Fault at constant-power load always fails (40 of 40 runs); most interesting contingency for datacenter is unreachable. LoadZIP has no current limit → constant-power load cannot self-limit.

**• Operational Limitations**
- Events are FAULT_ON/FAULT_OFF only (no generator trip, load step, true power-imbalance event). Monitors default to empty; post-#577 cases carry `"mon": []` on every bus unless explicitly enabled.
- GridKit portability: built on RHEL 8, runs unchanged on RHEL 9 (GLIBC ≤ 2.14, GLIBCXX ≤ 3.4.31 required; WECC240_validation passes in 0.64 s).

---
