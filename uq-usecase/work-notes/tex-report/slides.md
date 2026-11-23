# Datacenter Load and Transient Security: Pilot Study Slides

---

## Highlight Slide: Overview

**Title:** Datacenter Load and Transient Security: An Aleatoric UQ Pilot on GridKit

**Question:**
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
- Placement location, load size, and operating condition interact; quantifying these dependencies requires sampling-based uncertainty analysis

**Study Approach:**
- Pilot on 200-bus Illinois synthetic grid with real hourly load/wind data
- Apply sampling-based aleatoric UQ to quantify stability margins
- Establish and validate sampling methodology for 2000-bus grids

**[Figure: Illustration of grid topology with IBR + datacenter siting locations]**

---

## Slide 2: Grid, Setup, and ZIP Load Model

**Title:** ACTIVSg200 and Datacenter Representation

**Grid Characteristics:**
- **200 buses, 246 branches, 49 generators** (43 thermal, 6 wind)
- **109 loads** distributed across the network
- One-year synthetic hourly load/wind series (joint distribution from TAMU)

**Datacenter Model — ZIP Load 80/20/0:**
- 80% constant-power (P): typical of IT load
- 20% constant-current (I): cooling motors
- 0% impedance (Z): negligible

**Placement Scenarios:**
- **Baseline:** No datacenter
- **Sited:** Four candidate buses (6, 108, 120, 148) × three sizes (25, 50, 100 MW)

**Power Flow Workflow:**
1. DCOPF (dispatch): Who generates? (lossless, active power only)
2. ACPF (power flow): What are voltages and reactive power?
3. Dynamics (phasor simulation): Does system survive the fault?

---

## Slide 3: Critical Clearing Time (CCT) Measurement

**Title:** How We Measure Stability

**CCT Definition:**
The longest fault duration (in seconds) the system can withstand before a generator loses synchronism with the grid.

**Bisection Search Method:**
- **Fault window:** [0.05 s, 0.5 s] (protection systems clear real faults well before 0.5 s)
- **Stopping criterion:** Tolerance = 0.01 s
- **Output:** Last stable duration probed (conservative lower bound)

**Stability Test:**
- Pole slip detection: generator rotor angle relative to center of inertia
- If pole slip occurs → generator loses synchronism → system fails

**Censoring:**
- If system stable at 0.5 s cap → **right-censored run** (true CCT ≥ 0.5 s)
- Censored runs are very stable; we know lower bound only
- Handled via Kaplan–Meier estimator

**Result:** One CCT value per operating condition + placement + fault combination

---

## Slide 4: Sampling Campaigns Design

**Title:** Two Independent Aleatoric Campaigns

**Campaign Structure:**

| Dimension | Aleatoric-1 | Aleatoric-2 | Faultscreen-2 |
|-----------|------------|------------|--------------|
| **Operating conditions** | 37 sampled | 171 sampled (independent) | Stratified by net load |
| **Fault location** | Fixed (bus 108) | Fixed (bus 108) | Varied (17 locations) |
| **DC size** | 0, 25, 50, 100 MW | 0, 25, 50, 100 MW | Fixed 100 MW |
| **Purpose** | Pilot proof-of-concept | Large-sample validation | Fault impact ranking |

**Key Design Choice:**
- Aleatoric-1 and aleatoric-2 are **disjoint samples** (no overlap in operating conditions)
- Both use same (placement, fault) to show sampling convergence
- Faultscreen-2 varies fault to identify dominant contingencies

**Outcomes:**
- **37 matched conditions (aleatoric-1):** P(CCT < 0.25s) baseline = 0.49, datacenter = 0.73
- **171 conditions (aleatoric-2):** P(CCT < 0.25s) baseline = 0.515, datacenter = 0.731
- **9 of 180 runs had no outcome** (dispatch/power flow/dynamics failure)

---

## Slide 5: Key Finding — CCT Shift with Datacenter

**Title:** Datacenter Load Increases Failure Risk

**Greenwood Confidence Intervals (95%):**

**Baseline (no datacenter):**
- P(CCT < 0.25 s) = 0.49 (CI: 0.325–0.648)
- Wider interval, higher uncertainty at smaller samples

**With 100 MW Datacenter at Bus 108:**
- P(CCT < 0.25 s) = 0.73 (CI: 0.587–0.873)
- Narrower interval, more stable

**Paired Difference (171 conditions):**
- Shift: **+0.22** (95% CI: 0.15–0.28)
- Consistent: Same magnitude shift across both campaigns

**Physical Mechanism:**
- Constant-power load increases fault current
- Faster voltage recovery post-fault is offset by larger transient swings
- Net effect: Reduced margin to synchronism loss

**Interpretation:**
Placing 100 MW constant-power load at bus 108 reduces the probability that the system survives a 0.25 s fault by roughly 24 percentage points.

---

## Slide 6: Fault Screening — Which Faults Matter?

**Title:** Faultscreen-2 Results: Contingency Ranking

**Main Finding:**
13 of 17 fault locations were **censored in every baseline run** (all CCT ≥ 0.5 s).

**Four Separating Faults:**
Only four faults reliably separated operating conditions (produced uncensored runs):
- Bus 148 (most dynamic impact)
- Bus 108 (siting location; natural stress point)
- Two additional buses identified by sensitivity

**Insight:**
- Grid topology strongly constrains which faults matter
- Pre-screening by graph centrality or admittance could cut simulation cost
- Most contingencies are non-binding for this system

**Cost Implication:**
On a 2000-bus grid, contingency explosion makes pre-screening essential. Graph neural networks (Track 2, FY27) aim to solve this.

---

## Slide 7: Risk by Net Load Quintile

**Title:** How Operating Condition Shapes Risk

**Net Load Quintile Analysis:**
Stratified comparison of baseline vs. datacenter across five net-load bins.

**Key Observation:**
- **Q1 (lowest net load):** Censoring highest (~83% baseline), CCT shift modest
- **Q5 (highest net load):** Censoring lowest (~0%), CCT shift largest
- **Trend:** Risk increase from datacenter is **load-dependent**

**Mechanism:**
- High net load → system already stressed
- Datacenter siting adds extra stress → larger CCT drop
- Low net load → system slack → datacenter impact damped by voltage support

**Planning Implication:**
Datacenter placement timing matters. Placing during high-load periods is higher risk.

---

## Slide 8: Censoring and Right-Skewed Stability

**Title:** Why Censored Runs Tell Us About System Safety

**What Censoring Means:**
- **37 of 180 runs** in aleatoric campaigns were censored (system survived 0.5 s fault)
- These runs are very stable — true CCT ≥ 0.5 s is unknown but large
- Kaplan–Meier estimator correctly accounts for this information

**Why It Matters:**
- Censoring is **not uniform:** Baseline censors more than datacenter
  - Baseline: 83% at Q1, 0% at Q5
  - Datacenter: Lower censoring across all quintiles
- Dropping censored runs would **inflate measured failure risk** (biased estimate)

**Statistical Approach:**
- Use Kaplan–Meier (KM) estimator [16] to handle censoring
- KM survival curve estimates P(CCT ≥ t) correctly
- Greenwood variance formula [19] gives 95% confidence intervals

**Result:**
Properly accounting for censoring, the datacenter shift is robust and reproducible across both campaigns.

---

## Slide 9: Challenges, Limitations, and GridKit Findings

**Title:** What We Learned About the Tool and Method

**Challenges Encountered:**

1. **IDA Initial Condition Solver Conditioning:**
   - Some fault scenarios cause initialization failure
   - Discriminator: fault-on event (at t=0) fails more than fault-off
   - Likely interacts with constant-power load representation

2. **Constant-Power Load Has No Protective Devices:**
   - LoadZIP model in GridKit has no current limit or undervoltage trip
   - Physically realistic (IT load cannot self-limit)
   - But creates modeling gap: unconstrained reactive demand

3. **Fault Types Limited:**
   - Only FAULT_ON / FAULT_OFF available
   - No generator trip or load-step contingency
   - Limits representation of true power imbalance events

4. **Contingency Cost Scales:**
   - At 200 buses: dispatch dominates cost
   - At 2000 buses: dynamic simulation dominates
   - Brute-force sampling becomes infeasible

**Key GridKit Validation:**
- Successfully built and ran on RHEL 9 (portability confirmed)
- Timing: 0.64 s per case on single node
- Ready for 2000-bus deployment

---

## Slide 10: Next Steps — FY27 Planning

**Title:** Scaling to 2000-Bus Grids and Beyond

**Track 1: Better Sampling (Deterministic Methods)**
1. **Fault screening:** Pre-screen 2000-bus contingencies via graph methods; focus on few high-impact faults
2. **CCT acceleration:** Exploit neighbor CCT correlation; early termination detection
3. **Variance stratification:** Decompose variance by fault location, net load, placement size; use importance sampling

**Track 2: Learning-Based Surrogates (Graph Neural Networks)**
- **Goal:** Build fast neural surrogates from existing campaign data
- **Why graphs:** Grid admittance matrix is a weighted Laplacian; faults propagate via message-passing
- **Benefit:** Identify high-variance parameter regions; automate importance-sampling guidance; reduce per-sample cost

**Task Portfolio (FY27):**
1. **Task 4.1:** Conclude 200-bus pilot with expanded sampling and fault pre-screening
2. **Task 4.2:** Scale epistemic + aleatoric campaigns to 2000-bus Texas grid
3. **Task 4.3:** Integrate epistemic-aleatoric UQ on 200-bus case
4. **Task 4.4:** ZIP load model parameter estimation via Bayesian neural networks
5. **Task 4.5:** Digital Twin prototype for sequential parameter refinement
6. **Task 4.6 (Stretch):** Graph NN surrogates for adaptive sampling guidance

**Deliverables:** Methodology framework, 2000-bus validation, digital twin prototype

---
