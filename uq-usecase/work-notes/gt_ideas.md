# Graph Transformer Ideas: ML on the Gen-Parameter UQ Data

Exploratory research-planning note. This is an **idea survey**, not a fixed experiment
protocol — it maps how graph transformers (GTs) could be trained/validated on the data
produced by the epistemic UQ pipeline (`uq_setup.ipynb`), surveys the relevant phasor-dynamics
ML literature, lays out four candidate problems, and sketches the additional data runs we would
need. Model choices, splits, and metrics are kept deliberately open at this stage.

Cross-references:
- Data description and pipeline: `notes_2.md` (uq_setup section), `notebooks/uq_setup.ipynb`
- Motivation for fault-response variability: `uq_plan.md` (Task 0b)
- Aleatoric operating-point generation: `uq_plan.md` (Task 0 / aleatoric track)

---

## 1. Motivation and scope

The epistemic UQ track already generates large ensembles of GridKit dynamic simulations with
sampled generator parameters (currently inertia constant H). Each run is, in effect, a
**(graph, parameters, fault) → trajectory** mapping evaluated by an expensive ODE/DAE solver.
That is exactly the kind of structured, graph-native, physics-labeled dataset that graph
transformers are built to exploit.

In parallel, the **aleatoric track** is generating a large operating-point dataset — on the
order of **8760 hours × multiple years** of load/dispatch/wind realizations, each solved to a
valid AC power flow (see `uq_plan.md`, Task 0). Together the two tracks give us variation along
*both* axes: **operating point** (aleatoric) and **dynamic parameters** (epistemic). The full
map is **(graph, operating point, parameters, fault) → trajectory**, which is a much richer
supervision signal for graph learning than either track alone.

The question this note explores: **what practically interesting, publishable ML problems can we
pose on this data — using what we already have plus a couple of modest new runs — and where do
graph transformers specifically add value over simpler models?**

Scope of this note:
- Broad survey of candidate problems (A–E), presented even-handedly, across **both** the
  epistemic (parameter) and aleatoric (operating-point) tracks.
- Idea-level depth: framings, literature grounding, pros/cons, data needs. No committed model
  architecture, train/val split, or metric suite yet.
- Includes rough specs for the new data runs that some of these problems require.
- Treats the aleatoric operating-point dataset as a first-class data source, not just a
  confounder — it both broadens the training distribution and enables its own ML problems.

---

## 2. Data we already have

From the `uq_setup` pipeline, per sample:

**Epistemic track (dynamic parameters)** — from the `uq_setup` pipeline, per sample:

| Item | Detail |
|------|--------|
| Graph | Grid topology: buses = nodes, branches = edges (R, X, B, rating) |
| Node time series | Bus `Vr`/`Vi` (→ `Vm`/`Va`); generator `omega` (speed deviation), `delta` (rotor angle) |
| Sampled parameters | Genrou `H` for a fixed subset of generators (LHS; Gaussian σ=12% or uniform ±%) |
| Fault | A single fixed bus per case (currently one fault location per grid) |
| Grids | Hawaii40 (37 buses), Illinois / ACTIVSg200 (200 buses) |
| Scale to date | ~22,000 dynamic simulations total (epistemic runs only) |
| Serialization | Parquet, stacked or per-run (see `notes_2.md`) |

**Aleatoric track (operating points)** — currently being generated in parallel:

| Item | Detail |
|------|--------|
| Variation | Hourly load / dispatch / wind, ~8760 h × multiple years |
| Per point | A valid AC power-flow solution → `case.json` initial condition (see `uq_plan.md`, Task 0) |
| Fleet | Selected spinning gens thinned (lower inertia → more interesting fault response) |
| Label | Operating-point features (per-bus P/Q injections, dispatch, wind level, hour/season) |
| Role for ML | Steady-state graph state that *precedes* each dynamic run; conditions the transient |

Key properties that make this ML-ready:
- Every sample shares a **fixed topology** per grid (graph structure constant; node/edge
  features and labels vary). This is the easy regime for graph learning.
- Labels are **exact** (sampled H and known operating-point inputs), so supervised inverse
  problems have ground truth for free.
- Trajectories are physically consistent (same solver), so a surrogate has a well-defined target.
- The two tracks are **composable**: an aleatoric operating point can be paired with an
  epistemic parameter sample to form a joint (scenario × parameter) design, exercising the model
  over realistic operating conditions rather than a single static dispatch.

---

## 3. Literature landscape: GT/GNN in power-system dynamics

Established application buckets (roughly by maturity):

1. **Transient stability assessment (TSA).** Spatio-temporal GCN / graph-attention models
   classify stable vs. unstable post-fault from PMU-like signals. The most mature area; lots of
   prior art, so novelty bar is high.
2. **Dynamic parameter / state estimation.** GNNs to estimate inertia H and damping D from
   measured response. Very topical as system inertia falls with rising inverter-based resources;
   directly aligned with what we sample.
3. **PMU / sensor placement.** Optimal placement for observability and event localization;
   attention- and RL-based approaches emerging. Naturally couples to estimation quality.
4. **Fault / event localization.** Node-level "which bus faulted / where did the event
   originate" from measured transients. Clean supervised node-classification framing.
5. **Fast surrogates / emulation.** GNN and Graph-Neural-ODE trajectory predictors that replace
   or accelerate the dynamic solver, enabling cheap forward UQ and screening.

**Why a graph *transformer* rather than a plain GNN?**
- **Global attention** captures long-range electrical coupling (a fault propagates well beyond
  1–2 hops), which message-passing GNNs approximate only through many layers.
- **Positional/structural encodings** (Laplacian eigenvectors, random-walk PE) inject graph
  geometry, letting attention respect electrical distance.
- **Avoids over-smoothing** that limits deep MPNNs on larger grids like ACTIVSg200.
- Mature backbones exist: Graphormer, GraphGPS (hybrid MPNN + global attention), Exphormer,
  TokenGT.

---

## 4. Candidate problems (A–D)

| # | Problem | Input → Label | Needs new data? | Comment |
|---|---------|---------------|-----------------|---------|
| A | Inverse parameter inference (UQ of gen params) | trajectories → per-gen H (+ uncertainty) | No | Natural inverse of what we generate |
| B | Fast forward surrogate | (topology + operating point + H + fault loc) → trajectory / QoIs | No | Replaces the solver; enables cheap forward UQ |
| C | Fault localization | trajectories → faulted bus | Yes (fault-location sweep) | Clean node-classification task |
| D | Identifiability-driven sensor placement | inference-model attention → minimal sensor set to recover H within tolerance | No (builds on A) | Ties UQ and placement together; less crowded |
| E | Operating-point-conditioned dynamics | (topology + operating point + H + fault) → transient risk / QoIs | Uses aleatoric track | Screens which of 8760×yr operating points are dynamically fragile |

### A. Inverse parameter inference (UQ of generator parameters)
Learn the mapping from observed transients back to the generator parameters that produced them
(regress H per generator, ideally with calibrated uncertainty). This is the exact inverse of the
forward map we already sample, so ground truth is free.
- **Pros:** Directly uses existing data; ground-truth labels; strong story on falling-inertia
  relevance; identifiability analysis (which H are recoverable, and from where) is scientifically
  interesting.
- **Cons:** H may be weakly identifiable for generators electrically far from the fault; needs an
  uncertainty-aware head to be honest about that.

### B. Fast forward surrogate
Learn `(graph, operating point, H vector, fault location) → trajectory` or `→ QoIs` (frequency
nadir, ROCOF, max rotor-angle swing, voltage-recovery time). Replaces GridKit for screening and
forward UQ. Conditioning on the aleatoric operating point lets one surrogate cover the whole
8760×year envelope, not a single dispatch.
- **Pros:** Very practical (orders-of-magnitude speedup enables large forward UQ / optimization);
  QoI regression is a clean, low-risk first target; pairs naturally with A (forward vs. inverse)
  and E (operating-point conditioning).
- **Cons:** Full trajectory prediction is harder than QoI regression; must demonstrate speedup
  *and* accuracy to be compelling.

### C. Fault localization
Learn `trajectories → which bus faulted`. A canonical, well-scoped supervised node task.
- **Pros:** Clean labels; standard evaluation (accuracy / top-k / distance-to-true-bus); good for
  demonstrating GT's long-range attention.
- **Cons:** Requires a **fault-location sweep** run (current data has one fault bus per grid), so
  it is gated on new data generation.

### D. Identifiability-driven sensor placement
Use the inference model from A (e.g., its attention / sensitivity to node inputs) to find the
minimal set of measured buses needed to recover H within a tolerance.
- **Pros:** Novel angle that unifies UQ and sensor placement; high practical value (PMUs are
  expensive); builds directly on A with no new simulation data.
- **Cons:** Methodologically the most exploratory; needs a defensible link between model
  attention/sensitivity and true observability.

### E. Operating-point-conditioned dynamics (aleatoric screening)
Use the aleatoric operating-point dataset (8760 h × years of solved power flows) as graph state
feeding a dynamic prediction, and learn which operating points are dynamically fragile — i.e.
`(graph, operating point, H, fault) → transient risk / QoIs`, then rank/screen the year(s) of
operating points by post-fault severity.
- **Pros:** Directly exploits the large aleatoric dataset being generated now; strong operational
  story (screen a whole year of operating conditions cheaply for stability risk); combines
  aleatoric + epistemic uncertainty in one model; heavy use of the (scenario × parameter)
  composability noted in §2.
- **Cons:** Depends on the aleatoric track maturing and on pairing operating points with dynamic
  runs; label definition ("fragile") needs care; largest data-management effort of the five.

**Framing note.** A + B form the most natural forward/inverse pair on the *existing* epistemic
data and are the lowest-risk starting point; E extends B/A across the aleatoric operating-point
envelope as it comes online; D is the most novel stretch and reuses A's model; C is a clean
second supervised task that unlocks once the fault-location sweep exists. This note treats all
five as live options rather than committing to one.

---

## 5. Planned new data runs (rough specs)

Three extensions to the current runs would substantially widen the ML problem space. R1/R2 stay
within the epistemic (parameter) track; R3 draws on the aleatoric operating-point track being
generated in parallel.

### Run R1 — Fault-location sweep
Vary the **fault bus** across the grid, holding the H-sampling scheme fixed per fault.
- **Purpose:** Unlocks problem C; gives A/B/D generalization across fault locations; lets a
  surrogate condition on the fault node.
- **Rough spec:**
  - Hawaii40: sweep all 37 buses (or all valid fault buses).
  - Illinois / ACTIVSg200: a representative subset, e.g. ~40–60 buses spanning electrical
    distances / voltage levels (full 200-bus sweep later if warranted).
  - For each fault bus, a moderate H-sample count (e.g. a few hundred) — total sized to fit the
    same SLURM fan-out already used.
- **Label added:** fault-bus id per run (for C), fault node as a conditioning input (for B).

### Run R2 — Expanded perturbed-generator / parameter set
Sample parameters for **more generators** and optionally **more parameter types**.
- **Purpose:** Richer inference targets for A and D; probes identifiability as the unknown-vector
  grows.
- **Rough spec:**
  - Extend H sampling from the current fixed subset toward all (or most) generators.
  - Optionally add damping D and/or key machine reactances (e.g. `Xd'`) as additional sampled
    dimensions — start with H-only expanded, add others as a second phase.
  - Keep sampling distribution consistent with current practice (Gaussian σ≈12%) for comparability.
- **Label added:** larger H vector (and optionally D / reactances) per run.

### Run R3 — Aleatoric operating-point coupling
Pair the aleatoric operating points (8760 h × years, each a solved AC power flow) with dynamic
runs, so each dynamic simulation starts from a realistic hourly operating condition.
- **Purpose:** Unlocks problem E; broadens the training distribution for B (surrogate valid
  across the operating envelope, not one dispatch); lets A/D be tested for robustness to
  operating-point drift.
- **Rough spec:**
  - Sample a representative subset of operating points (e.g. stratified by season / load level /
    wind level) rather than all 8760×years dynamically at first.
  - For each selected operating point, run a fault (fixed or from R1's sweep) with a modest
    H-sample count; scale up as the aleatoric dataset and SLURM budget allow.
  - Reuse the Task-0 `.m → case.json` pipeline (offline gens absent from JSON) for initial
    conditions.
- **Label added:** operating-point features (per-bus P/Q, dispatch, wind level, hour/season) as
  conditioning inputs; post-fault QoIs / fragility label for E.

**Note on composition:** R1 (fault location), R2 (parameter set), and R3 (operating point) are
orthogonal axes. Full-factorial coverage is infeasible, so realistic designs will sample the
joint space (e.g. LHS over parameters × stratified operating points × sampled fault buses).

---

## 6. Dataset design (survey level)

Not committing to specifics yet, but the shape is clear:

- **Graph:** buses = nodes; branches = edges with features (R, X, B, rating). Generators either
  attached as node features on their bus, or represented as separate generator-nodes linked to
  their bus (to be decided).
- **Node features:**
  - Static: base voltage, bus type (PV/PQ/slack), generator-present flag, nominal H / S rating.
  - Operating point (aleatoric): per-bus P/Q injection, generator dispatch, wind level, and
    time context (hour / season) for the pre-fault steady state.
  - Dynamic: windows of the measured trajectory (`Vm`, `Va`, and for generator buses `omega`,
    `delta`).
- **Positional / structural encodings:** Laplacian eigenvectors or random-walk PE to give the
  transformer graph geometry.
- **Labels (per problem):** H vector (A); trajectory or QoIs (B); fault-bus id (C); minimal
  sensor set / observability target (D); post-fault fragility / QoIs conditioned on operating
  point (E).
- **Splits (ideas):** by sample index within a grid; a cross-grid generalization test (train
  Hawaii → test Illinois) for topology transfer; and, for the aleatoric axis, a temporal split
  (train on some years/seasons, test on held-out ones) to check operating-point generalization.

---

## 7. Models and baselines (survey level)

- **Graph transformer backbones:** Graphormer, GraphGPS (MPNN + global attention hybrid),
  Exphormer, TokenGT.
- **Temporal handling:** trajectories via time-as-tokens, a temporal transformer head, or
  windowed features — to be chosen per problem.
- **Baselines (to show GT actually helps):** GCN / GAT / MPNN (graph but no global attention);
  and non-graph MLP / CNN on stacked per-node signals (no graph structure). The GT should beat
  both to justify the complexity.

Evaluation *directions* (not fixed metrics yet):
- A: per-generator H error + calibration; identifiability vs. graph distance from fault.
- B: trajectory RMSE / QoI error; wall-clock speedup vs. GridKit.
- C: localization accuracy / top-k / distance-to-true-bus.
- D: sensors-vs-accuracy trade-off curve; recovered-H error at fixed sensor budget.
- E: fragility ranking quality (e.g. how well the top-k riskiest operating points are recovered);
  QoI error across the operating envelope.

---

## 8. Phased roadmap

1. **Prototype (small):** pick one problem (A or B is lowest-risk) on Hawaii40 with existing
   data; establish the graph dataset builder and a GT-vs-baseline harness.
2. **Scale:** move to Illinois / ACTIVSg200; confirm the approach holds at 200 buses.
3. **New data:** generate R1 (fault sweep) — enables C and improves A/B/D generalization; then R2
   (expanded gens/params) for richer inference; and R3 (aleatoric operating-point coupling) as
   the parallel operating-point dataset matures.
4. **Operating-point conditioning:** fold the aleatoric operating points into B and take on E
   (screen operating conditions by dynamic fragility).
5. **Cross-grid:** train on one grid, test on another, to probe topology generalization.
6. **Stretch:** D (identifiability-driven sensor placement) on top of the matured A model.

---

## 9. Open questions and risks

- **Identifiability.** How recoverable is H for generators far from the fault, and from limited
  sensors? This is the crux for A and D — and scientifically the most interesting part.
- **Simulation cost.** R1/R2/R3 add solver runs; the joint (parameter × operating point × fault)
  space is huge, so sweeps must be sampled and sized to the existing SLURM budget.
- **Aleatoric–epistemic coupling.** Pairing operating points with parameter samples and faults
  raises data-management and bookkeeping complexity (which `case.json` × which H sample × which
  fault); the composition needs a clean indexing scheme.
- **Fixed vs. variable topology.** Everything here assumes fixed topology per grid; cross-grid
  transfer and (later) topology perturbations are harder and open. Note the aleatoric track
  thins the generator fleet (offline gens absent from JSON), so effective topology can shift
  between operating points — the model must handle generator on/off status.
- **Point vs. probabilistic inference.** A and D are most valuable with calibrated uncertainty,
  which raises the modeling bar versus plain regression.
- **Novelty positioning.** TSA (bucket 1) and basic estimation (bucket 2) are crowded; the more
  defensible contributions are identifiability-aware inference (A), the surrogate's forward-UQ
  utility (B), operating-point fragility screening at scale (E), and the placement angle (D).
