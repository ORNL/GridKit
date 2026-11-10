# IDA initial-condition failures under a constant-power datacenter

**Status:** hunting wrapped up. Main hypothesis supported (§11): the IC failures match
the fold sets exactly on all three dispatches tested (§10.8). Not filed as a GridKit issue.
**Owner:** isatkaus.

This file holds every finding on IDA initial-condition (IC) failures in the
datacenter pilot. The plan doc and the pilot report link here and carry only
summaries.

## 0. Read first

**Where the hunt stands.** The main IDA initial-condition (IC) failures in the 0.8
(and 1.0) campaigns are explained by working hypothesis H1 (§6, §11): a fault-on
**fold** (saddle-node, "nose") of the network algebra. There is no fault-on equilibrium
on the branch reached from the pre-fault state, so IDA has no consistent initial
condition to find. The hunt results are in §10, the hypothesis and its tests in §11, and
every script and output path in §12. Terms are defined in §9.

**Verified (tool output, §10):**
- At 200 MW the folded buses are **exactly** the IDACalcIC failures at alphaP 0.2 (3
  buses), 0.8 (18) and 1.0 (21). No exceptions in either direction.
- The IDASolve failures (2, 3 and 2 buses) are not folds at 200 MW. Three of them
  (buses 2, 134, 142 at 0.8) fold at 250 MW: integration failures near a fold.
- **Cross-dispatch (§10.8):** the same exact match holds on h2151 (15 folds = 15 IDACalcIC)
  and h4770 (16 = 16). Only the IDA IC failures are tested this way.
- Folds grow with alphaP and with MW, and shrink as the datacenter moves towards
  constant impedance (alphaZ). Fold sets are nested in all three cases.
- Bus 148 folds at 0.8 and 200 MW but not with constant impedance.

**Working hypothesis H1 (supported, not proven).** A constant-power share moves the
fault-on nose so that a bolted fault (X = 0.01 pu) at more buses lies beyond it. At
those buses no post-fault steady state exists on the reachable branch, so IDA's IC
solve fails. This is a physical limit of the model at that fault, not an IDA defect.

**What it means.** (a) At those buses the post-fault system has no steady state the
solver can reach. That is a censored outcome, not solver noise, and should be modelled
as such (see report §7.5). (b) The X = 0.01 bolted-fault convention places many buses
close to the nose; a finite-impedance fault would change the set. (c) The datacenter's
constant-power share is what moves bus 148 across the nose, which is the reason for
high-P modelling in the first place.

**Not resolved.** (i) Bus 146 has a converged PF fault-on solution (§7 item 7) but
folds on the continuation branch. The PF root may sit on another branch; not checked.
(ii) The IDASolve failures are not explained by H1. The three at 0.8 on h1879 are
near-folds at 200 MW (folds at 250 MW). Bus 117 fails IDASolve on h2151 and h4770 where
it is not a fold at all (§10.8). The user has asked not to chase these. (iii) The X = 0.01 convention and the lack of an exciter limiter in the modelremain caveats (§10.9).

## 1. The problem

The datacenter is modelled as a `LoadZIP` with `alphaP` (constant power) and
`alphaI` (constant current), with `alphaZ = 1 − alphaI − alphaP` derived by GridKit
(LoadZIP README). With a bolted-like fault shunt on a bus, the IDA IC solve fails
(`IDACalcIC`, or `IDASolve` later in the run) more often as `alphaP` rises:

| campaign (sited rows) | alphaP | ok | sim_failed | of which IDACalcIC / IDASolve |
|---|---|---|---|---|
| `faultscreen1` | 0.2 | 552 | 77 | 62 / 15 |
| `faultscreen1_ap08` | 0.8 | 432 | 197 | 176 / 21 |
| `faultscreen1_ap1` | 1.0 | 429 | 200 | 191 / 9 |
| `aleatoric1` | 0.2 | 478 | 0 | — |
| `aleatoric1_ap08` | 0.8 | 417 | 61 | 57 / 4 |
| `aleatoric1_ap1` | 1.0 | 410 | 68 | 63 / 5 |

Classes come from each condition's `sim.log` in scratch. Every sim failure in these
six campaigns is an IDA error. `aleatoric1` has no failures at 0.2.

Records do not carry the IDA class; the class was read from `sim.log` files.

## 2. Reproduction tools

- `scripts/ida_zoom.py`: IDA only. Stages `bus` (BusFault at every bus, controls with
  no fault), `solver` (integrator and IC options), `mw` (datacenter MW sweep). Takes
  `--src` (a failing dispatch case dir), `--alphas` (default `0.2,1.0`), `--fault-bus`,
  `--buses`. Writes `/scratch/isatkaus/gridkit-campaigns/ida_zoom/<name>/results.csv`.
- `scripts/pm_fault_check.py`: PowerModels.jl AC power flow on the same faulted case,
  datacenter MW sweep, warm and flat starts. Outputs under `ida_zoom/pf_mw/`.
  Julia via `pm_utils.resolve_julia()`.
- Julia: `pm_utils.resolve_julia()` uses `/nopt/nrel/apps/cpu_stack/software/julia/julia-1.12.1/usr/bin/julia`.
  `module avail julia` lists only 1.12.6; PATH may pick `~/bin/julia` (1.12.5).

## 3. Verified results (IDA, one run per bus)

Source: `ida_zoom/ic181`, `ic181_solver`, `ic181_mw`. Fault at bus 181 unless
stated.

- **Bus sweep (`ic181`, 200 buses).** Fault applied at t = 1.0 s (`event`) and at
  t = 0 (`closed`) give identical failure sets. `alphaP` 0.2: 49 failing buses
  (IDACalcIC 40, IDASolve 9). `alphaP` 1.0: 153 failing buses (IDACalcIC 125,
  IDASolve 28). The 0.2 set is a subset of the 1.0 set. With an **open** fault device
  there are no failures at either `alphaP`. That control is not clean: `event` mode
  differs from it in the event schedule as well as the fault state. See §3a.
- **Solver controls (`ic181_solver`, buses 1, 66, 94, 108, 148, 181).** Changing the
  consistent-IC type, tolerances or maximum order turns some IDACalcIC failures into
  IDASolve failures but rescues no `alphaP` 1.0 failure. Zeroing the datacenter removes
  every failure. At `alphaP` 0 (constant current) buses 1, 94, 108, 148 and 181 pass.
  Bus 66 at `alphaP` 0 gave `other`, which did not reproduce on rerun. Treated as
  transient, not explained.
- **Fault impedance.** At bus 181 and `alphaP` 1.0, fault X = 0.1 pu passes. Buses 1
  and 148 still fail at X = 0.1.

### 3a. Event-schedule discriminator (`discrim181`, `ida_discriminate.py`)

Fault at bus 181 case (dispatch `faultscreen1_ap1/h6483_n21_dc148-200/b181_d0.1`),
`alphaP` 1.0, 200 MW datacenter at bus 148. Result per bus (status from `sim.log`):

| mode | fault state at t = 0 | event schedule | tmax | bus 1 | bus 148 | bus 181 | bus 94 |
|---|---|---|---|---|---|---|---|
| `none` (no BusFault) | — | — | 0.2 | ok | ok | ok | ok |
| `open` | open | none | 0.2 | ok | ok | ok | ok |
| `open_noev_t12` | open | none | 1.2 | ok | ok | ok | ok |
| `open_ev_t02` | open | `fault_on` at 1.0 | 0.2 | IDACalcIC | IDACalcIC | IDACalcIC | ok |
| `event` | open | `fault_on` at 1.0 | 1.2 | IDACalcIC | IDACalcIC | IDACalcIC | ok |
| `closed` | closed | none | 0.2 | IDACalcIC | IDACalcIC | IDACalcIC | ok |

Verified from code (`application/PhasorDynamics/DynamicSimulation.cpp`, lines 51–69):
the run calls `initializeSimulation(0.0)`, then `runSimulation(event.time)` (run to
t = 1.0 whatever `tmax` says), sets the fault on, then calls
`initializeSimulation(event.time)`. So `open_ev_t02` fails at the **t = 1.0
re-initialization with the fault closed**. Its t = 0 solve is the same as in
`open_noev_t12`, which passes. The `IDACalcIC` log line does not print the time.

Consequence: the failing solve is a fault-on algebraic solve in both `closed` (at
t = 0) and `event` (at t = 1.0). The "inactive device fails at t = 0" framing in the
report's §8 draft is wrong: the failing solve is the post-fault one.

Also from code (`GridKit/Solver/Dynamic/Ida.cpp`, lines 231–246): `initializeSimulation`
calls `IDAReInit` with the current state `yy_`, `yp_` before `IDACalcIC`. At the event
the current state is the pre-fault state, so the IC guess is the pre-fault state.
- **Datacenter MW (`ic181_mw`, buses 1, 94, 148, 181; 25 to 600 MW).** Largest MW with
  an OK initial condition, `alphaP` 0.2 / 1.0: bus 1 100 / 25; bus 148 100 / 25;
  bus 181 250 / 150; bus 94 400 / 250.

## 4. PowerModels.jl cross-check (`ida_zoom/pf_mw/`)

Same fault shunt (G and B from the BusFault stamp, on the base `ac.m`). Datacenter P
and Q set to the MW level at bus 148, constant power (equivalent to `alphaP` 1.0).
PM.jl does **not** enforce generator Q limits here: the bus-181, 150 MW point has 26
Q-limit violations and `vmin` 0.22 pu.

| fault bus | PF converges (warm) | IDA IC OK, `alphaP` 1.0 | comparison |
|---|---|---|---|
| 94 | 25 to 250 MW (flat: 25 to 200) | 25 to 250 MW; IC at 300 | matches |
| 181 | 25 to 150 MW (flat: 25 to 125) | 25 to 150 MW; IC at 200 and above (175 not run) | matches to grid resolution |
| 1 | none (infeasible at 25 MW) | 25 MW OK, IC at 50 and above | **mismatch** |
| 148 | none (infeasible at 25 MW) | 25 MW OK, IC at 50 and above | **mismatch** |

The no-fault control converges to 250 MW (warm) and fails at 300 MW.

Note: PF does not depend on `alphaP`. It tests only the constant-power limit.

## 5. Campaign evidence

- **`cct181` (fault 181, 14 datacenter buses).** At `alphaP` 0.2, 9 of 224 dynamic
  runs fail: bus 42 at 200 MW (8 runs) and bus 39 at 200 MW (1 run). At `alphaP` 1.0,
  44 of 224 fail: buses 39 and 42, all 8 runs each at 100 MW and at 200 MW. So failures
  are not purely `alphaP`-driven. Per-record IDA class is not stored (scratch batch
  logs: 57 IDACalcIC and 3 IDASolve messages).
- **`faultscreen1_ap1`, `aleatoric1_ap1`.** The 51 and 322 dispatch-infeasible
  conditions match the 0.2 baselines. `alphaP` does not change which conditions are
  dispatched.
- **Same-run CCT.** At faults 94 and 68, 16 `alphaP` pairs gave identical CCT to four
  decimals. Pole slip occurs within the first second.
- **80/20/0 re-runs (`faultscreen1_ap08`, `aleatoric1_ap08`).** Failures at 0.8 are
  concentrated at:
  - `faultscreen1`: fault buses 1, 8, 141 and 148 (37 each at 0.8, 3 each at 0.2 for
    buses 1 and 8); datacenter buses 108 (112 failures) and 148 (85).
  - `aleatoric1`: fault 66 only, datacenter buses 108 and 120, at 200 and 300 MW.
  No condition that fails at 0.2 succeeds at 0.8. The CCT shift on conditions that
  succeed at both is at most 0.007 s (faultscreen) and 0.042 s (aleatoric).

## 6. What is established and what is not

**Verified:**
- Every sim failure in the sited campaigns is an IDA error (IDACalcIC or IDASolve).
- Failures need a fault-on algebraic solve: the fault is closed at t = 0, or it closes
  at an event (§3a). An open fault with no event gives none.
- Failures disappear when the datacenter is zeroed, and at `alphaP` 0 for the buses
  tested.
- IDA IC type, tolerances and max order do not rescue `alphaP` 1.0 failures.
- At buses 94 and 181, the PF limit matches the IDA IC limit (grid resolution).
- At buses 1 and 148, PF is infeasible at 25 MW where IDA IC succeeds.
- The PF and IDA checks use the same dispatch (`ac.m`, `case.json` built from it;
  base load 27.6 MW at bus 148 plus the sweep MW).

**Working hypothesis 1 (not confirmed):** at buses 94 and 181 the IDA IC failure
coincides with the loss of a constant-power AC solution on the faulted network. That
is a loadability limit, not an IDA defect. Falsified by a converged PF at an
IDA-failing point with Q limits respected, or a PF failure at an IDA-OK point on the
same network.

**Working hypothesis 2 (not confirmed):** the failing solves are event-time
re-initializations, and the IC guess is the pre-fault state (§3a, `Ida.cpp`). The
pre-fault state can be far from a post-fault solution, so Newton may fail from it even
when a fault-on solution exists. That would be an IDA initial-guess issue, not a
missing equilibrium. Hypotheses 1 and 2 predict different things: 1 predicts failures
at fault-on points that PF also rejects; 2 predicts failures at points PF accepts.

**Not established:**
- Whether the failures at buses 1 and 148 are physical or a PF setup issue. PM.jl does
  not enforce Q limits, and those buses are infeasible in PF even at 25 MW.
- Whether the 0.8 failures are IDA defects or physical loadability limits. Now
  supported as physical (H1, §10–§11), but not proven: the fold scan and the IDA sweep
  agree at 200 MW on h1879, and on h2151 and h4770 as well (§10.8).
- Whether the 0.8 failure set is nested in the 1.0 set. It is, on h1879 at 200 MW
  (§7 item 5).
- The cause of the bus 66 `other` error at `alphaP` 0.

**Conclusion, as far as the data go (superseded by §0 and §11):** the 0.8 failures are
IDA-reported non-convergence that coincides with a fold of the fault-on algebra at
every bus tested on h1879. The data do not show an IDA defect.

## 7. Proposed next tests

1. ~~**PF on the exact dispatch.**~~ **Done, operating point ruled out.** The PF check
   uses the dispatch directory's `ac.m`. The IDA `case.json` is built from that same
`ac.m` (`campaign.py` line 956). Bus 148 has PD 227.6 MW in the `h6483` dispatch
`ac.m` (used for the §4 and §7 item 2 checks); the 200 MW datacenter is split out, leaving a
27.6 MW base load. The h1879 reference case used for the hunt (§10) differs: bus 148
PD is 238.296 MW there (38.3 MW host plus 200 MW), verified in `ac.m`. The PF and IDA start from the same dispatch, and the bus 1 and 148 PF mismatch
is not an operating-point mismatch. It is still unexplained.
2. **PF with Q limits enforced (PV to PQ switching). Done, inconclusive.**
   `scripts/ida_pvpq.jl` (PowerModels.jl, separate from `pm_solve.jl`) with driver
   `scripts/ida_pvpq.py`. Output `ida_zoom/pvpq/results.csv`. Warm and flat starts give
   the same result. The fault-181 dispatch, datacenter at bus 148, MW grid 25 to 600.
   - **Buses 1 and 148:** plain PF (no Q limits) is `LOCALLY_INFEASIBLE` at every MW from
     25 to 600, including 25 MW where IDA IC is `ok`. Switching never runs there, so this
     test cannot say whether Q limits remove the infeasibility. Neither bus has a
     generator. Bus 1 is a load-only PQ bus, connected by a line to bus 119 in `ac.m`.
   - **Buses 94 and 181 and no fault:** the plain PF boundary matches IDA to grid
     resolution. Bus 94: PF ok through 250 MW, fails from 300; IDA ok through 250, IC from
     300. Bus 181: PF ok through 150 MW, fails from 175; IDA ok through 150, IC from 200.
   - **Switching:** with no fault, 25 to 100 MW, plain PF converges with 6 Q-limit
     violations. Switching converges in 2 iterations with 7 generators switched and 0
     residual violations. At bus 94, 25 to 100 MW, the same switching makes PF
     `LOCALLY_INFEASIBLE` at iteration 1, while IDA IC is `ok`.
   - **Reading:** the switch is one-way (no PQ to PV reversal), so bus 94 infeasibility
     after switching is an approximation artifact and cannot be trusted. The test does
     not discriminate H1 from H2. Plain PF also disagrees with IDA at buses 1 and 148 in
     the direction H2 does not predict (PF fails where IDA IC succeeds). That mismatch is
     still unexplained. Fault-shunt conversion to PM is one candidate, not yet checked;
     the test 3 result below locates the failure threshold in fault strength; the PF vs
     IDA mismatch at buses 1 and 148 is still open.
3. **Fault-strength sweep, closed mode. Done (job 19794871).** `scripts/ida_xsweep.py`,
   buses 1, 148 and 94 at `alphaP` 1.0 and 200 MW, fault X from 1.0 down to 0.01 pu.
   Output `ida_zoom/xsweep181/results.csv`. The reference case's fault is X = 0.01 pu,
   R = 0 (`case.json`), so the sweep includes the reference case.
   - Buses 1 and 148: ok at X = 1.0 and 0.5. IDASolve at X = 0.2 (IC found, integration
     fails). IDACalcIC at X = 0.1, 0.05, 0.02, 0.01.
   - Bus 94: ok at every X, down to 0.01.
   - Reading: the failure set is a monotone threshold in fault strength at buses 1 and
     148. IDASolve at X = 0.2 means an IC was found there, so the first failure is in the
     integration, not the IC. A monotone threshold is what hypothesis 1 predicts, but it
     also fits hypothesis 2 (a larger jump from the pre-fault guess for a stronger fault).
     Not decisive. A finer sweep between 0.5 and 0.1 would locate the threshold.
4. **Event-time reproducer. Done (job 19794872).** `scripts/ida_event_repro.py`:
   self-contained case dirs, event time t in {0.5, 1.0, 2.0} at buses 1 and 94. Output
   `ida_zoom/event_repro/results.csv`.
   - Bus 1: IDACalcIC at every event time (0.5, 1.0, 2.0). Bus 94: ok at every event time.
   - No-fault controls (`noev`) at buses 1 and 94: ok. Note this control removes the
     fault entirely; the fault-present-from-t=0 control is the `closed` mode (IC fails
     at bus 1, test 3).
   - Reading: failure does not depend on the event time, so the event schedule is not
     the cause. The pre-fault state is essentially steady at these times, so this is
     consistent with both hypotheses and does not discriminate them. It is still a
     self-contained reproducer for GridKit devs.
5. **0.8 bus sweep. Done (job 19794873).** `ida_zoom.py --src
   faultscreen1_ap08/h1879_n26_dc148-200/b1_d0.1 --stage bus --alphas 0.2,0.8,1.0`, all
   200 fault buses. Output `ida_zoom/ap08_bus_h1879/results.csv`.
   - `event` and `closed` modes give the same failing sets. `open` has no failures.
   - Failing buses (`event`): 5 at alphaP 0.2, 21 at 0.8, 23 at 1.0.
   - Nesting: every 0.2 failure is a 0.8 failure, and every 0.8 failure is a 1.0 failure.
     No exceptions. Two buses fail at 1.0 but not at 0.8.
   - Buses 1 and 148: ok at 0.2, IDACalcIC at 0.8 and 1.0.
   - Reading: the 0.8 failures are a subset of the 1.0 failures at this dispatch. This
     supports "the extra 0.8 failures are the 1.0 failures, not new ones", as far as IDA
     IC goes at this one dispatch. The nesting was checked only on h1879; the 0.8 IC
     concordance on two other dispatches is in §10.8.
6. **Mechanism tests on `h1879` (`ida_mechanism.py`, job 19794925, 913 runs, event
   mode). Done.** Output `ida_zoom/mechanism/results.csv`. Fault X is the per-unit
   fault reactance on the 100 MVA base (BusFault README: `B = −X/(R²+X²)`; at the
   reference X = 0.01, R = 0 the fault is a shunt of about −j100 pu).
   - **No datacenter** (`dc_zero`, all 200 buses): only buses 64 and 146 fail
     (IDACalcIC). Caveat: the dispatch is unchanged, so the generators still carry the
     200 MW datacenter setpoint and the slack absorbs the difference.
   - **MW sweep** (failing buses from alphaP 1.0, 23 buses; counts out of 23):

     | DC MW | 25 | 50 | 100 | 150 | 200 |
     |---|---|---|---|---|---|
     | alphaP 0.2 | 2 | 2 | 2 | 2 | 5 |
     | alphaP 0.8 | 2 | 2 | 7 | 16 | 21 |
     | alphaP 1.0 | 2 | 2 | 13 | 17 | 23 |

     At 25–50 MW only buses 64 and 146 fail at every alphaP. Above that, the 0.8 curve
     rises at lower MW than the 0.2 curve.
   - **X sweep at 200 MW** (X = 0.5, 0.3, 0.2, 0.15, 0.1, 0.05, 0.02, 0.01): alphaP 0.2
     fails only at X = 0.01 (5 buses) and X = 0.02 (1 bus). alphaP 0.8 fails at
     X = 0.01 (21 buses) and X = 0.02 (11 buses). Ten buses fail at 0.8 with X = 0.02
     but not at 0.2 (7, 8, 9, 86, 101, 110, 124, 141, 148, 193).
   - **Failures at both 0.2 and 0.8** (the 5 buses at alphaP 0.2: 64, 101, 141, 146, 193):
     - 64 and 146 fail with no datacenter, at every alphaP, from 25 MW. They are not
       driven by P. Their cause is the fault or network (or the IC guess), not the
       datacenter.
     - 101, 141 and 193 pass with no datacenter and fail at 200 MW at alphaP 0.2. So the
       datacenter matters for them. At alphaP 0.8 they already fail at 100 MW (buses 101,
       141 and 193 are in the 0.8 list at 100 MW), so alphaP also lowers their threshold.
   - **Reading (working hypothesis, not verified).** The data are consistent with
     constant-power load shrinking the voltage margin of the fault-on state: a higher
     alphaP makes fails appear at lower MW and milder fault strength. That is what a
     fold, or the absence of a fault-on equilibrium (H1), would produce. Falsifier: a
     fault-on solution that exists at a failing (bus, MW, alphaP) point, with a seed that
     IDA can use. Tested in item 7 below.
7. **Seed test (`ida_seed_test.py`, `ida_zoom/seed_test/`). Done, on the `h1879`
   dispatch, `closed` mode (fault on at t = 0), alphaP irrelevant at dc_zero.** The IDA
   guess at t = 0 is the bus `init` (Vr, Vi) in `case.json` (`Bus::initialize`). For each
   bus, the fault-on AC PF was solved with PowerModels (`pm_solve.jl`, fault shunt from
   `pm_fault_check.build_case`). The solved voltages were written into a copy of
   `case.json`, and the run was repeated with that seed. Output (`IDACalcIC` = IC fail):

   | fault bus | DC | PF fault-on converged | IDA, original guess | IDA, PF seed |
   |---|---|---|---|---|
   | 64 | 0 MW | no | IDACalcIC | not run |
   | 146 | 0 MW | yes | IDACalcIC | IDACalcIC |
   | 1 | 0 MW | yes | ok | ok |
   | 94 | 0 MW | yes | ok | ok |
   | 64 | 200 MW | no | IDACalcIC | not run |
   | 146 | 200 MW | yes | IDACalcIC | IDACalcIC |
   | 1 | 200 MW | no | IDACalcIC | not run |
   | 94 | 200 MW | yes | ok | ok |

   - Bus 146 fails from a converged PF state too. So the guess is not the only cause
     there. Caveat on what the seed changes: the script writes only bus `init` Vr/Vi.
     Machine states are not stored in `case.json`. `Genrou::initialize()`
     (`GenrouImpl.hpp` ~line 378) derives delta, Eqp, Edp and the flux states from the
     seeded Vr/Vi and the pre-fault `p0`/`q0`, and `SystemModel::initialize()` calls it
     at t = 0. So the machine states do follow the seed, but from pre-fault injections.
     The fault-on PF's generator Q (and slack P) are not passed to the machines. The seed
     is therefore not consistent with the PF's generator state, and this does not separate
     H1 from H2 for bus 146.
   - Bus 64 (both DC) and bus 1 (200 MW): the PF has no fault-on solution, and IDA fails.
     Bus 1 (0 MW), 94 (both DC): PF converges, IDA ok. Bus 146: PF converges, IDA fails.
     So PF infeasibility agreed with IDA failure in 3 of 3 cases, and PF feasibility
     agreed with IDA ok in 3 of 5. This is a small sample, and the PF is not validated
     here (the earlier h6483 mismatch at buses 1 and 148 is still unexplained).
   - Next: a seed consistent with the machine states IDA derives. Set `p0`/`q0` to the
     fault-on PF injections (or solve the network algebra with the machine states from
     the pre-fault state held fixed), then rerun IDA. Until then, H1 vs H2 is open for
     bus 146.

## 8. Draft GridKit issue (not filed; superseded)

> **Superseded.** This draft predates the hunt. Its "needs before filing" list is partly
> done (§10–§11). Do not file from it without revising §0 and §11 first.

Title: IDA `IDACalcIC` fails for a constant-power `LoadZIP` under a bolted-like fault
shunt, where a consistent power flow may exist.

The draft in the report (§8) describes an inactive device failing at t = 0. §3a shows the
failing solve is the event-time one with the fault closed, so that draft is superseded.
Needs before filing: tests 1 and 4 above, with the minimal case reproduced from a clean
case directory.

## 9. Terms, with a picture

- **Network algebra.** The power-flow equations F(x) = 0 at every bus: the current
  injections balance. The unknowns are the bus voltages.
- **Fixed-state algebra.** The network algebra with the generator states held at their
  pre-fault values (`gen_emf` in `ida_fixed_state_test.py`, which mirrors the
  event-time re-initialisation). It isolates the network from generator dynamics.
- **Pre-fault state.** The operating point before the fault (fault strength s = 0). The
  hunt cases are balanced to max |F| = 4.07e-5 (`hunt/h1879/cases`, from `ida_hunt_matrix.py`).
- **Fault shunt and fault strength s.** The bolted BusFault is a shunt admittance
  Y_f = 1/(R + jX) at the faulted bus, with R = 0 and X = 0.01 pu. The continuation
  uses s·Y_f (docstring of `ida_fixed_state_test.py`). Effective reactance is X/s, so a
  larger s is a stronger fault. s = 1 is the target fault.
- **Continuation.** Stepping s from 0 to 1 with Newton at each step, starting from the
  previous solution, so it follows one branch: the one connected to the pre-fault state.
  If the step shrinks below a minimum, it stops at **s_stop** (`reached = False`).
- **Upper and lower branch.** For a given fault the algebra can have two nearby roots.
  The upper branch is high-voltage and operable; the lower is low-voltage and
  unstable. IDA starts on the upper branch.
- **Fold, nose, saddle-node.** The point where the upper and lower branches merge. There
  the Jacobian of the algebra is singular, and beyond it no root exists on that branch.
  "Nose" is the common name from P–V curves; "fold" and "saddle-node bifurcation" name
  the same point. A "saddle point" is a different object (an equilibrium of a different
  stability type) and is not what we mean here.
- **X_crit.** X/s_stop: the effective fault reactance at the fold. A bolted fault with
  effective X below X_crit has no upper-branch root. Verified: bus 148 at alphaP 0.8,
  200 MW, X_crit = 0.027236 (`alpha_sweep.csv`).
- **Warm vs flat start.** Warm = Newton starts from the pre-fault solution; flat = starts
  from 1∠0. They can land on different roots (`ida_mw_balanced.py`).
- **ZIP load and alphaP / alphaI / alphaZ.** alphaP is the constant-power share, alphaI
  constant-current, alphaZ constant-impedance; alphaZ = 1 − alphaP − alphaI. The
  campaign's 0.2 case is P/I/Z = 20/80/0; the 0.8 case is 80/20/0.
- **IDACalcIC and IDASolve.** IDACalcIC: IDA could not compute a consistent initial
  condition at the event time. IDASolve: IC found but the integration then failed.
- **Q limits and PV→PQ switching.** A generator at its reactive limit is held at that
  limit and its bus becomes PQ; it switches back to PV when the voltage recovers.
  PowerModels' PF does not do this; `q_limit_pf.py` is an outer loop that does (§10.5).
- **Censored outcome.** A run that has no steady state to report, rather than a random
  failure. A fold gives a censored outcome in this sense.

Schematic only, not to scale. The curve is one branch of V(s); the dashed part is where
the fold removes the root:

```
 V
 1.0 |*
     | *               upper branch (operable): IDA starts at s = 0
     |  *
     |   *
     |    *
 Vn  |     o  <- nose = fold = saddle-node, at s = s_stop  (X_crit = X/s_stop)
     |    *
     |   *             lower branch (unstable): not operable
     |  *
 0.0 |*
     +-------------------------------> s  (fault strength; 0 = no fault)
     0            s_stop     1 (target, X = 0.01 pu)
                  fold       beyond the nose: no upper-branch root (IDACalcIC)
```

In this schematic s = 1 lies beyond s_stop, so there is no upper-branch root at the
target fault. That is the situation the IDACalcIC failures are in. When s_stop > 1 the
continuation reaches the target and the bus is not folded.

## 10. Method and results

### 10.1 Method

- **Fold test.** `ida_fixed_state_test.py` and the `ida_fold_alpha.py`, `ida_zip_sweep.py`
  and `ida_mw_balanced.py` drivers. Fixed-state network algebra, warm start from the
  pre-fault state, continuation in s from 0 to 1, target X = 0.01 pu, R = 0. Output:
  `reached`, `s_stop`, `X_crit = X/s_stop`. A **fold** is `reached = False`.
- **IDA test.** `ida_zoom.py --stage bus`, event mode, one run per bus, same campaign
  case (`b1_d0.1/case.json`). Status per bus: ok, IDACalcIC or IDASolve.
- **Concordance.** Same bus, same alphaP, same case: does the fold set equal the set of
  IDACalcIC failures?
- **Case.** ACTIVSg200, reference dispatch `h1879_n26_dc148-200` (datacenter at bus 148,
  200 MW).   Balanced hunt cases from `ida_hunt_matrix.py` and `ida_mw_balanced.py` (pre-fault
  max |F| = 4.07e-5).- **Fixed parameters.** Fault X = 0.01 pu, R = 0 (default; not varied, per request).
  alphaP = 0.8 is the pilot level; 1.0 is shown as the upper bound only.

### 10.2 alphaP at 200 MW (fixed-state folds vs IDA)

| alphaP | folds (fixed-state) | IDACalcIC | IDASolve | ok | fold = IDACalcIC? |
|---|---|---|---|---|---|
| 0.2 | 3 (64, 141, 146) | 3 (64, 141, 146) | 2 (101, 193) | 195 | exact |
| 0.8 | 18 | 18 | 3 (2, 134, 142) | 179 | exact |
| 1.0 | 21 | 21 (incl. 2, 134, 142) | 2 (116, 131) | 177 | exact |

- Sources: `fixed_state/alpha_sweep.csv`; `ap08_bus_h1879/results.csv` (event mode).
- X_crit over folded buses: 0.2: 0.0101 (bus 141) to 0.0210 (bus 64). 0.8: 0.0103
  (bus 117) to 0.0293 (bus 141). 1.0: 0.0113 (bus 134) to 0.0359 (bus 141).
- The fold sets are nested as alphaP rises (3 ⊂ 18 ⊂ 21).
- At 0.8 the three IDASolve buses (2, 134, 142) are **not folds at 200 MW**. At 250 MW
  they are: X_crit 0.01905 (bus 2), 0.01418 (bus 134), 0.01891 (bus 142). They are near
  a fold, not on one at 200 MW.
- At 0.2 the IDASolve buses 101 and 193 are not folds at 200 MW.

### 10.3 alphaZ (constant-impedance share), 200 MW

| alphaZ | 0 | 0.25 | 0.5 | 0.75 |
|---|---|---|---|---|
| folds | 21 | 16 | 10 | 2 (64, 146) |

- alphaZ = 0 is alphaP = 1.0 and gives 21, as in §10.2: consistent.
- At alphaZ 0.5 the folds are 7, 8, 64, 86, 101, 110, 141, 146, 148, 193. Nested in
  alphaZ.
- Source: `fixed_state/zip_sweep.csv`; log `zip_sweep_slurm.log` (job 19795920).

### 10.4 Megawatts (datacenter size), alphaP 0.8, unlimited Q

| MW | 100 | 150 | 200 | 250 | 300 | 350 | 400 |
|---|---|---|---|---|---|---|---|
| folds | 7 | 14 | 18 | 23 | 25 | 30 | 39 |

- Nested across MW. The grid is 50 MW.
- Bus 146 folds at all MW (X_crit 0.011567 at 100 MW to 0.012003 at 400 MW); bus 64 at
  all MW (0.02095 to 0.02135).
- Source: `mw_balanced/fold_mw.csv` (job 19795961).

### 10.5 Q limits (PV→PQ outer loop, `q_limit_pf.py`)

| MW | 100 | 150 | 200 | 250 | 300 | 350 | 400 |
|---|---|---|---|---|---|---|---|
| folds, unlimited Q | 7 | 14 | 18 | 23 | 25 | 30 | 39 |
| folds, Q-limited | 7 | 14 | 18 | 23 | 26 | 33 | 42 |
| buses added by Q limits | — | — | — | — | 29 | 128, 135, 140 | 113, 187, 198 |

- The unlimited set is a subset of the Q-limited set at every MW (checked from the CSVs).
- At 200 MW one unit is switched to PQ; the fold set is unchanged.
- At 400 MW twelve units are switched to PQ (converged, balance |F| about 6e-5).
- Source: `hunt/h1879/qlim_grid.csv` (job 19797220); recount `q_limit/fold_recount/`
  (job 19796511).

### 10.6 Dispatch robustness, 200 MW, alphaP 0.8, unlimited Q

| dispatch | folds | vs h1879 |
|---|---|---|
| h1879 (reference) | 18 | — |
| h2151 | 15 | subset; missing 64, 117, 146 |
| h4770 | 16 | subset; missing 117, 146 |

- At 400 MW: h2151 36 (unlimited) / 37 (Q-limited); h4770 38 / 39; h1879 39 / 42.
- Bus 146 is **not** a fold on h2151 or h4770, so "bus 146 folds in every case" is wrong.
  It folds on h1879 in every mode, including constant-Z.
- Bus 64 folds on h1879 and h4770, not on h2151.
- Sources: `hunt/h2151/h2151.csv`, `hunt/h4770/h4770.csv` (jobs 19796643, 19796644).

### 10.7 Constant-impedance control, h1879

- alphaZ = 1 (alphaP = 0, the ZIP is pure impedance at the datacenter).
- Folds at 200 and 400 MW, unlimited or Q-limited: buses 64 and 146 only. Bus 148 does
  not fold. Consistent with the alphaZ sweep at 0.75 (§10.3).
- Source: `hunt/h1879/zconst.csv` (job 19796642).
- Caveat: the generators keep their 200 MW setpoint in this control, so it is not a
  pure "no datacenter" case.

### 10.8 Cross-dispatch IDA check

- Job 19797339, `ida_xdisp.sh`: IDA bus sweep at alphaP 0.8, 200 MW, event mode, for
  h1879 (control), h2151 and h4770 (balanced hunt cases `hunt/<d>/cases/unlim_mw200.json`).
  All 200 buses per dispatch; closed and open modes were also run by the stage and are
  not used here.
- Prediction under H1: the IDACalcIC set equals the fold set of the same dispatch
  (§10.6).
- Outputs: `ida_zoom/xdisp_{h1879,h2151,h4770}/results.csv`; SLURM log
  `ida_zoom/xdisp/slurm_19797339.log`.

| dispatch | IDACalcIC | fold set (200 MW, unlimited Q) | equal? | IDASolve | ok |
|---|---|---|---|---|---|
| h1879 (control) | 18 | 18 | exact | 2, 134, 142 | 179 |
| h2151 | 15 | 15 | exact | 2, 117, 134, 142 | 181 |
| h4770 | 16 | 16 | exact | 2, 117, 134, 142 | 180 |

- **IC concordance holds on all three dispatches.** Every IDACalcIC bus is a fold, and
  every fold is an IDACalcIC bus. Prediction P5 holds for the IC failures. Falsifier F3
  (IC failure at a fold-free bus) was not observed.
- The h1879 control reproduces the earlier campaign-case sweep (§10.2) exactly: the
  balanced hunt case and the campaign case give the same IC and IDASolve buses at 0.8.
- **Caveat: IDASolve at a fold-free bus.** Bus 117 is an IDASolve failure on h2151 and
  h4770, but both reach s = 1 there, so the fold-free algebra has a root at the target
  fault. This failure is an integration failure at a bus with no fold on that dispatch.
  H1 does not explain it. It is one of the "couple of IDASolve fails" the user asked not
  to chase (§0 item ii). Bus 117 is a fold on h1879 (X_crit 0.010314, the 0.8 minimum).
  On h2151 and h4770 it is not a fold at all, so it is not a near-fold in the same sense.

### 10.9 Caveats on these results

- X = 0.01 pu bolted fault only. Finite-impedance faults would move every X_crit.
- Warm-start roots only. A fold means no root on the continuation branch; it does not
  rule out a root on another branch (see bus 146 in §11.3).
- MW grid is 50 MW; alphaP grid is 0.2, 0.8, 1.0; alphaZ grid is 0.25.
- Exciter and field limiters (SexsPti) and OEL/UEL are not modelled in the fold tests.
- Q limit tolerance 0.1 MVAr and a 30-iteration cap are assumed harmless (not varied).

## 11. Hypothesis H1 and what it means

### 11.1 Statement

H1: at alphaP 0.8 (and 1.0), a constant-power datacenter at bus 148 moves the fault-on
fold of the network algebra so that a bolted fault at each folded bus lies beyond the
nose. At those buses no operable post-fault root exists on the branch reached from the
pre-fault state. IDA's IC failure is then a physical no-steady-state outcome, not an IDA
defect.

### 11.2 Predictions and status

| # | Prediction | Test | Status |
|---|---|---|---|
| P1 | IDA IC failures = folds, per bus | §10.2 (3 alphaP levels) | Supported at h1879 200 MW (exact at 0.2, 0.8, 1.0) |
| P2 | Folds grow with alphaP (nested) | §10.2 | Supported (3 ⊂ 18 ⊂ 21) |
| P3 | Folds shrink as load becomes constant-Z | §10.3 | Supported (21, 16, 10, 2) |
| P4 | Folds grow with MW (nested) | §10.4 | Supported (7 to 39) |
| P5 | Fold sets are dispatch-dependent but of the same type | §10.6, §10.8 | Supported: IC = fold exactly on h1879, h2151, h4770 |
| P6 | Q limits can only add folds | §10.5 | Supported at all MW (superset) |
| F1 | Falsifier: a fold at a bus with IDA IC ok at the same alphaP and MW | Concordance | Not observed at 200 MW |
| F2 | Falsifier: an operable upper-branch root at an IDA-failing bus | Seed test (bus 146) | Open (§11.3) |
| F3 | Falsifier: IDA IC fails at a fold-free bus, on another dispatch | §10.8 | Not observed (IC only; IDASolve bus 117 is a separate case, §10.8) |

### 11.3 Alternatives

- **H2, IDA initial guess.** Predicts failures where a fault-on root exists. Against it:
  the exact fold match at all three alphaP. For it: the bus 146 seed test. PF converges
  there, yet the continuation folds. The PF root may be on another branch. Unresolved
  (§0). H2 is not excluded for bus 146 or the three IDASolve buses at 0.8.
- **H3, solver settings.** Not tested by the fold scans, which run no IDA. The fold
  result is therefore independent of IDA tolerances. IDA-side settings are a separate
  question and were not varied in this hunt.
- **H4, PF vs IDA mismatch at buses 1 and 148 (§4).** Unresolved for the h6483
  dispatch. At 200 MW on h1879 (alphaP 0.8), buses 1 and 148 are both in the fold set,
  so the fold reading is consistent with the IDA result there. The fold set at h6483 was
  not computed, so this is not verified for that dispatch.

### 11.4 What would falsify H1

- A fold at a bus where IDA IC is ok at the same case, alphaP and MW.
- An operable upper-branch root (for example from a different continuation path, such as
  MW continuation or a different seed) at a bus where IDA IC fails.
- Cross-dispatch disagreement in §10.8 (IDACalcIC set differs from the fold set). None
  observed on the three dispatches tested.

### 11.5 What it means

1. **Pilot campaign at 0.8 and 1.0.** At the folded buses the IDA IC failure is a
   physical no-steady-state outcome. Those runs are censored (§0), not solver noise to
   repair. The censoring guidance in the pilot report (§7.5) applies.
2. **Siting.** A higher constant-power share leaves fewer buses where a bolted fault
   stays on the operable branch: 3 → 18 → 21 folds across alphaP 0.2 → 0.8 → 1.0 at
   200 MW. In practice the datacenter can be placed at fewer buses where the modelled
   fault is simulable. This is the "fewer places" consequence the user anticipated.
3. **High-P modelling.** alphaP is the lever; 0.8 is the pilot level. At 250 MW and
   above, Q limits change the fold set (§10.5), so a 0.8 campaign at those MW should
   enforce Q limits (or report X_crit).
4. **Report X_crit, not counts.** The fold count depends on the fixed X = 0.01 pu. The
   per-bus X_crit is the more useful output and can be compared across fault
   definitions.
5. **Not a GridKit defect (under H1).** A clearer message at IC time ("no fault-on
   equilibrium on this branch") would make the censored outcome visible. That is a
   feature request, not a bug, and it is optional.
6. **Industry practice.** Planning PFs enforce generator Q limits by PV→PQ switching.
   Dynamic studies rely on exciter and limiter models (OEL/UEL) and do not switch buses
   in the PF. That is why the Q-limited fold recount matters for siting, and why the
   exciter and limiter models are a caveat (§10.9).

### 11.6 What is not claimed

- Proof for all dispatches, all X and R, all load mixes and all machine states.
- That a fold rules out every root. It rules out roots on the branch reached by
  continuation from the pre-fault state.
- That the three IDASolve buses at 0.8 are explained at 200 MW. They are near-folds at
  250 MW.

## 12. Experiment index

Scripts are in `/home/isatkaus/gridkit/uq-usecase/scripts/`. Outputs are under
`/scratch/isatkaus/gridkit-campaigns/ida_zoom/` unless stated. Python:
`/home/isatkaus/flash-w-dom/conda-envs/h-py312-basic/bin/python`. Julia via
`pm_utils.resolve_julia()` (1.12.1). Heavy runs go to SLURM (account `scidac`, partition
`short`). Login-node parallel runs need `OMP_NUM_THREADS=1` (and the BLAS/MKL variants),
or they oversubscribe.

| # | Test | Script(s) | Output | SLURM job | Result |
|---|---|---|---|---|---|
| 1 | Fault-strength sweep | `ida_xsweep.py` | `ida_zoom/xsweep181/results.csv` | 19794871 | §7 item 3 |
| 2 | Event-time reproducer | `ida_event_repro.py` | `ida_zoom/event_repro/results.csv` | 19794872 | §7 item 4 |
| 3 | 0.8 bus sweep (campaign case) | `ida_zoom.py --stage bus` | `ida_zoom/ap08_bus_h1879/results.csv` | 19794873 | §10.2, §5 |
| 4 | Mechanism tests | `ida_mechanism.py` | `ida_zoom/mechanism/results.csv` | 19794925 | §7 item 6 |
| 5 | Seed test (PF seed vs continuation) | `ida_seed_test.py` | `ida_zoom/seed_test/` | not recorded | §7 item 7, §11.3 |
| 6 | PV→PQ PF (PowerModels) | `ida_pvpq.py`, `ida_pvpq.jl` | `ida_zoom/pvpq/results.csv` | not recorded | §7 item 2, §4 |
| 7 | Fold vs alphaP | `ida_fold_alpha.py`, `ida_fold_alpha.sh` | `fixed_state/alpha_sweep.csv`; log `fold_alpha_slurm.log` | not recorded | §10.2 |
| 8 | Fold vs alphaZ | `ida_zip_sweep.py`, `ida_zip_sweep.sh` | `fixed_state/zip_sweep.csv`; log `zip_sweep_slurm.log` | 19795920 | §10.3 |
| 9 | Fold vs MW (balanced cases) | `ida_mw_balanced.py`, `ida_mw_balanced.sh` | `mw_balanced/fold_mw.csv` | 19795961 | §10.4 |
| 10 | Q-limit outer loop | `q_limit_pf.py` (module) | `q_limit/mw200/`, `q_limit/mw400/` | inline | §10.5 |
| 11 | Q-limited fold recount | `ida_q_fold_recount.py`, `ida_q_fold_recount.sh` | `q_limit/fold_recount/fold_qlim_mw.csv` | 19796511 | §10.5 |
| 12 | Hunt: h1879 Q-limit grid | `ida_hunt_matrix.py`, `ida_hunt_matrix.sh --dispatch h1879 --modes qlim` | `hunt/h1879/qlim_grid.csv` | 19797220 (19796641 crashed, fixed) | §10.5 |
| 13 | Hunt: constant-Z | `ida_hunt_matrix.sh --loads zconst` | `hunt/h1879/zconst.csv` | 19796642 | §10.7 |
| 14 | Hunt: h2151 dispatch | `ida_hunt_matrix.sh --dispatch h2151` | `hunt/h2151/h2151.csv` | 19796643 | §10.6 |
| 15 | Hunt: h4770 dispatch | `ida_hunt_matrix.sh --dispatch h4770` | `hunt/h4770/h4770.csv` | 19796644 | §10.6 |
| 16 | Hunt smoke test | `ida_hunt_matrix.py` (login node) | `hunt/h1879/smoke.csv` | none | §10.1 |
| 17 | Cross-dispatch IDA check | `ida_xdisp.sh` (calls `ida_zoom.py --stage bus`) | `xdisp/`, `xdisp_{h1879,h2151,h4770}/results.csv`; log `xdisp/slurm_19797339.log` | 19797339 | §10.8 (done) |

Earlier analysis scripts (§2–§5): `ida_discriminate.py`, `ida_closed_check.py`,
`ida_closed_variants.py`, `ida_mw_balanced.py` (also row 9), `ida_fixed_state_test.py`
(core algebra and continuation; used by rows 7–9 and the hunt).

**Reproduce a row.** Run `sbatch scripts/<script>.sh` for rows with a wrapper, or
`python scripts/<script>.py` on the login node for small runs. For the hunt:
`sbatch scripts/ida_hunt_matrix.sh --dispatch h1879 --modes unlim qlim --mw 200 --tag <name>`.
Check the output CSVs listed above against the tables in §10.

**Caveat on job IDs.** The job IDs above come from session notes. They do not appear in
the log files under `ida_zoom/`, so confirm each one with `sacct` before the report is
finalised. The rows marked "not recorded" have no ID in the notes.
