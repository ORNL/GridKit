# Work Notes 4 — Epistemic UQ runs: full inventory and report metrics (2026-09-11)

Companion to [`notes_3.md`](notes_3.md) (Hawaii v6–v10 hypothesis designs). This note
compiles **every** epistemic run completed to date into one table, and the headline
metrics (data volume, sample count, throughput) needed for the report.

All figures below are read directly from the run store at
`/kfs2/projects/scidac/scidac-data/gridkit-runs/` (per-run `meta.yml`, `du -sh`, and
SLURM `.out`/`.err` file timestamps), verified in this session on 2026-09-11.

---

## Summary comparison table — all epistemic runs (v1 → v10 + Morris)

| Case | Run | N samples | Distribution | Gens sampled | Purpose / hypothesis | Disk size |
|---|---|---:|---|---:|---|---:|
| Hawaii40 | v1 | 20 | uniform ±10% | 4 | initial bus-diverse pilot | 267M |
| Hawaii40 | v2 | 1,000 | uniform ±10% | 4 | bus-diverse, scale-up test | 15G |
| Hawaii40 | v3 | 1,000 | normal σ=12% | 4 | switch to Gaussian marginal | 12G |
| Hawaii40 | v4 | 4,000 | normal σ=12% | 2 | narrowed gen set | 46G |
| Hawaii40 | **v5** | **16,000** | normal σ=12% | 4 | production epistemic (hop distances 1–3+) | **182G** |
| Hawaii40 | v6 | 4,000 | uniform ±20% | 4 | smallest H (weakest intrinsic response) | 46G |
| Hawaii40 | v7 | 4,000 | uniform ±20% | 4 | farthest from fault (electrical periphery) | 46G |
| Hawaii40 | v8 | 4,000 | uniform ±20% | 4 | smallest mva + dispatch (weak coupling) | 46G |
| Hawaii40 | v9 | 4,000 | uniform ±20% | 4 | largest H×mva (positive control) | 46G |
| Hawaii40 | v10 | 4,000 | uniform ±20% | 4 | small p0 / lightly dispatched | 46G |
| Hawaii40 | ian-csv | 1,000 | Morris (external design) | 39 (full fleet) | fleet-wide sensitivity screening | 12G |
| Illinois (ACTIVSg200) | v1 | 1,000 | normal σ=12% | 4 | initial large-case pilot | 36G |
| Illinois (ACTIVSg200) | **v2** | **4,000** | normal σ=12% | 4 | production epistemic (hop distances 4–9) | **144G** |

Not counted toward the sample/data totals below: `threebusbasic-uq/` and
`threebusbasic-v1/` (77M each) — small early dev/smoke-test runs on the 3-bus toy case,
no parquet results, not part of the epistemic UQ campaign proper.

See [`cases/hawaii.md`](../cases/hawaii.md) and [`cases/illinois.md`](../cases/illinois.md)
for per-run generator selection and version-history detail; [`notes_3.md`](notes_3.md)
for the v6–v10 hypothesis rationale.

---

## Headline metrics for the report

- **Total dynamic simulations run to date: 48,020** (sum of the N-samples column above).
- **Total data volume: ~678 GB** (sum of the disk-size column above; dominated by
  hawaii-v5 at 182G and illinois-v2 at 144G, which together are ~48% of the total).
- **Two grids covered**: Hawaii40 (37 bus) and Illinois/ACTIVSg200 (200 bus).
- Per-run data lands as one Parquet file per sample (`run_NNN.parquet`) — 48,020 Parquet
  files currently on disk (one-to-one with the sample count above; verified by
  `find <run>/runs -name '*.parquet' | wc -l` matching `meta.yml` `n_samples` for every
  run listed).

## Bragging points (for the report)

- **48,020 dynamic simulations, ~678 GB of results**, generated to date across two grids
  (Hawaii40, Illinois/ACTIVSg200).
- **200-bus Illinois grid**: 4,000 full nonlinear dynamic simulations completed in
  under 2 minutes of compute per 1,000-sample chunk, across 4 nodes × 4 workers
  (16-way parallel).
- **Our interactive `m_viz` workflow for dynamic-case visualization on the geographic
  grid layout already scales to 10,000-bus grids**, using scalable binning to keep
  rendering fast.
- **Next FY**: scaling this same pipeline to interconnect-size grids (10K–70K buses —
  WECC, Eastern Interconnect), which will require considerably more than 4 nodes per
  sweep.

## Throughput: parallel SLURM pipeline timing

Verified fact: figures below are from `sacct -j <jobid> --format=JobID,JobName,Elapsed,...`
against Kestrel's job-accounting history (queried directly by job ID; the history goes
back well past these July 2026 jobs). Each chunk runs on 1 SLURM node/task, fanning out
**4 parallel workers via `xargs -P 4`** inside the node; 4 chunks submitted per run =
**16-way parallelism** across 4 nodes.

**Hawaii40 v5 — 16,000 samples** (job IDs 15345297–15345301):

| Job | Role | Elapsed | Start | End |
|---|---|---|---|---|
| 15345297–300 | 4× simulation chunks (4,000 samples/node, 4 workers/node) | 00:05:10–00:05:15 | 10:50:25 | 10:55:35–10:55:40 |
| 15345301 | collect (parquet write, `collect.sh`) | 00:06:16 | 10:55:48 | 11:02:04 |
| **End-to-end total** | | | 10:50:25 | 11:02:04 (**~11.6 min**) |

**16,000 dynamic simulations completed in ~5 min 15 sec of SLURM `Elapsed` compute time**
across 4 nodes × 4 workers (16-way parallel) — ≈ 3,050 samples/minute aggregate, ≈ 190
samples/minute per worker. Collection into 16,000 individual Parquet files added another
~6 min 16 sec. Full run (submit → all data collected): **~11.6 minutes**.

**Illinois v2 — 4,000 samples** (200-bus case; job IDs 15338688–91, 15338753–56):

| Job | Role | Elapsed | State | Start | End |
|---|---|---|---|---|---|
| 15338688 | chunk 0, attempt 1 | 00:00:01 | **FAILED** | 15:40:09 | 15:40:10 |
| 15338689 | chunk 1 | 00:01:41 | COMPLETED | 15:40:12 | 15:41:53 |
| 15338690 | chunk 2 | 00:01:40 | COMPLETED | 15:41:31 | 15:43:11 |
| 15338691 | chunk 3 | 00:01:43 | COMPLETED | 15:41:51 | 15:43:34 |
| 15338753 | chunk 0, attempt 2 (resubmitted) | 00:01:47 | COMPLETED | 15:57:21 | 15:59:08 |

No `collect.sh` SLURM job exists for illinois-v2 (`COLLECT_IN_SLURM=False`; collected
manually from the notebook instead).

Each successful 1,000-sample chunk (4 workers/node) completed in **~1 min 40 sec –
1 min 47 sec** of `Elapsed` compute — i.e. once running, a 200-bus, 1,000-sample chunk is
about **20x faster in wall-clock than the raw sample count would suggest naively**
(≈ 590 samples/minute aggregate per chunk, ≈ 145 samples/minute per worker). The one
chunk-0 failure (1 second, immediate) was resubmitted independently ~17 minutes later
(unrelated queue/scheduling gap, not a measure of compute time); the 4 chunks' actual
compute windows sum to under 7 minutes total for all 4,000 samples.

**Bottom line for the report**: the HPC-parallelized pipeline (LHS sampling → per-sample
case patching → SLURM node/worker fan-out → parallel Parquet collection) runs a
16,000-sample production epistemic sweep (Hawaii40) in **under 12 minutes end-to-end**,
and a 4,000-sample sweep on the much larger 200-bus Illinois case in **under 2 minutes of
compute per 1,000-sample chunk**. The entire 48,020-sample, ~678 GB campaign to date
represents cumulative compute that would have taken far longer run serially on a single
core.

---

## Caveats / what would sharpen these numbers further

1. **Illinois v2 chunk 0 resubmission** — the first attempt failed in 1 second (job
   15338688); root cause not re-diagnosed here. The resubmitted attempt (15338753)
   completed normally in 1:47, in line with the other 3 chunks, so it does not affect
   the per-chunk throughput figure, only the total wall-clock-from-first-submit window.
2. Disk sizes are `du -sh` on live Lustre storage; figures can drift slightly if any
   collection job is still appending files. All numbers above were captured in one
   session (2026-09-11) so they are mutually consistent.
3. If exact core-hour/node-hour billing numbers are needed for the report (rather than
   wall-clock elapsed), `sacct --format=...,AllocCPUS,CPUTime` or `sreport` can pull that
   directly using the same job IDs listed above.
3. Disk sizes are `du -sh` on live Lustre storage; figures can drift slightly if any
   collection job is still appending files. All numbers above were captured in one
   session (2026-09-11) so they are mutually consistent.
