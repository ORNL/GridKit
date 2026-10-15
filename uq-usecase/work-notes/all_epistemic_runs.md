# All Epistemic Runs — Inventory & Metrics

**Live table of all completed epistemic UQ sweeps.** Update this table whenever a new
epistemic run completes. Reference for sample counts, disk usage, and throughput metrics.

Data verified against run store at `/kfs2/projects/scidac/scidac-data/gridkit-runs/`
(per-run `meta.yml`, `du -sh`, and SLURM accounting).

---

## Summary comparison table — all epistemic runs

| Case | Run | N samples | Distribution | Gens sampled | Purpose / hypothesis | Disk size |
|---|---|---:|---|---:|---|---:|
| Hawaii40 | hawaii-v1 | 20 | uniform ±10% | 4 | initial bus-diverse pilot | 267M |
| Hawaii40 | hawaii-v2 | 1,000 | uniform ±10% | 4 | bus-diverse, scale-up test | 15G |
| Hawaii40 | hawaii-v3 | 1,000 | normal σ=12% | 4 | switch to Gaussian marginal | 12G |
| Hawaii40 | hawaii-v4 | 4,000 | normal σ=12% | 2 | narrowed gen set | 46G |
| Hawaii40 | **hawaii-v5** | **16,000** | normal σ=12% | 4 | production epistemic (hop distances 1–3+) | **182G** |
| Hawaii40 | hawaii-v6 | 4,000 | uniform ±20% | 4 | smallest H (weakest intrinsic response) | 46G |
| Hawaii40 | hawaii-v7 | 4,000 | uniform ±20% | 4 | farthest from fault (electrical periphery) | 46G |
| Hawaii40 | hawaii-v8 | 4,000 | uniform ±20% | 4 | smallest mva + dispatch (weak coupling) | 46G |
| Hawaii40 | hawaii-v9 | 4,000 | uniform ±20% | 4 | largest H×mva (positive control) | 46G |
| Hawaii40 | hawaii-v10 | 4,000 | uniform ±20% | 4 | small p0 / lightly dispatched | 46G |
| Hawaii40 | hawaii-ian-csv | 1,000 | Morris (external design) | 39 (full fleet) | fleet-wide sensitivity screening | 12G |
| Illinois (ACTIVSg200) | illinois-ian-csv | 1,025 | Morris (external design) | 39 (full fleet) | fleet-wide sensitivity screening (200-bus) | 37G |
| Illinois (ACTIVSg200) | illinois-v1 | 1,000 | normal σ=12% | 4 | initial large-case pilot | 36G |
| Illinois (ACTIVSg200) | **illinois-v2** | **4,000** | normal σ=12% | 4 | production epistemic (hop distances 4–9) | **144G** |
| Illinois (ACTIVSg200) | illinois-v6-il | 4,000 | uniform ±20% | 4 | smallest H (weakest intrinsic response) | 144G |
| Illinois (ACTIVSg200) | illinois-v7-il | 4,000 | uniform ±20% | 4 | farthest from fault (electrical periphery) | 144G |
| Illinois (ACTIVSg200) | illinois-v8-il | 4,000 | uniform ±20% | 4 | smallest mva + dispatch (weak coupling) | 144G |
| Illinois (ACTIVSg200) | illinois-v9-il | 4,000 | uniform ±20% | 4 | largest H×mva (positive control) | 144G |
| Illinois (ACTIVSg200) | illinois-v10-il | 4,000 | uniform ±20% | 4 | small p0 / lightly dispatched | 144G |

**Not included in totals below:** `threebusbasic-uq/` and `threebusbasic-v1/` (77M each)
— small early dev/smoke-test runs on the 3-bus toy case, no parquet results.

**Related documentation:**
- Per-run generator selection and version history: [`cases/hawaii.md`](../cases/hawaii.md), [`cases/illinois.md`](../cases/illinois.md)
- Hypothesis rationale for Hawaii v6–v10: [`notes_3.md`](notes_3.md)

---

## Headline metrics

- **Total dynamic simulations: 69,045** (sum of N-samples column above)
- **Total data volume: ~1.4 TB** (sum of disk-size column; Hawaii-v5 + Illinois-v2 + v6-v10-il = ~65% of total)
- **Grids covered: 2** (Hawaii40: 37 bus; Illinois/ACTIVSg200: 200 bus)
- **Data format:** One Parquet file per sample (~69K files on disk)

---

## Throughput & HPC Performance

### Hawaii40 v5 (16,000 samples, 4 nodes × 4 workers = 16-way parallelism)

| Job | Role | Elapsed | Start | End |
|---|---|---|---|---|
| 15345297–300 | 4× chunks (4,000 samples/node) | 00:05:10–00:05:15 | 10:50:25 | 10:55:35–10:55:40 |
| 15345301 | collect (Parquet write) | 00:06:16 | 10:55:48 | 11:02:04 |
| **End-to-end** | | | **10:50:25** | **11:02:04 (~11.6 min)** |

**16,000 simulations in ~5.2 min compute** (3,050 samples/min aggregate, 190/min per worker).
Collection added ~6.3 min. Full run: **~11.6 minutes**.

### Illinois v2 (4,000 samples, 200-bus case, 4 chunks)

| Job | Role | Elapsed | State | Start | End |
|---|---|---|---|---|---|
| 15338688 | chunk 0, attempt 1 | 00:00:01 | **FAILED** | 15:40:09 | 15:40:10 |
| 15338689 | chunk 1 | 00:01:41 | COMPLETED | 15:40:12 | 15:41:53 |
| 15338690 | chunk 2 | 00:01:40 | COMPLETED | 15:41:31 | 15:43:11 |
| 15338691 | chunk 3 | 00:01:43 | COMPLETED | 15:41:51 | 15:43:34 |
| 15338753 | chunk 0, resubmit | 00:01:47 | COMPLETED | 15:57:21 | 15:59:08 |

Each chunk (1,000 samples, 4 workers): **~1 min 40–47 sec** (~590 samples/min aggregate,
145/min per worker). Chunk-0 first attempt failed (scheduling), resubmitted successfully.

---

## How to update this table

1. Run an epistemic sweep via `uq_setup.ipynb`
2. Once all simulations + collection complete, note:
   - Run ID (e.g., `v11`, `v2-retry`, `morris-design-2`)
   - N_SAMPLES from `meta.yml`
   - Distribution method (uniform, normal, Morris, etc.)
   - Number of generators sampled
   - Purpose/hypothesis (one sentence)
   - Disk size: `du -sh <UQ_RUN_ROOT>`
3. Add a row to the table above
4. Commit this file to git

---

## Appendix: Bragging points (for reports/proposals)

- **49,045 dynamic simulations, ~715 GB of results** across Hawaii40 and Illinois grids.
- **200-bus Illinois grid:** 4,000 nonlinear simulations in under 2 minutes compute
  per 1,000-sample chunk (4 nodes × 4 workers); 1,025-sample Morris design in ~37 min total.
- **Interactive visualization at 10K-bus scale:** `m_viz` workflow already scales with
  geographic binning.
- **Next FY target:** Interconnect-size grids (10K–70K buses: WECC, Eastern Interconnect).

---

**Last updated:** 2026-09-28
