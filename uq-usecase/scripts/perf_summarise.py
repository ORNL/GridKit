"""Summarise a perf_sample.sh CSV into CPU, memory and wall-time numbers.

    python perf_summarise.py <perf_samples.csv> [--workers N] [--cores C] [--out summary.json]

A worker is a python process whose parent is not python (the sbatch shell).
Julia children are attributed to the worker that spawned them. Peaks come from
the sampled timeline, so a process that lives less than one sample interval is
missed. CPU totals use the last cumulative cputimes value seen per pid.
"""

import argparse
import json
from pathlib import Path

import pandas as pd


def summarise(df: pd.DataFrame, workers: int | None, cores: int | None) -> dict:
    df = df.copy()
    df["rss_gb"] = df["rss_kb"] / 1024**2

    t0, t1 = df["t"].min(), df["t"].max()
    wall_s = int(t1 - t0)

    # Worker roots: python processes whose parent is not a python process.
    py = df[df["comm"].str.startswith("python")]
    py_pids = set(py["pid"])
    roots = set(py.loc[~py["ppid"].isin(py_pids), "pid"])
    tree_pid = {}
    for pid in roots:
        tree_pid[pid] = pid
    for _, r in df[df["comm"].str.startswith("julia")].iterrows():
        if r["ppid"] in roots:
            tree_pid[r["pid"]] = r["ppid"]
    df["worker"] = df["pid"].map(tree_pid)

    per_t = df.groupby("t")["rss_gb"].sum()
    per_worker_t = df.dropna(subset=["worker"]).groupby(["t", "worker"])["rss_gb"].sum()
    per_worker_peak = per_worker_t.groupby("worker").max()

    last_cpu = df.groupby("pid")["cpu_s"].max()
    cpu_by_worker = df.dropna(subset=["worker"]).groupby("worker")["pid"].apply(
        lambda p: float(last_cpu.loc[p.unique()].sum())
    )

    total_cpu_s = float(last_cpu.sum())
    n_workers = int(len(roots))
    out = {
        "samples": int(df["t"].nunique()),
        "wall_s": wall_s,
        "workers_seen": n_workers,
        "workers_requested": workers,
        "peak_node_rss_gb": round(float(per_t.max()), 2),
        "peak_worker_tree_rss_gb": round(float(per_worker_peak.max()), 3),
        "median_worker_tree_peak_rss_gb": round(float(per_worker_peak.median()), 3),
        "julia_processes_seen": int(df.loc[df["comm"].str.startswith("julia"), "pid"].nunique()),
        "total_cpu_s": round(total_cpu_s, 1),
        "mean_cpu_per_worker_s": round(float(cpu_by_worker.mean()), 1) if n_workers else None,
    }
    if cores:
        out["allocated_cores"] = cores
        out["cpu_utilisation_of_allocation"] = round(total_cpu_s / (wall_s * cores), 3) if wall_s else None
    if n_workers and wall_s:
        out["cpu_utilisation_per_worker"] = round(float(cpu_by_worker.mean()) / wall_s, 3)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--workers", type=int)
    ap.add_argument("--cores", type=int)
    ap.add_argument("--out")
    a = ap.parse_args()
    df = pd.read_csv(a.csv)
    s = summarise(df, a.workers, a.cores)
    print(json.dumps(s, indent=2))
    if a.out:
        Path(a.out).write_text(json.dumps(s, indent=2) + "\n")


if __name__ == "__main__":
    main()
