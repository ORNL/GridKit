"""Step statistics of line_distributed runs without NumPy (runs inside the container).

    python3 steps_stats.py <window_start_s> <window_end_s> <run dir> [<run dir> ...] [--trace <run dir>]
"""

import csv
import json
import statistics
import sys


def stats(path, lo, hi):
    data = json.load(open(path + "/input.json"))
    tau = min(d for d in data["line"]["H"]["delays"] if d > 0)
    rows = [(float(r["t"]), float(r["h"])) for r in csv.DictReader(open(path + "/steps.csv"))]
    window = [h for t, h in rows if lo <= t <= hi]
    above = sum(h for h in window if h > tau) / sum(window)
    return tau, statistics.median(window), max(window), above, len(rows)


def main(argv):
    lo, hi = float(argv[0]), float(argv[1])
    trace = None
    paths = []
    rest = argv[2:]
    while rest:
        item = rest.pop(0)
        if item == "--trace":
            trace = rest.pop(0)
        else:
            paths.append(item)
    print(f"{'run':34s} {'tau us':>8s} {'median h us':>12s} {'h/tau':>7s} {'max h/tau':>10s} {'time>tau %':>11s} {'steps':>7s}")
    for path in paths:
        tau, med, mx, above, n = stats(path, lo, hi)
        print(f"{path.split('/')[-1]:34s} {tau * 1e6:8.1f} {med * 1e6:12.1f} {med / tau:7.2f} {mx / tau:10.2f} {100 * above:11.1f} {n:7d}")
    if trace:
        rows = list(csv.DictReader(open(trace + "/steps.csv")))
        selected = [r for r in rows if lo <= float(r["t"]) <= hi]
        print(f"=== step trace {trace}, every 10th accepted step: t ms : h us / history error")
        print(" ".join(f"{float(r['t']) * 1e3:.2f}:{float(r['h']) * 1e6:.0f}/{float(r['history_error']):.2f}" for r in selected[::10]))


if __name__ == "__main__":
    main(sys.argv[1:])
