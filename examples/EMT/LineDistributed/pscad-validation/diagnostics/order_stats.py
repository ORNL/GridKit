"""BDF order and step statistics of line_distributed runs (no NumPy; runs in the container).

    python3 order_stats.py <window_start_s> <window_end_s> <run dir> [<run dir> ...]
"""

import csv
import json
import statistics
import sys
from collections import Counter


def main(argv):
    lo, hi = float(argv[0]), float(argv[1])
    print(f"{'run':34s} {'tau us':>7s} {'median h us':>12s} {'h/tau':>6s} {'time>tau %':>10s} {'steps':>7s}  order share in window (order: % of steps)")
    for path in argv[2:]:
        data = json.load(open(path + "/input.json"))
        tau = min(d for d in data["line"]["H"]["delays"] if d > 0)
        rows = list(csv.DictReader(open(path + "/steps.csv")))
        window = [r for r in rows if lo <= float(r["t"]) <= hi]
        hs = [float(r["h"]) for r in window]
        above = sum(h for h in hs if h > tau) / sum(hs)
        orders = Counter(int(r["order"]) for r in window if "order" in r and r["order"])
        share = ", ".join(f"{k}: {100 * v / len(window):.0f}%" for k, v in sorted(orders.items()))
        print(f"{path.split('/')[-1]:34s} {tau * 1e6:7.1f} {statistics.median(hs) * 1e6:12.1f} {statistics.median(hs) / tau:6.2f} {100 * above:10.1f} {len(rows):7d}  {share}")


if __name__ == "__main__":
    main(sys.argv[1:])
