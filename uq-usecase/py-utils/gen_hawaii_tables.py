"""
Generate hawaii_m_tables.md from Hawaii40_20231026.m using read_matpower_case.
Usage: python3 py-utils/gen_hawaii_tables.py
Output: cases/hawaii_m_tables.md
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from m_viz_utils import read_matpower_case

M_PATH = Path(
    "/Users/isatkaus/projects/scidac/scidac-data/Hawaii40/raw-tamu-data/Hawaii40_20231026.m"
)
OUT_PATH = Path(__file__).parent.parent / "cases" / "hawaii_m_tables.md"

cd = read_matpower_case(M_PATH)
bus = cd.bus.copy()
gen = cd.gen.copy()
branch = cd.branch.copy()
genfuel = cd.genfuel
gencost = cd.gencost

# attach fuel to gen table
if genfuel is not None:
    gen["fuel"] = list(genfuel)

# ── summary ─────────────────────────────────────────────────────────────────
pq = int((bus.BUS_TYPE == 1).sum())
pv = int((bus.BUS_TYPE == 2).sum())
sl = int((bus.BUS_TYPE == 3).sum())
load_buses = int((bus.PD > 0).sum())
total_pd = bus.PD.sum()
total_qd = bus.QD.sum()
v_min = bus.VM.min()
v_max = bus.VM.max()

n_branches = len(branch)
tap_nonunity = int(((branch.TAP != 0) & (branch.TAP != 1)).sum())
status0_branches = (
    int((branch.BR_STATUS == 0).sum()) if "BR_STATUS" in branch.columns else 0
)

n_gen = len(gen)
online_gen = gen[gen.GEN_STATUS == 1]
offline_gen = gen[gen.GEN_STATUS == 0]
n_online = len(online_gen)
n_offline = len(offline_gen)
total_pg = online_gen.PG.sum()
total_qg = online_gen.QG.sum()
total_pmax = online_gen.PMAX.sum()

if "fuel" in gen.columns:
    fuel_counts = online_gen.groupby("fuel").size().to_dict()
    fuel_pg = online_gen.groupby("fuel")["PG"].sum().to_dict()
    fuel_str = ", ".join(f"{v} {k}" for k, v in sorted(fuel_counts.items()))
else:
    fuel_str = "n/a"
    fuel_counts = {}
    fuel_pg = {}

# ── bus_name lookup ──────────────────────────────────────────────────────────
if "bus_name" in bus.columns:
    bus_name_col = "bus_name"
else:
    bus_name_col = None


def bname(row):
    if bus_name_col and bus_name_col in row.index:
        return row[bus_name_col]
    return ""


# ── build markdown ───────────────────────────────────────────────────────────
lines = []

lines.append(f"# MATPOWER base case tables ({M_PATH.name})\n")
lines.append(
    "Parsed via `read_matpower_case` from `matpowercaseframes`. Values are as stored in the\n"
    "`.m` file: powers in MW/MVAr, angles in degrees, voltages in pu or kV.\n"
)

# tap info
if tap_nonunity == 0:
    tap_note = (
        "All branches have `TAP=0` or `TAP=1` and `SHIFT=0` "
        "(no off-nominal tap transformers, no phase shifters)."
    )
else:
    tap_note = f"{tap_nonunity} branches have off-nominal tap ratios."
lines.append(tap_note + "\n")

# branch status
br_status_note = (
    f"All {n_branches} branches have `BR_STATUS=1` (in-service)."
    if status0_branches == 0
    else f"{status0_branches} branches have `BR_STATUS=0` (out-of-service)."
)
lines.append(br_status_note + "\n")
lines.append(
    "All summary numbers below are computed by parsing the tables in this file.\n"
)

# summary table
lines.append("## Summary statistics\n")
lines.append("| Quantity | Value |")
lines.append("| :--- | :--- |")
lines.append(
    f"| **Buses** | {len(bus)} total: {pq} PQ (type 1), {pv} PV (type 2), {sl} slack (type 3) |"
)
lines.append(f"| **Load buses** (PD > 0) | {load_buses} buses |")
lines.append(f"| **Total PD** | {total_pd:.2f} MW |")
lines.append(f"| **Total QD** | {total_qd:.2f} MVAr |")
lines.append(
    f"| **Base-case V range** | [{v_min:.3f}, {v_max:.3f}] pu (VM column from `mpc.bus` as stored in the .m file) |"
)
lines.append(f"| **Branches** | {n_branches} (all in-service) |")
lines.append(
    f"| **Generators total** | {n_gen}: {n_online} online, {n_offline} offline |"
)
offline_fuel = (
    ", ".join(sorted(offline_gen["fuel"].unique()))
    if "fuel" in gen.columns and n_offline > 0
    else "n/a"
)
lines.append(f"| **Offline gens** | {n_offline} ({offline_fuel}) |")
lines.append(f"| **Online fuel mix** | {fuel_str} ({n_online} total) |")
lines.append(f"| **Total PG (online)** | {total_pg:.2f} MW |")
lines.append(f"| **Total QG (online)** | {total_qg:.2f} MVAr |")
lines.append(f"| **Total PMAX (online)** | {total_pmax:.2f} MW |")
lines.append(
    f"| **Power imbalance PG - PD** | {total_pg - total_pd:.2f} MW (losses + slack dispatch) |"
)
lines.append("")

# ── mpc.bus table ────────────────────────────────────────────────────────────
lines.append("## mpc.bus (37 buses)\n")
lines.append(
    "`BUS_TYPE`: 1=PQ, 2=PV, 3=slack. `BASE_KV`: 69 or 138 kV.\n"
    "`VM`/`VA`: pre-computed PF solution stored in the .m file (PowerWorld export).\n"
)

bus_cols = [
    "BUS_I",
    "BUS_TYPE",
    "PD",
    "QD",
    "GS",
    "BS",
    "BUS_AREA",
    "VM",
    "VA",
    "BASE_KV",
    "ZONE",
    "VMAX",
    "VMIN",
]
if bus_name_col:
    bus_cols.append(bus_name_col)

header = "| " + " | ".join(bus_cols) + " |"
sep = "| " + " | ".join(":---" for _ in bus_cols) + " |"
lines.append(header)
lines.append(sep)
for _, row in bus.iterrows():
    vals = []
    for c in bus_cols:
        v = row[c]
        if isinstance(v, float):
            vals.append(f"{v:.3g}" if c not in ("VM", "VA") else f"{v:.3f}")
        else:
            vals.append(str(v))
    lines.append("| " + " | ".join(vals) + " |")
lines.append("")

# ── mpc.gen table ────────────────────────────────────────────────────────────
lines.append("## mpc.gen (45 generators)\n")
lines.append(
    "`gen_row` is the 1-based row index in `mpc.gen`. "
    "Rows marked `*` have `GEN_STATUS = 0` (offline).\n"
)

lines.append(
    "| gen_row | bus | bus_name | PG (MW) | QG (MVAr) | PMAX (MW) | PMIN (MW) | QMAX (MVAr) | QMIN (MVAr) | MBASE (MVA) | status | fuel |"
)
lines.append(
    "|--------:|----:|:---------|--------:|----------:|----------:|----------:|------------:|------------:|------------:|:------:|:-----|"
)

# build bus_name lookup dict
if bus_name_col:
    bname_dict = dict(zip(bus.BUS_I, bus[bus_name_col]))
else:
    bname_dict = {}

for i, (_, row) in enumerate(gen.iterrows(), 1):
    bus_i = int(row.GEN_BUS)
    bname_str = bname_dict.get(bus_i, "")
    status_str = "✓" if row.GEN_STATUS == 1 else "✗ off"
    fuel_str2 = row["fuel"] if "fuel" in row.index else ""
    pg = f"{row.PG:.2f}"
    qg = f"{row.QG:.2f}"
    pmax = f"{row.PMAX:.2f}"
    pmin = f"{row.PMIN:.2f}"
    qmax = f"{row.QMAX:.2f}"
    qmin = f"{row.QMIN:.2f}"
    mbase = f"{row.MBASE:.2f}" if "MBASE" in row.index else "100.00"
    lines.append(
        f"| {i:>3} | {bus_i} | {bname_str} | {pg} | {qg} | {pmax} | {pmin} | {qmax} | {qmin} | {mbase} | {status_str} | {fuel_str2} |"
    )
lines.append("")

# ── mpc.branch table ─────────────────────────────────────────────────────────
lines.append("## mpc.branch (57 branches)\n")
lines.append(
    "All branches in-service (`BR_STATUS=1`). "
    "Parallel circuits are common in Hawaii40.\n"
    "`TAP=0` means plain line (unity tap by MATPOWER convention); "
    "`TAP=1.0` means explicit unity transformer.\n"
)

br_cols = [
    "F_BUS",
    "T_BUS",
    "BR_R",
    "BR_X",
    "BR_B",
    "RATE_A",
    "TAP",
    "SHIFT",
    "BR_STATUS",
]
header = "| br_row | " + " | ".join(br_cols) + " |"
sep = "| ------:|" + "".join(" ---: |" for _ in br_cols)
lines.append(header)
lines.append(sep)
for i, (_, row) in enumerate(branch.iterrows(), 1):
    vals = []
    for c in br_cols:
        v = row[c]
        if isinstance(v, float):
            if c in ("BR_R", "BR_X", "BR_B"):
                vals.append(f"{v:.6f}")
            elif c in ("RATE_A",):
                vals.append(f"{v:.1f}")
            elif c in ("TAP", "SHIFT"):
                vals.append(f"{v:.5f}")
            else:
                vals.append(f"{v:.2f}")
        else:
            vals.append(str(v))
    lines.append("| " + str(i) + " | " + " | ".join(vals) + " |")
lines.append("")

# ── write ────────────────────────────────────────────────────────────────────
OUT_PATH.write_text("\n".join(lines) + "\n")
print(f"Written: {OUT_PATH}")
print(f"  buses={len(bus)}, gens={len(gen)}, branches={len(branch)}")
