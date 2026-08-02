#!/usr/bin/env python3
"""Render the tower cross-section of an overhead-line 1.0 document."""

import json
import sys
from collections import OrderedDict
from pathlib import Path

import yaml

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

PHASE_COLOR = {
    "a": "#C1443C",
    "b": "#D99A2B",
    "c": "#2E6F9E",
    "n": "#7A8894",
    "g": "#3F4A56",
}
PHASE_NAME = {
    "a": "phase a",
    "b": "phase b",
    "c": "phase c",
    "n": "neutral",
    "g": "grounded (shield)",
}

INK = "#1B2027"
MUTED = "#8A93A0"
FAINT = "#DCE1E8"
PAPER = "#FFFFFF"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Latin Modern Roman", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.edgecolor": FAINT,
    "text.color": INK,
    "axes.labelcolor": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "figure.facecolor": PAPER,
    "axes.facecolor": PAPER,
    "savefig.facecolor": PAPER,
})


def attachments(doc, path):
    """Attachment coordinates of the document's tower type."""
    towers = dict((doc.get("catalog") or {}).get("towers") or {})
    for inc in doc.get("include", []):
        catalog = yaml.safe_load((path.parent / inc).read_text(encoding="utf-8"))
        towers.update(catalog.get("towers") or {})
    return towers[doc["tower"]]["attachments"]


def groups(doc, points):
    """Group conductors by (phase, circuit), preserving file order."""
    out = OrderedDict()
    for cond in doc["conductors"]:
        point = points[cond["at"]]
        key = (cond["phase"], cond.get("circuit", 1))
        out.setdefault(key, []).append((point["x"], point["height"]))
    return out


def cluster(values, tol):
    """Collapse near-equal values into representative levels."""
    levels = []
    for v in sorted(values):
        if levels and abs(v - levels[-1][-1]) <= tol:
            levels[-1].append(v)
        else:
            levels.append([v])
    return [sum(g) / len(g) for g in levels]


def render(path, out_path):
    path = Path(path)
    doc = json.loads(path.read_text(encoding="utf-8"))
    unit = "m"
    grp = groups(doc, attachments(doc, path))
    multi_circuit = max(k[1] for k in grp) > 1

    xs = [x for pts in grp.values() for x, _ in pts]
    ys = [y for pts in grp.values() for _, y in pts]
    half = max(abs(min(xs)), abs(max(xs)))
    half = half + max(0.22 * half, 0.05 * max(ys))
    ytop = max(ys) * 1.10
    ybot = -0.05 * max(ys)
    span_x, span_y = 2 * half, ytop - ybot

    bundles = [pts for pts in grp.values() if len(pts) > 1]
    scale = min(8.6 / span_y, 6.4 / span_x)

    if bundles:
        fig, (ax, bx) = plt.subplots(
            1, 2,
            figsize=(span_x * scale + 3.5, span_y * scale + 1.5),
            gridspec_kw={"width_ratios": [1, 0.42]},
        )
    else:
        fig, ax = plt.subplots(figsize=(span_x * scale + 1.9, span_y * scale + 1.5))
        bx = None

    # ---- ground ---------------------------------------------------------
    ax.add_patch(Rectangle((-half, ybot), 2 * half, -ybot,
                           facecolor="#F1F3F6", edgecolor="none", zorder=0))
    ax.axhline(0, color="#AAB2BD", lw=1.1, zorder=1)

    # ---- height leaders -------------------------------------------------
    for lvl in cluster(ys, tol=0.02 * span_y):
        at_level = [x for x, y in zip(xs, ys) if abs(y - lvl) <= 0.02 * span_y]
        ax.plot([-half, min(at_level) - 0.03 * span_x], [lvl, lvl],
                color=FAINT, lw=0.8, ls=(0, (1, 3)), zorder=1)
        ax.text(-half + 0.015 * span_x, lvl + 0.008 * span_y,
                f"{lvl:g}", fontsize=7.5, color=MUTED, va="bottom", ha="left")

    # ---- centreline -----------------------------------------------------
    ax.plot([0, 0], [0, ytop * 0.985], color=FAINT, lw=0.9,
            ls=(0, (6, 4)), zorder=1)
    ax.text(0, ytop * 0.995, "centreline", fontsize=7, color=MUTED,
            ha="center", va="bottom")

    # ---- conductors -----------------------------------------------------
    for (phase, circuit), pts in grp.items():
        col = PHASE_COLOR[phase]
        gx = [p[0] for p in pts]
        gy = [p[1] for p in pts]
        ax.scatter(gx, gy, s=46, color=col, edgecolor=PAPER,
                   linewidth=1.1, zorder=4)
        label = f"{phase}$_{circuit}$" if (multi_circuit and phase in "abc") else phase
        ax.text(sum(gx) / len(gx), max(gy) + 0.022 * span_y, label,
                fontsize=10.5, color=col, ha="center", va="bottom", zorder=5)

    # ---- frame ----------------------------------------------------------
    ax.set_xlim(-half, half)
    ax.set_ylim(ybot, ytop * 1.045)
    ax.set_aspect("equal")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_position(("outward", 8))
    ax.spines["bottom"].set_color(FAINT)
    ax.set_yticks([])
    ax.tick_params(axis="x", length=3, labelsize=8)
    ax.set_xlabel(f"$x$  [{unit}]", fontsize=9, labelpad=4)
    ax.text(-half, ytop * 1.115, doc.get("name", "Tower geometry"),
            fontsize=12.5, ha="left", va="bottom")
    n_cond = len(doc["conductors"])
    ax.text(-half, ytop * 1.075,
            f"{n_cond} conductors  ·  heights in {unit} above ground at base",
            fontsize=8, color=MUTED, ha="left", va="bottom")

    seen = []
    for (phase, _) in grp:
        if phase not in seen:
            seen.append(phase)
    handles = [plt.Line2D([], [], marker="o", ls="", markersize=5.5,
                          markerfacecolor=PHASE_COLOR[p], markeredgecolor=PAPER,
                          label=PHASE_NAME[p]) for p in seen]
    leg = ax.legend(handles=handles, loc="lower right", frameon=False,
                    fontsize=8, handletextpad=0.4, labelspacing=0.35,
                    borderaxespad=0.2)
    for t in leg.get_texts():
        t.set_color(MUTED)

    # ---- bundle detail panel -------------------------------------------
    if bx is not None:
        pts = max(bundles, key=len)
        cx = sum(p[0] for p in pts) / len(pts)
        cy = sum(p[1] for p in pts) / len(pts)
        rel = [(p[0] - cx, p[1] - cy) for p in pts]
        pair = min(
            ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** 0.5
            for i, a in enumerate(rel) for b in rel[i + 1:]
        )
        r = max(max(abs(dx), abs(dy)) for dx, dy in rel) * 1.9

        import math
        hull = sorted(rel, key=lambda p: math.atan2(p[1], p[0]))
        bx.plot([p[0] for p in hull] + [hull[0][0]],
                [p[1] for p in hull] + [hull[0][1]],
                color=FAINT, lw=0.9, ls=(0, (3, 3)), zorder=1)
        bx.scatter([p[0] for p in rel], [p[1] for p in rel], s=46,
                   color="#4A5560", edgecolor=PAPER, linewidth=1.1, zorder=3)
        bx.scatter([0], [0], s=10, marker="+", color=MUTED, zorder=2)
        bx.set_xlim(-r, r)
        bx.set_ylim(-r, r)
        bx.set_aspect("equal")
        bx.set_xticks([])
        bx.set_yticks([])
        for side in bx.spines.values():
            side.set_visible(False)
        bx.set_anchor("N")
        bx.set_title("bundle detail", fontsize=9, color=INK, pad=8)
        bx.text(0, -r * 0.88,
                f"{len(pts)} sub-conductors  ·  {pair:g} {unit} apart",
                fontsize=7.5, color=MUTED, ha="center", va="center")

    fig.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight", pad_inches=0.32)
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    for p in sys.argv[1:]:
        src = Path(p)
        out_dir = src.parent / "figures"
        out_dir.mkdir(exist_ok=True)
        stem = src.name.rsplit(".line.json", 1)[0]
        print(render(p, out_dir / f"{stem}.png"))
