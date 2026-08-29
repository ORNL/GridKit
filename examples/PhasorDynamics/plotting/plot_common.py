"""Shared input handling and figure style for PhasorDynamics plots."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
MPL_CACHE = REPO_ROOT / "build/examples/PhasorDynamics/plotting/.matplotlib"
MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_MUTED = "#52514e"
GRID = "#c9c8c3"
BLUE = "#2a78d6"
BLUE_RAMP = [
    "#cde2fb",
    "#9ec5f4",
    "#6da7ec",
    "#3987e5",
    "#2a78d6",
    "#256abf",
    "#1c5cab",
    "#184f95",
    "#104281",
    "#0d366b",
]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def resolve_path(parent: Path, value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else parent / path


def load_study(study_path: Path) -> tuple[dict[str, Any], dict[str, Any], Path]:
    study_path = study_path.resolve()
    study = load_json(study_path)
    if "system_model_file" not in study:
        raise ValueError(f"{study_path} has no system_model_file")
    case_path = resolve_path(study_path.parent, study["system_model_file"])
    if not case_path.is_file():
        raise FileNotFoundError(f"case file not found: {case_path}")
    return study, load_json(case_path), study_path


def configured_output(study: dict[str, Any], study_path: Path, key: str) -> Path:
    value = study.get(key)
    if not value:
        raise ValueError(f"{study_path} has no {key}")
    return resolve_path(study_path.parent, value)


def event_times(study: dict[str, Any]) -> list[tuple[float, str]]:
    return [
        (float(event["time"]), str(event.get("type", "")))
        for event in study.get("events", [])
    ]


def case_name(case: dict[str, Any], study: dict[str, Any]) -> str:
    header = case.get("header", {})
    return str(header.get("case_name") or Path(study["system_model_file"]).stem)


def default_plot_path(study_path: Path, study: dict[str, Any], suffix: str) -> Path:
    stem = Path(study["system_model_file"]).stem
    return study_path.parent / f"{stem}_{suffix}.png"


def apply_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Liberation Serif", "Nimbus Roman", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "axes.linewidth": 0.6,
            "axes.edgecolor": INK_MUTED,
            "axes.labelcolor": INK,
            "axes.facecolor": SURFACE,
            "figure.facecolor": SURFACE,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "grid.color": GRID,
            "grid.linewidth": 0.4,
            "grid.alpha": 0.6,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "xtick.color": INK_MUTED,
            "ytick.color": INK_MUTED,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
        }
    )
