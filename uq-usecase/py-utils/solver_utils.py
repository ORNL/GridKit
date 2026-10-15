"""Build and edit GridKit `solver.json` files.

The shipped Illinois solver faults bus 2 at t=10 s and runs to t=20 s, which
wastes most of the run on pre-fault steady state. Typical UQ settings are a
fault at 1 s cleared at 1.1 s with tmax=10 s.

`element_id` in a solver event is the **0-based ordinal of the BusFault device**
in `case.json`, not a bus number. GridKit assigns it at load time
(`gridkit_fault_id = gridkit_fault_indices_.size()` in `SystemModelImpl.hpp`).
For Illinois the ordering is not bus order -- `fault_1` is bus 2 and `fault_2`
is bus 1 -- so always resolve it with `fault_element_id()` rather than by
arithmetic on the bus number.

`apply_solver_overrides` deliberately mirrors the `events` semantics already used
by `gridkit_utils.make_run_dir` (full replacement, or selective patching by
`index`) so the epistemic and aleatoric tracks behave identically. That logic is
inline in `make_run_dir`; it is duplicated rather than refactored here to avoid
disturbing the working epistemic pipeline.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pandas as pd

# Keys that only make sense for GridKit's own regression tests.
_TEST_ONLY_KEYS = ("reference_file", "error_tolerance", "output_file")

__all__ = [
    "list_faults",
    "fault_element_id",
    "apply_solver_overrides",
    "make_solver",
    "validate_solver",
    "describe_solver",
    "set_fault_impedance",
    "write_solver",
]


def _load(obj, what: str) -> dict:
    if isinstance(obj, (str, Path)):
        with open(obj) as fh:
            return json.load(fh)
    if isinstance(obj, dict):
        return copy.deepcopy(obj)
    raise TypeError(f"{what} must be a path or dict, got {type(obj).__name__}")


def list_faults(case) -> pd.DataFrame:
    """Every `BusFault` in the case with the `element_id` a solver event needs."""
    case = _load(case, "case")
    names = {b["number"]: b.get("name", "") for b in case.get("buses", [])}
    rows = []
    for i, d in enumerate(
        d for d in case.get("devices", []) if d.get("class") == "BusFault"
    ):
        bus = d["ports"]["bus"]
        p = d.get("params", {})
        rows.append(
            {
                "element_id": i,
                "fault_id": d.get("id"),
                "bus": bus,
                "bus_name": names.get(bus, ""),
                "R": p.get("R"),
                "X": p.get("X"),
            }
        )
    return pd.DataFrame(rows)


def fault_element_id(case, bus: int | str) -> int:
    """Resolve a bus number or bus name to the solver `element_id`."""
    f = list_faults(case)
    if isinstance(bus, str):
        hit = f[f["bus_name"] == bus]
        if hit.empty:
            hit = f[f["bus_name"].str.contains(bus, case=False, na=False)]
    else:
        hit = f[f["bus"] == int(bus)]
    if hit.empty:
        raise KeyError(f"no BusFault for bus {bus!r}")
    if len(hit) > 1:
        raise KeyError(
            f"bus {bus!r} matches {len(hit)} faults: {hit['fault_id'].tolist()}"
        )
    return int(hit["element_id"].iloc[0])


def apply_solver_overrides(solver: dict, overrides: dict | None) -> dict:
    """Patch a solver dict. Mirrors `gridkit_utils.make_run_dir` semantics.

    `events` may be a complete replacement list, or a list of patches each
    carrying an `index` key naming the event to modify.
    """
    solver = copy.deepcopy(solver)
    if not overrides:
        return solver
    for key, val in overrides.items():
        if key == "events" and val and all("index" not in e for e in val):
            solver["events"] = copy.deepcopy(val)
        elif key == "events" and val:
            for patch in val:
                idx = patch["index"]
                for k, v in patch.items():
                    if k != "index":
                        solver["events"][idx][k] = v
        else:
            solver[key] = val
    return solver


def make_solver(
    base_solver,
    case=None,
    *,
    fault_bus: int | str | None = None,
    fault_time: float | None = None,
    clearing_time: float | None = None,
    fault_duration: float | None = None,
    tmax: float | None = None,
    dt_monitor: float | None = None,
    system_model_file: str | None = None,
    strip_test_keys: bool = True,
    overrides: dict | None = None,
) -> dict:
    """Return a solver dict with the fault and timing you actually want.

    Only what you pass is changed. `fault_bus` accepts a bus number or name and
    needs `case` to resolve the `element_id`. Give either `clearing_time` or
    `fault_duration`. `overrides` is applied last as a raw escape hatch.
    """
    solver = _load(base_solver, "base_solver")
    if strip_test_keys:
        for k in _TEST_ONLY_KEYS:
            solver.pop(k, None)

    if system_model_file is not None:
        solver["system_model_file"] = system_model_file
    if tmax is not None:
        solver["tmax"] = float(tmax)
    if dt_monitor is not None:
        solver["dt_monitor"] = float(dt_monitor)

    events = solver.get("events", [])
    on = next((e for e in events if e.get("type") == "fault_on"), None)
    off = next((e for e in events if e.get("type") == "fault_off"), None)

    if fault_bus is not None:
        if case is None:
            raise ValueError("fault_bus needs case= to resolve the element_id")
        eid = fault_element_id(case, fault_bus)
        for e in (on, off):
            if e is not None:
                e["element_id"] = eid

    if fault_time is not None and on is not None:
        on["time"] = float(fault_time)
    if clearing_time is not None and fault_duration is not None:
        raise ValueError("give clearing_time or fault_duration, not both")
    if off is not None:
        if clearing_time is not None:
            off["time"] = float(clearing_time)
        elif fault_duration is not None and on is not None:
            off["time"] = float(on["time"]) + float(fault_duration)

    if on is not None and off is not None and off["time"] <= on["time"]:
        raise ValueError(f"fault clears at {off['time']} but starts at {on['time']}")
    if "tmax" in solver and off is not None and solver["tmax"] <= off["time"]:
        raise ValueError(
            f"tmax={solver['tmax']} is not after fault clearing at {off['time']}"
        )

    # Validate last: raw overrides are applied after the structured arguments.
    return validate_solver(apply_solver_overrides(solver, overrides))


def validate_solver(solver: dict) -> dict:
    """Raise if the event timing is self-inconsistent. Returns the solver."""
    events = solver.get("events", [])
    on = next((e for e in events if e.get("type") == "fault_on"), None)
    off = next((e for e in events if e.get("type") == "fault_off"), None)
    if on is not None and off is not None:
        if off["time"] <= on["time"]:
            raise ValueError(
                f"fault_off at {off['time']}s is not after fault_on at {on['time']}s"
            )
        if on["element_id"] != off["element_id"]:
            raise ValueError(
                f"fault_on targets element_id {on['element_id']} but fault_off "
                f"targets {off['element_id']}; the fault would never clear"
            )
    last = max((e["time"] for e in events), default=None)
    if last is not None and "tmax" in solver and solver["tmax"] <= last:
        raise ValueError(
            f"tmax={solver['tmax']}s is not after the last event at {last}s"
        )
    return solver


def describe_solver(solver: dict, case=None) -> str:
    """One-line human summary, so a wrong fault location is visible."""
    bits = [f"tmax={solver.get('tmax')}s", f"dt_monitor={solver.get('dt_monitor')}"]
    faults = list_faults(case) if case is not None else None
    for e in solver.get("events", []):
        where = ""
        if faults is not None:
            hit = faults[faults["element_id"] == e.get("element_id")]
            if not hit.empty:
                where = f" bus {int(hit['bus'].iloc[0])} ({hit['bus_name'].iloc[0]})"
        bits.append(
            f"{e.get('type')}@{e.get('time')}s [element_id={e.get('element_id')}{where}]"
        )
    return "  ".join(bits)


def set_fault_impedance(
    case, bus: int | str, R: float | None = None, X: float | None = None
) -> dict:
    """Patch a `BusFault`'s R/X. These live in `case.json`, not `solver.json`."""
    case = _load(case, "case")
    eid = fault_element_id(case, bus)
    fault_devs = [d for d in case["devices"] if d.get("class") == "BusFault"]
    p = fault_devs[eid].setdefault("params", {})
    if R is not None:
        p["R"] = float(R)
    if X is not None:
        p["X"] = float(X)
    return case


def write_solver(solver: dict, path: str | Path, indent: int = 2) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(solver, indent=indent))
    return path
