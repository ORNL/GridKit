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
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

# Keys that only make sense for GridKit's own regression tests.
# Keys that make the shipped solver a regression test: they grade the run
# against a stored reference trace. Any perturbed scenario "fails" that
# comparison by construction, so they are stripped for scenario use.
# `output_file` is deliberately NOT here: without it DynamicSimulation
# writes no monitor csv at all.
_TEST_ONLY_KEYS = (
    "reference_file",
    "error_tolerance",
    "error_type",
    "abs_err_threshold",
)

__all__ = [
    "list_faults",
    "fault_element_id",
    "set_fault_bus",
    "enable_monitors",
    "BUS_MONITORS",
    "GENROU_MONITORS",
    "apply_solver_overrides",
    "make_solver",
    "validate_solver",
    "describe_solver",
    "set_fault_impedance",
    "write_solver",
]

# Valid "mon" entries, from the GridKit enums: BusMonitorableVariables in
# Bus/BusData.hpp and GenrouMonitorableVariables in
# SynchronousMachine/GENROU/GenrouData.hpp. An unknown name is logged as an
# error by the parser and then ignored, so typos cost you the channel silently.
BUS_MONITORS = ("Vr", "Vi")
GENROU_MONITORS = ("ir", "ii", "p", "q", "delta", "omega", "speed")


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


def _resolve_bus(case: dict, bus: int | str) -> int:
    """Bus number from a number or a bus name."""
    if not isinstance(bus, str):
        return int(bus)
    buses = case.get("buses", [])
    hits = [b["number"] for b in buses if b.get("name") == bus]
    if not hits:
        hits = [
            b["number"] for b in buses if bus.lower() in str(b.get("name", "")).lower()
        ]
    if len(hits) != 1:
        raise KeyError(f"bus name {bus!r} matched {len(hits)} buses")
    return int(hits[0])


def set_fault_bus(case, bus: int | str, fault_id: str | None = None) -> dict:
    """Move a `BusFault` to another bus. Returns a modified copy.

    Relocating one device is preferred over shipping a fault per bus: every
    extra `BusFault` is another component GridKit assembles, and the pre-#577
    case carried 200 of them. `fault_id` is only needed when a case has more
    than one fault.
    """
    case = _load(case, "case")
    faults = [d for d in case.get("devices", []) if d.get("class") == "BusFault"]
    if not faults:
        raise KeyError("case has no BusFault device to move")
    if fault_id is not None:
        sel = [d for d in faults if d.get("id") == fault_id]
        if not sel:
            raise KeyError(
                f"no BusFault with id {fault_id!r}; have "
                f"{[d.get('id') for d in faults]}"
            )
        target = sel[0]
    elif len(faults) == 1:
        target = faults[0]
    else:
        raise KeyError(
            f"{len(faults)} BusFaults present, pass fault_id= to pick one: "
            f"{[d.get('id') for d in faults]}"
        )
    target["ports"]["bus"] = _resolve_bus(case, bus)
    return case


def enable_monitors(
    case,
    bus: Sequence[str] | None = ("Vr", "Vi"),
    genrou: Sequence[str] | None = ("delta", "omega"),
) -> dict:
    """Set which channels land in the monitor csv. Returns a modified copy.

    The post-#577 ACTIVSg200 case monitors `omega` only: every bus carries
    `"mon": []`, so no bus voltage reaches the csv and nothing PMU-observable
    can be computed. Pass None to leave a group untouched.

    Valid names are BUS_MONITORS and GENROU_MONITORS.
    """
    case = _load(case, "case")
    for group, names, valid in (
        ("bus", bus, BUS_MONITORS),
        ("genrou", genrou, GENROU_MONITORS),
    ):
        if names is None:
            continue
        bad = sorted(set(names) - set(valid))
        if bad:
            raise ValueError(f"invalid {group} monitor(s) {bad}; valid: {list(valid)}")
    if bus is not None:
        for b in case.get("buses", []):
            b["mon"] = list(bus)
    if genrou is not None:
        for d in case.get("devices", []):
            if d.get("class") == "Genrou":
                d["mon"] = list(genrou)
    return case


def apply_solver_overrides(
    solver: dict, overrides: dict | None, strip_test_keys: bool = True
) -> dict:
    """Patch a solver dict. Mirrors `gridkit_utils.make_run_dir` semantics.

    `events` may be a complete replacement list, or a list of patches each
    carrying an `index` key naming the event to modify.

    `strip_test_keys` removes the validation-only keys the shipped solver
    carries. Leave it on for scenario runs: with `reference_file` present
    DynamicSimulation grades the run against a stored trace and exits
    nonzero whenever the case has been perturbed.
    """
    solver = copy.deepcopy(solver)
    if strip_test_keys:
        for k in _TEST_ONLY_KEYS:
            solver.pop(k, None)
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
