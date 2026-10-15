"""Scenario perturbations on MATPOWER frames.

Operates on the frame dict produced by `aux_to_m.aux_to_frames` (bootstrap) or
`aux_to_m.m_to_frames` (any `.m`, including Sienna / PCM output), so the same
functions apply regardless of where the case came from.

Every function returns a new frame dict; inputs are not mutated.

Typical scenario loop::

    frames = m_to_frames(base_m)
    frames = scale_load(frames, 0.92)
    frames = set_gen_output(frames, {"65_1": 31.4}, fuel="wind")
    write_matpower_m(frames, "scenario_0042.m")
    # -> PM.jl solve -> build_case_from_m -> DynamicSimulation
"""

from __future__ import annotations

import copy
import warnings

import numpy as np
import pandas as pd

__all__ = [
    "copy_frames",
    "gen_index",
    "scale_load",
    "set_bus_load",
    "set_gen_status",
    "set_gen_output",
    "summarize",
]


def copy_frames(frames: dict) -> dict:
    return {
        k: (v.copy() if isinstance(v, (pd.DataFrame, pd.Series)) else copy.copy(v))
        for k, v in frames.items()
    }


def gen_index(frames: dict) -> pd.Series:
    """`"{bus}_{rank}"` per `mpc.gen` row, matching the GridKit Genrou id."""
    bus = frames["gen"]["GEN_BUS"].astype(int)
    multi = bus.value_counts()
    if (multi > 1).any():
        warnings.warn(
            f"buses with >1 generator: {multi[multi > 1].to_dict()}; ids are "
            "positional, so .m row order must match the base case.",
            stacklevel=2,
        )
    return bus.astype(str) + "_" + (bus.groupby(bus).cumcount() + 1).astype(str)


# --------------------------------------------------------------------------
# load
# --------------------------------------------------------------------------


def scale_load(frames: dict, factor, buses=None) -> dict:
    """Scale `PD`/`QD`. `factor` is a scalar or a `{bus: factor}` mapping."""
    f = copy_frames(frames)
    bus = f["bus"]
    bus_i = bus["BUS_I"].astype(int)
    if isinstance(factor, dict):
        mult = bus_i.map(factor).fillna(1.0).to_numpy(dtype=float)
    else:
        mult = np.full(len(bus), float(factor))
        if buses is not None:
            mult = np.where(bus_i.isin(list(buses)), mult, 1.0)
    bus["PD"] = bus["PD"].to_numpy() * mult
    bus["QD"] = bus["QD"].to_numpy() * mult
    return f


def set_bus_load(
    frames: dict,
    pd_mw: dict,
    qd_mvar: dict | None = None,
    keep_power_factor: bool = True,
) -> dict:
    """Set `PD` at specific buses.

    With `keep_power_factor` (default) `QD` is scaled by the same ratio as `PD`,
    which is what an hourly load timeseries normally implies. Pass `qd_mvar` to
    set reactive power explicitly instead.
    """
    f = copy_frames(frames)
    bus = f["bus"]
    bus_i = bus["BUS_I"].astype(int)
    for b, p_new in pd_mw.items():
        m = bus_i == int(b)
        if not m.any():
            raise KeyError(f"bus {b} not in case")
        p_old = float(bus.loc[m, "PD"].iloc[0])
        bus.loc[m, "PD"] = float(p_new)
        if qd_mvar is not None and b in qd_mvar:
            bus.loc[m, "QD"] = float(qd_mvar[b])
        elif keep_power_factor and p_old != 0.0:
            bus.loc[m, "QD"] = float(bus.loc[m, "QD"].iloc[0]) * float(p_new) / p_old
    return f


# --------------------------------------------------------------------------
# generators
# --------------------------------------------------------------------------


def _headroom_redistribute(
    gen: pd.DataFrame, deficit_mw: float, exclude: np.ndarray, slack_bus: set[int]
) -> dict:
    """Share `deficit_mw` across online gens in proportion to PMAX - PG."""
    on = (gen["GEN_STATUS"] == 1).to_numpy() & ~exclude
    if not on.any():
        raise ValueError("no online generators left to absorb the deficit")
    # The slack already balances losses; leave it out of the explicit share.
    pick = on & ~gen["GEN_BUS"].astype(int).isin(slack_bus).to_numpy()
    if not pick.any():
        pick = on
    headroom = (gen["PMAX"].to_numpy() - gen["PG"].to_numpy())[pick].clip(min=0.0)
    total = headroom.sum()
    if total <= 0:
        raise ValueError(
            f"no headroom to absorb {deficit_mw:.1f} MW; online units are at PMAX"
        )
    gen.loc[pick, "PG"] = gen.loc[pick, "PG"].to_numpy() + headroom * (
        deficit_mw / total
    )
    return {
        "absorbed_mw": float(deficit_mw),
        "n_units": int(pick.sum()),
        "headroom_mw": float(total),
        "capped": bool(deficit_mw > total),
    }


def set_gen_status(
    frames: dict, status: dict, redistribute: str | None = "headroom"
) -> dict:
    """Switch generators on/off. `status` maps gen id or bus number to 0/1.

    Switching a unit off leaves its output unserved. With
    `redistribute="headroom"` the lost MW is shared across the remaining online
    units in proportion to their spare capacity, which keeps the slack from
    absorbing everything. Use `redistribute=None` to let the PF slack take it.
    """
    f = copy_frames(frames)
    gen = f["gen"]
    gid = gen_index(f)
    bus_i = gen["GEN_BUS"].astype(int)
    slack = set(f["bus"].loc[f["bus"]["BUS_TYPE"] == 3, "BUS_I"].astype(int))

    touched = np.zeros(len(gen), dtype=bool)
    lost = 0.0
    for key, val in status.items():
        m = (
            (gid == str(key)).to_numpy()
            if isinstance(key, str)
            else (bus_i == int(key)).to_numpy()
        )
        if not m.any():
            raise KeyError(f"generator {key!r} not found")
        touched |= m
        if int(val) == 0:
            lost += float(gen.loc[m, "PG"].sum())
            gen.loc[m, ["PG", "QG"]] = 0.0
            gen.loc[m, "GEN_STATUS"] = 0
        else:
            gen.loc[m, "GEN_STATUS"] = 1

    report = {
        "switched": {str(k): int(v) for k, v in status.items()},
        "unserved_mw": lost,
    }
    if lost > 0 and redistribute == "headroom":
        report["redistribution"] = _headroom_redistribute(gen, lost, touched, slack)
    elif lost > 0:
        report["redistribution"] = "left to PF slack"
    f["_perturb_report"] = report
    return f


def set_gen_output(
    frames: dict,
    pg_mw: dict | None = None,
    fuel: str | None = None,
    scale: float | None = None,
) -> dict:
    """Set or scale `PG` on selected generators (wind, PV, any non-dispatchable).

    Either pass `pg_mw` as `{gen_id or bus: MW}`, or `scale` together with
    `fuel` to scale every unit of that fuel type. Units whose `PMIN == PMAX`
    (common for non-dispatchables) have both limits moved with `PG`, otherwise
    the result is infeasible for a later OPF.
    """
    f = copy_frames(frames)
    gen = f["gen"]
    gid = gen_index(f)
    bus_i = gen["GEN_BUS"].astype(int)

    if scale is not None:
        if fuel is None or "genfuel" not in f:
            raise ValueError("scale= requires fuel= and a genfuel table")
        sel = (
            pd.Series(list(f["genfuel"]), index=gen.index)
            .astype(str)
            .str.lower()
            .eq(fuel.lower())
        )
        pg_mw = {g: float(p) * scale for g, p, s in zip(gid, gen["PG"], sel) if s}

    for key, p_new in (pg_mw or {}).items():
        m = (
            (gid == str(key)).to_numpy()
            if isinstance(key, str)
            else (bus_i == int(key)).to_numpy()
        )
        if not m.any():
            raise KeyError(f"generator {key!r} not found")
        fixed = np.isclose(gen.loc[m, "PMIN"], gen.loc[m, "PMAX"])
        gen.loc[m, "PG"] = float(p_new)
        if fixed.all():
            gen.loc[m, ["PMIN", "PMAX"]] = float(p_new)
        else:
            gen.loc[m, "PMAX"] = np.maximum(gen.loc[m, "PMAX"], float(p_new))
            gen.loc[m, "PMIN"] = np.minimum(gen.loc[m, "PMIN"], float(p_new))
    return f


def summarize(frames: dict, label: str = "") -> pd.Series:
    """Quantities worth eyeballing before writing a scenario `.m`."""
    bus, gen = frames["bus"], frames["gen"]
    on = gen["GEN_STATUS"] == 1
    slack = set(bus.loc[bus["BUS_TYPE"] == 3, "BUS_I"].astype(int))
    slack_pg = float(gen.loc[on & gen["GEN_BUS"].astype(int).isin(slack), "PG"].sum())
    slack_max = float(gen.loc[gen["GEN_BUS"].astype(int).isin(slack), "PMAX"].sum())
    return pd.Series(
        {
            "label": label,
            "total_PD_MW": float(bus["PD"].sum()),
            "total_QD_MVAr": float(bus["QD"].sum()),
            "total_PG_MW": float(gen.loc[on, "PG"].sum()),
            "gens_online": int(on.sum()),
            "PG_minus_PD_MW": float(gen.loc[on, "PG"].sum() - bus["PD"].sum()),
            "slack_PG_MW": slack_pg,
            "slack_PMAX_MW": slack_max,
            "slack_headroom_MW": slack_max - slack_pg,
            "n_PG_over_PMAX": int(
                (gen.loc[on, "PG"] > gen.loc[on, "PMAX"] + 1e-6).sum()
            ),
        }
    )
