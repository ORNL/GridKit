"""Utilities for loading MATPOWER .m cases in m_viz.ipynb."""

from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
from pathlib import Path
import re
import shlex
from typing import Any

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from matpowercaseframes import CaseFrames


@dataclass
class MatpowerCaseData:
    """Normalized MATPOWER case tables and lightweight metadata."""

    case_name: str
    source_path: str
    bus: pd.DataFrame
    gen: pd.DataFrame
    branch: pd.DataFrame
    genfuel: pd.Index | None = None
    gencost: pd.DataFrame | None = None
    bus_name: pd.Index | None = None  # mpc.bus_name cell array, one entry per bus row
    base_mva: float | None = None
    version: str | None = None


@dataclass
class GeoPrepResult:
    """Result bundle for geo enrichment and split-point preprocessing."""

    case_data: MatpowerCaseData
    bus_geo_df: pd.DataFrame
    split_applied: bool
    n_unique_bus_locations_before: int
    n_unique_bus_locations_after: int


def read_matpower_case(m_file_path: str | Path) -> MatpowerCaseData:
    """Load a MATPOWER .m file and return normalized case tables.

    Uses the standard ``matpowercaseframes`` package for the numeric matrices
    (bus, gen, branch, gencost) and a built-in curly-brace parser for the
    string cell-arrays that the standard package does not handle
    (``mpc.genfuel``, ``mpc.gentype``, ``mpc.bus_name``).
    """

    source_path = Path(m_file_path).expanduser().resolve()
    if not source_path.exists():
        raise FileNotFoundError(f"MATPOWER case file not found: {source_path}")

    case = CaseFrames(str(source_path))

    bus_df = normalize_bus_df(case.bus)
    gen_df = normalize_gen_df(case.gen)
    branch_df = normalize_branch_df(case.branch)

    # Parse curly-brace cell arrays directly from the .m file text.
    m_text = source_path.read_text(encoding="utf-8", errors="replace")
    genfuel = _parse_m_cell_array(m_text, "genfuel")
    bus_name_list = _parse_m_cell_array(m_text, "bus_name")

    genfuel_index = pd.Index(genfuel, name="genfuel") if genfuel else None
    bus_name_index = pd.Index(bus_name_list, name="bus_name") if bus_name_list else None

    # Attach bus_name as a column on bus_df so name-based lookups work.
    if bus_name_index is not None and len(bus_name_index) == len(bus_df):
        bus_df = bus_df.copy()
        bus_df["bus_name"] = bus_name_index.values

    return MatpowerCaseData(
        case_name=getattr(case, "name", source_path.stem),
        source_path=str(source_path),
        bus=bus_df,
        gen=gen_df,
        branch=branch_df,
        genfuel=genfuel_index,
        gencost=_copy_df_attr(case, "gencost"),
        bus_name=bus_name_index,
        base_mva=_coerce_optional_float(getattr(case, "baseMVA", None)),
        version=_coerce_optional_str(getattr(case, "version", None)),
    )


def normalize_bus_df(bus_df: pd.DataFrame) -> pd.DataFrame:
    """Return bus table indexed by BUS_I with stable numeric typing."""

    df = bus_df.copy()
    if "BUS_I" not in df.columns:
        df = df.reset_index()
    if "BUS_I" not in df.columns:
        raise ValueError("MATPOWER bus table is missing required column 'BUS_I'")

    df["BUS_I"] = _to_int_series(df["BUS_I"], "BUS_I")
    df = df.set_index("BUS_I", drop=False).sort_index()
    return df


def normalize_gen_df(gen_df: pd.DataFrame) -> pd.DataFrame:
    """Return generator table with numeric GEN_BUS and a stable row index."""

    df = gen_df.copy().reset_index(drop=False)
    if "GEN_BUS" not in df.columns:
        raise ValueError(
            "MATPOWER generator table is missing required column 'GEN_BUS'"
        )

    df["GEN_BUS"] = _to_int_series(df["GEN_BUS"], "GEN_BUS")
    df = df.rename(columns={df.columns[0]: "gen_row"})
    # Name the GEN_STATUS column (col 8 in MATPOWER gen, index 7 after gen_row prepend)
    if "GEN_STATUS" not in df.columns and len(df.columns) > 8:
        df = df.rename(columns={df.columns[8]: "GEN_STATUS"})
    df = df.set_index("gen_row", drop=False)
    return df


def normalize_branch_df(branch_df: pd.DataFrame) -> pd.DataFrame:
    """Return branch table with numeric endpoint columns and a stable row index."""

    df = branch_df.copy().reset_index(drop=False)
    required_cols = ["F_BUS", "T_BUS"]
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(
            f"MATPOWER branch table is missing required columns: {', '.join(missing)}"
        )

    df["F_BUS"] = _to_int_series(df["F_BUS"], "F_BUS")
    df["T_BUS"] = _to_int_series(df["T_BUS"], "T_BUS")
    df = df.rename(columns={df.columns[0]: "branch_row"})
    df = df.set_index("branch_row", drop=False)
    return df


def summarize_case(case_data: MatpowerCaseData) -> pd.Series:
    """Return a concise summary for quick notebook display."""

    return pd.Series(
        {
            "case_name": case_data.case_name,
            "source_path": case_data.source_path,
            "base_mva": case_data.base_mva,
            "version": case_data.version,
            "n_buses": len(case_data.bus),
            "n_generators": len(case_data.gen),
            "n_branches": len(case_data.branch),
            "n_unique_gen_buses": int(case_data.gen["GEN_BUS"].nunique()),
            "n_unique_branch_from_buses": int(case_data.branch["F_BUS"].nunique()),
            "n_unique_branch_to_buses": int(case_data.branch["T_BUS"].nunique()),
        }
    )


def validate_case_tables(case_data: MatpowerCaseData) -> dict[str, Any]:
    """Validate key MATPOWER table relationships for notebook sanity checks."""

    bus_ids = set(case_data.bus["BUS_I"].tolist())
    gen_bus_ids = set(case_data.gen["GEN_BUS"].tolist())
    from_bus_ids = set(case_data.branch["F_BUS"].tolist())
    to_bus_ids = set(case_data.branch["T_BUS"].tolist())

    missing_gen_buses = sorted(gen_bus_ids - bus_ids)
    missing_from_buses = sorted(from_bus_ids - bus_ids)
    missing_to_buses = sorted(to_bus_ids - bus_ids)

    required_checks = {
        "bus_has_BUS_I": "BUS_I" in case_data.bus.columns,
        "gen_has_GEN_BUS": "GEN_BUS" in case_data.gen.columns,
        "branch_has_F_BUS": "F_BUS" in case_data.branch.columns,
        "branch_has_T_BUS": "T_BUS" in case_data.branch.columns,
    }

    return {
        "required_columns": required_checks,
        "all_required_columns_present": all(required_checks.values()),
        "missing_gen_bus_ids": missing_gen_buses,
        "missing_branch_from_bus_ids": missing_from_buses,
        "missing_branch_to_bus_ids": missing_to_buses,
        "gen_buses_are_valid": len(missing_gen_buses) == 0,
        "branch_from_buses_are_valid": len(missing_from_buses) == 0,
        "branch_to_buses_are_valid": len(missing_to_buses) == 0,
    }


def validation_report_df(case_data: MatpowerCaseData) -> pd.DataFrame:
    """Render validation output as a two-column DataFrame."""

    report = validate_case_tables(case_data)
    rows = []
    for key, value in report.items():
        rows.append({"check": key, "value": value})
    return pd.DataFrame(rows)


def create_bus_geo_df_from_gic(file_path: str | Path) -> pd.DataFrame:
    """Parse TAMU .gic and return bus-indexed geographic table with lat/lon."""

    gic_path = Path(file_path).expanduser().resolve()
    if not gic_path.exists():
        raise FileNotFoundError(f"GIC file not found: {gic_path}")

    with open(gic_path, "r") as f:
        lines = f.readlines()

    data_lines = []
    for line in lines[1:]:
        if line.strip().startswith("0"):
            break
        data_lines.append(line)

    name_str = "\n".join([line.split("'")[1] for line in data_lines])
    data_str = "".join(
        [line.split("'")[0] + line.split("'")[-1] for line in data_lines]
    )

    substation_df = pd.DataFrame(
        {"index": range(1, len(data_lines) + 1), "name": name_str.split("\n")}
    )
    substation_df.set_index("index", inplace=True)

    df_data = pd.read_csv(
        StringIO(data_str),
        sep=r"\s+",
        header=None,
        names=["unknown", "lat", "lon", "value"],
    )
    substation_df = pd.concat([substation_df, df_data], axis=1)

    start_idx = None
    for i, line in enumerate(lines):
        if "Begin Bus Substation Data" in line:
            start_idx = i + 1
            break

    if start_idx is None:
        raise ValueError(
            "Could not find 'Begin Bus Substation Data' section in GIC file"
        )

    data_lines = []
    for line in lines[start_idx:]:
        if line.strip().startswith("0"):
            break
        data_lines.append(line)

    bus_substation_map = {}
    for line in data_lines:
        parts = line.split()
        if len(parts) >= 2:
            bus_id = int(parts[0])
            substation_idx = int(parts[1])
            bus_substation_map[bus_id] = substation_idx

    bus_geo_df = pd.DataFrame(
        substation_df.loc[[v for v in bus_substation_map.values()], :]
    )
    bus_geo_df["sub_id"] = bus_geo_df.index
    bus_geo_df.index = list(bus_substation_map.keys())
    bus_geo_df.index.name = "BUS_I"
    return bus_geo_df


def create_bus_geo_df_from_aux(file_path: str | Path) -> pd.DataFrame:
    """Parse PowerWorld AUX and return bus-indexed geographic table with lat/lon.

    Expected mapping path:
    - `Substation` table provides substation `Number`, `Latitude`, `Longitude`
    - `Bus` table provides bus `Number` and `SubNumber`
    """

    aux_path = Path(file_path).expanduser().resolve()
    if not aux_path.exists():
        raise FileNotFoundError(f"AUX file not found: {aux_path}")

    lines = aux_path.read_text().splitlines()
    substation_df = _parse_aux_table(lines, "Substation")
    bus_df = _parse_aux_table(lines, "Bus")

    required_sub_cols = {"Number", "Latitude", "Longitude"}
    required_bus_cols = {"Number", "SubNumber"}
    missing_sub = sorted(required_sub_cols - set(substation_df.columns))
    missing_bus = sorted(required_bus_cols - set(bus_df.columns))
    if missing_sub:
        raise ValueError(
            f"AUX Substation table is missing required columns: {', '.join(missing_sub)}"
        )
    if missing_bus:
        raise ValueError(
            f"AUX Bus table is missing required columns: {', '.join(missing_bus)}"
        )

    sub = substation_df[["Number", "Name", "Latitude", "Longitude"]].copy()
    sub["Number"] = pd.to_numeric(sub["Number"], errors="coerce")
    sub["Latitude"] = pd.to_numeric(sub["Latitude"], errors="coerce")
    sub["Longitude"] = pd.to_numeric(sub["Longitude"], errors="coerce")
    sub = sub.dropna(subset=["Number", "Latitude", "Longitude"])
    sub["sub_id"] = sub["Number"].astype(int)

    bus = bus_df[["Number", "SubNumber"]].copy()
    bus["Number"] = pd.to_numeric(bus["Number"], errors="coerce")
    bus["SubNumber"] = pd.to_numeric(bus["SubNumber"], errors="coerce")
    bus = bus.dropna(subset=["Number", "SubNumber"])
    bus["BUS_I"] = bus["Number"].astype(int)
    bus["sub_id"] = bus["SubNumber"].astype(int)

    merged = bus.merge(
        sub[["sub_id", "Name", "Latitude", "Longitude"]],
        on="sub_id",
        how="left",
    )
    if merged[["Latitude", "Longitude"]].isna().any().any():
        missing_bus_ids = (
            merged[merged[["Latitude", "Longitude"]].isna().any(axis=1)]["BUS_I"]
            .astype(int)
            .tolist()
        )
        raise ValueError(
            "AUX bus->substation merge has missing coordinates for BUS_I "
            f"(first 20): {missing_bus_ids[:20]}"
        )

    bus_geo_df = merged.rename(
        columns={"Name": "name", "Latitude": "lat", "Longitude": "lon"}
    )[["BUS_I", "name", "sub_id", "lat", "lon"]]
    bus_geo_df = (
        bus_geo_df.drop_duplicates(subset=["BUS_I"]).set_index("BUS_I").sort_index()
    )
    return bus_geo_df


def create_bus_geo_df(file_path: str | Path) -> pd.DataFrame:
    """Parse supported geo source (.gic or .aux) into bus-indexed lat/lon."""

    geo_path = Path(file_path).expanduser().resolve()
    suffix = geo_path.suffix.lower()

    if suffix == ".gic":
        return create_bus_geo_df_from_gic(geo_path)
    if suffix == ".aux":
        return create_bus_geo_df_from_aux(geo_path)

    # Fallback for uncommon extensions: try GIC first, then AUX parser.
    try:
        return create_bus_geo_df_from_gic(geo_path)
    except Exception:
        return create_bus_geo_df_from_aux(geo_path)


def split_points(bus_df: pd.DataFrame, radius: float = 0.006) -> pd.DataFrame:
    """Return bus_df with colocated lat/lon points spread on small circles."""

    out_df = bus_df.copy()
    lat_name, lon_name = _infer_lat_lon_cols(out_df)

    unique_locs = (
        out_df[[lat_name, lon_name]]
        .groupby([lat_name, lon_name])
        .size()
        .reset_index(name="n")
    )

    for _, loc_row in unique_locs.iterrows():
        if int(loc_row["n"]) <= 1:
            continue

        lat = loc_row[lat_name]
        lon = loc_row[lon_name]
        colocated_idx = out_df[
            (out_df[lat_name] == lat) & (out_df[lon_name] == lon)
        ].index
        n_buses = len(colocated_idx)
        theta = 2 * np.pi / n_buses

        for k, bus_idx in enumerate(colocated_idx):
            x_new = lon + radius * np.cos(k * theta)
            y_new = lat + radius * np.sin(k * theta)
            out_df.loc[bus_idx, lat_name] = y_new
            out_df.loc[bus_idx, lon_name] = x_new

    return out_df


def attach_geo_to_case(
    case_data: MatpowerCaseData,
    geo_file_path: str | Path,
    split_colocated: bool = True,
    split_radius: float = 0.006,
    gen_fanout: bool = False,
    gen_fanout_radius: float = 0.003,
    gen_fanout_singles: bool = False,
) -> GeoPrepResult:
    """Attach geo lat/lon to bus/gen/branch tables, optionally applying split-points."""

    bus_geo_df = create_bus_geo_df(geo_file_path)

    bus_df = case_data.bus.copy()
    if "BUS_I" not in bus_df.columns:
        raise ValueError("case_data.bus must include BUS_I")

    bus_df["lat"] = bus_df["BUS_I"].map(bus_geo_df["lat"])
    bus_df["lon"] = bus_df["BUS_I"].map(bus_geo_df["lon"])

    if bus_df[["lat", "lon"]].isna().any().any():
        missing_ids = bus_df[bus_df[["lat", "lon"]].isna().any(axis=1)][
            "BUS_I"
        ].tolist()
        raise ValueError(
            f"Missing GIC coordinates for BUS_I ids (first 20): {missing_ids[:20]}"
        )

    n_before = int(bus_df[["lat", "lon"]].drop_duplicates().shape[0])
    if split_colocated:
        bus_df = split_points(bus_df, radius=split_radius)
    n_after = int(bus_df[["lat", "lon"]].drop_duplicates().shape[0])

    gen_df = case_data.gen.copy()
    bus_indexed = bus_df.set_index("BUS_I")
    gen_df["bus_lat"] = gen_df["GEN_BUS"].map(bus_indexed["lat"])
    gen_df["bus_lon"] = gen_df["GEN_BUS"].map(bus_indexed["lon"])
    gen_df["lat"] = gen_df["bus_lat"]
    gen_df["lon"] = gen_df["bus_lon"]
    if gen_fanout:
        gen_df = split_generator_points(
            gen_df, radius=gen_fanout_radius, fan_singles=gen_fanout_singles
        )

    branch_df = case_data.branch.copy()
    bus_lat = bus_df.set_index("BUS_I")["lat"]
    bus_lon = bus_df.set_index("BUS_I")["lon"]
    branch_df["from_lat"] = branch_df["F_BUS"].map(bus_lat)
    branch_df["from_lon"] = branch_df["F_BUS"].map(bus_lon)
    branch_df["to_lat"] = branch_df["T_BUS"].map(bus_lat)
    branch_df["to_lon"] = branch_df["T_BUS"].map(bus_lon)

    if all(col in branch_df.columns for col in ["PF", "QF", "PT", "QT", "RATE_A"]):
        branch_df["loading_pct"] = compute_branch_loading(branch_df)

    geo_case = MatpowerCaseData(
        case_name=case_data.case_name,
        source_path=case_data.source_path,
        bus=bus_df,
        gen=gen_df,
        branch=branch_df,
        genfuel=case_data.genfuel.copy() if case_data.genfuel is not None else None,
        gencost=case_data.gencost.copy() if case_data.gencost is not None else None,
        base_mva=case_data.base_mva,
        version=case_data.version,
    )

    return GeoPrepResult(
        case_data=geo_case,
        bus_geo_df=bus_geo_df,
        split_applied=split_colocated,
        n_unique_bus_locations_before=n_before,
        n_unique_bus_locations_after=n_after,
    )


def compute_branch_loading(branch_df: pd.DataFrame) -> pd.Series:
    """Compute branch loading (%) from PF/QF/PT/QT and RATE_A columns."""

    sf = np.sqrt(branch_df["PF"] ** 2 + branch_df["QF"] ** 2)
    st = np.sqrt(branch_df["PT"] ** 2 + branch_df["QT"] ** 2)
    smax = np.maximum(sf, st)
    loading_pct = 100 * smax / branch_df["RATE_A"]
    loading_pct = loading_pct.where(branch_df["RATE_A"] != 0, np.nan)
    return loading_pct


def lookup_fault_bus(
    bus_df: pd.DataFrame,
    fault_bus: int | str,
) -> pd.DataFrame:
    """Return the row(s) of *bus_df* matching *fault_bus*.

    - ``int``: exact match on ``BUS_I``.
    - ``str``: exact match on ``bus_name``, then case-insensitive substring fallback.

    Returns an empty DataFrame when no match is found.
    """
    if isinstance(fault_bus, int):
        return bus_df[bus_df["BUS_I"] == fault_bus]

    if "bus_name" not in bus_df.columns:
        print(
            "WARNING: bus_name column not available — cannot look up fault bus by name."
        )
        return pd.DataFrame()

    row = bus_df[bus_df["bus_name"] == fault_bus]
    if row.empty:
        row = bus_df[bus_df["bus_name"].str.contains(fault_bus, case=False, na=False)]
    return row


def plot_grid(
    case_data: MatpowerCaseData,
    map_style: str = "carto-positron",
    zoom: float | None = None,
    branch_bins: int = 10,
    n_hover_points: int = 10,
    show_loading: bool = True,
    show_gen_connectors: bool = False,
    marker_scale: float = 1.0,
) -> go.Figure:
    """Geo plot: branches, buses (PD), generators (fuel-colored when available).

    Parameters
    ----------
    show_loading:
        True  — color branches by loading_pct bins (viridis). Requires ``loading_pct``
                on case_data.branch (set by attach_geo_to_case when PF/QF/PT/QT/RATE_A
                are present). Falls back silently to gray when the column is missing.
        False — uniform gray branches regardless of whether loading_pct exists.
    """

    required_bus_cols = ["lat", "lon", "PD", "BUS_I"]
    missing_bus_cols = [c for c in required_bus_cols if c not in case_data.bus.columns]
    if missing_bus_cols:
        raise ValueError(f"Missing required bus columns: {missing_bus_cols}")

    required_branch_cols = [
        "from_lat",
        "from_lon",
        "to_lat",
        "to_lon",
        "F_BUS",
        "T_BUS",
    ]
    missing_branch_cols = [
        c for c in required_branch_cols if c not in case_data.branch.columns
    ]
    if missing_branch_cols:
        raise ValueError(f"Missing required branch columns: {missing_branch_cols}")

    required_gen_cols = ["lat", "lon", "GEN_BUS"]
    missing_gen_cols = [c for c in required_gen_cols if c not in case_data.gen.columns]
    if missing_gen_cols:
        raise ValueError(f"Missing required gen columns: {missing_gen_cols}")

    fig = go.Figure()

    has_loading = show_loading and "loading_pct" in case_data.branch.columns

    if has_loading:
        loading_pct = case_data.branch["loading_pct"]
        colorscale = _viridis_colorscale()

        for i in range(branch_bins):
            lo = i * 100.0 / branch_bins
            hi = (i + 1) * 100.0 / branch_bins
            idx = loading_pct[(loading_pct >= lo) & (loading_pct < hi)].index
            if len(idx) == 0:
                continue

            bin_df = case_data.branch.loc[idx]
            lat_vals = _branch_lat_array(bin_df, n_hover_points)
            lon_vals = _branch_lon_array(bin_df, n_hover_points)
            hover_vals = _branch_hover_array(bin_df, n_hover_points)
            frac = (i + 0.5) / branch_bins
            color = _sample_viridis(colorscale, frac)

            fig.add_trace(
                go.Scattermap(
                    mode="lines",
                    lat=lat_vals,
                    lon=lon_vals,
                    hovertext=hover_vals,
                    hoverinfo="text",
                    line=dict(color=color, width=2.5),
                    showlegend=False,
                )
            )

        # branches with nan loading (RATE_A == 0 or missing)
        nan_idx = loading_pct[loading_pct.isna()].index
        if len(nan_idx) > 0:
            nan_df = case_data.branch.loc[nan_idx]
            lat_vals = _branch_lat_array(nan_df, n_hover_points)
            lon_vals = _branch_lon_array(nan_df, n_hover_points)
            hover_vals = _branch_hover_array(nan_df, n_hover_points)
            fig.add_trace(
                go.Scattermap(
                    mode="lines",
                    lat=lat_vals,
                    lon=lon_vals,
                    hovertext=hover_vals,
                    hoverinfo="text",
                    line=dict(color="rgba(120, 120, 120, 0.4)", width=1.5),
                    name="branches (no rating)",
                    showlegend=True,
                )
            )

        # dummy colorbar trace for loading scale
        tickvals = [i / branch_bins for i in range(branch_bins + 1)]
        ticktext = [f"{int(v * 100)}%" for v in tickvals]
        dcolorsc = _discrete_colorscale(
            tickvals,
            [
                _sample_viridis(colorscale, (i + 0.5) / branch_bins)
                for i in range(branch_bins)
            ],
        )
        fig.add_trace(
            go.Scattermap(
                lat=[None],
                lon=[None],
                mode="markers",
                marker=dict(
                    colorscale=dcolorsc,
                    showscale=True,
                    cmin=0,
                    cmax=1,
                    colorbar=dict(
                        title=dict(text="Line<br>loading"),
                        thickness=15,
                        tickvals=tickvals,
                        ticktext=ticktext,
                        outlinewidth=0,
                        x=1.17,
                    ),
                ),
                hoverinfo="none",
                showlegend=False,
            )
        )
    else:
        # fallback: no loading data — uniform gray branches
        lat_vals = _branch_lat_array(case_data.branch, n_hover_points)
        lon_vals = _branch_lon_array(case_data.branch, n_hover_points)
        hover_vals = _branch_hover_array(case_data.branch, n_hover_points)
        fig.add_trace(
            go.Scattermap(
                mode="lines",
                lat=lat_vals,
                lon=lon_vals,
                hovertext=hover_vals,
                hoverinfo="text",
                line=dict(color="rgba(80, 80, 80, 0.45)", width=1.5),
                name="branches",
                showlegend=True,
            )
        )

    # buses
    bus_sizes = (
        (np.sqrt(case_data.bus["PD"].clip(lower=0.01).values) * 2.2 + 3) * marker_scale
    ).tolist()
    has_bus_name = "bus_name" in case_data.bus.columns
    has_json_bus = "json_bus_number" in case_data.bus.columns
    bus_hover = [
        (
            f"<b>BUS_I: {row.BUS_I}</b><br>"
            + (f"{row.bus_name}<br>" if has_bus_name else "")
            + f"PD: {row.PD:.2f} MW<br>QD: {row.QD:.2f} MVAr"
            + (f"<br>json: number={row.json_bus_number}" if has_json_bus else "")
        )
        for row in case_data.bus.itertuples()
    ]
    fig.add_trace(
        go.Scattermap(
            name="buses (PD)",
            mode="markers",
            lat=case_data.bus["lat"].tolist(),
            lon=case_data.bus["lon"].tolist(),
            hovertext=bus_hover,
            hoverinfo="text",
            marker=dict(
                color=case_data.bus["PD"].values,
                colorscale="Blues",
                colorbar=dict(
                    title=dict(text="Load<br>(PD in MW)"),
                    x=1.02,
                    thickness=15,
                    outlinewidth=0,
                ),
                size=bus_sizes,
                opacity=0.85,
            ),
            showlegend=True,
        )
    )

    # optional generator fan-out connectors
    if show_gen_connectors:
        gen_conn_lat, gen_conn_lon, gen_conn_hover = _build_gen_connector_line_arrays(
            case_data.gen
        )
        if gen_conn_lat:
            fig.add_trace(
                go.Scattermap(
                    mode="lines",
                    lat=gen_conn_lat,
                    lon=gen_conn_lon,
                    hoverinfo="none",
                    line=dict(color="rgba(60, 60, 60, 0.35)", width=1.0),
                    name="fake gen connections",
                    showlegend=True,
                )
            )

    # generators
    palette = [
        "#1f77b4",
        "#d62728",
        "#2ca02c",
        "#ff7f0e",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
    ]
    fuel_color_map = {
        "coal": "#222222",
        "nuclear": "#a50026",
        "hydro": "#4393c3",
        "wind": "#4dac26",
        "solar": "#f4a582",
        "ng": "#fd8d3c",
        "gas": "#fd8d3c",
        "oil": "#8c6d31",
        "geothermal": "#762a83",
        "biomass": "#74c476",
        "other": "#969696",
        "sync_cond": "#d4b9da",
    }
    if case_data.genfuel is not None and len(case_data.genfuel) == len(case_data.gen):
        gen_plot_df = case_data.gen.copy()
        gen_plot_df["genfuel"] = case_data.genfuel.values
        palette_fuels = [
            f
            for f in sorted(gen_plot_df["genfuel"].unique())
            if f not in fuel_color_map
        ]
        for i, fuel in enumerate(sorted(gen_plot_df["genfuel"].unique())):
            gdf = gen_plot_df[gen_plot_df["genfuel"] == fuel]
            color = (
                fuel_color_map[fuel]
                if fuel in fuel_color_map
                else palette[palette_fuels.index(fuel) % len(palette)]
            )
            gen_sizes = (
                (np.sqrt(gdf["PG"].clip(lower=0.01).values) * 1.8 + 4) * marker_scale
                if "PG" in gdf.columns
                else np.full(len(gdf), 6.0 * marker_scale)
            )
            _has_json_gen = "json_gen_id" in gdf.columns
            _has_json_gen_bus = "json_gen_bus_number" in gdf.columns
            gen_hover = [
                f"<b>GEN_BUS: {row.GEN_BUS}</b><br>"
                f"gen_row: {getattr(row, 'gen_row', 'n/a')}<br>"
                f"PG: {row.PG:.2f} MW<br>QG: {row.QG:.2f} MVAr<br>Fuel: {fuel}"
                + (
                    (
                        f"<br>json: number={getattr(row, 'json_gen_bus_number', row.GEN_BUS)}, device_id={row.json_gen_id}"
                        if row.json_gen_id is not None
                        else "<br><b>\u26a0 offline \u2014 not in JSON (GEN_STATUS=0)</b>"
                    )
                    if _has_json_gen
                    else ""
                )
                for row in gdf.itertuples()
            ]
            fig.add_trace(
                go.Scattermap(
                    name=f"gen: {fuel}",
                    mode="markers",
                    lat=gdf["lat"].tolist(),
                    lon=gdf["lon"].tolist(),
                    hovertext=gen_hover,
                    hoverinfo="text",
                    marker=dict(size=gen_sizes.tolist(), color=color, opacity=0.9),
                    showlegend=True,
                )
            )
    else:
        fig.add_trace(
            go.Scattermap(
                name="generators",
                mode="markers",
                lat=case_data.gen["lat"].tolist(),
                lon=case_data.gen["lon"].tolist(),
                hovertext=[
                    (
                        (
                            f"GEN_BUS: {row.GEN_BUS}<br>gen_row: {row.gen_row}"
                            if "gen_row" in case_data.gen.columns
                            else f"GEN_BUS: {row.GEN_BUS}"
                        )
                        + (
                            (
                                f"<br>json: number={getattr(row, 'json_gen_bus_number', row.GEN_BUS)}, device_id={row.json_gen_id}"
                                if row.json_gen_id is not None
                                else "<br><b>\u26a0 offline \u2014 not in JSON (GEN_STATUS=0)</b>"
                            )
                            if "json_gen_id" in case_data.gen.columns
                            else ""
                        )
                    )
                    for row in case_data.gen.itertuples()
                ],
                hoverinfo="text",
                marker=dict(size=6 * marker_scale, color="#d62728", opacity=0.9),
                showlegend=True,
            )
        )

    center_lon = float(case_data.bus["lon"].mean())
    center_lat = float(case_data.bus["lat"].mean())
    if zoom is None:
        zoom = 7 if len(case_data.bus) <= 400 else 5

    fig.update_layout(
        title=f"{case_data.case_name}: grid map"
        + (" (loading)" if has_loading else ""),
        map={
            "style": map_style,
            "center": {"lon": center_lon, "lat": center_lat},
            "zoom": zoom,
        },
        legend={"orientation": "h", "y": -0.15, "x": 0, "itemsizing": "constant"},
        margin=dict(r=120),
        height=650,
        width=780,
    )
    return fig


def plot_first_geo(
    case_data: MatpowerCaseData,
    map_style: str = "carto-positron",
    zoom: float | None = None,
    show_gen_connectors: bool = False,
) -> go.Figure:
    """Alias for plot_grid(show_loading=False). Kept for backward compatibility."""

    return plot_grid(
        case_data,
        map_style=map_style,
        zoom=zoom,
        show_loading=False,
        show_gen_connectors=show_gen_connectors,
    )


def _infer_lat_lon_cols(df: pd.DataFrame) -> tuple[str, str]:
    if "lat" in df.columns and "lon" in df.columns:
        return "lat", "lon"
    if "latitude" in df.columns and "longitude" in df.columns:
        return "latitude", "longitude"
    if "Latitude" in df.columns and "Longitude" in df.columns:
        return "Latitude", "Longitude"
    raise ValueError("Could not find latitude/longitude columns in dataframe")


def _parse_aux_table(lines: list[str], table_name: str) -> pd.DataFrame:
    """Parse a PowerWorld AUX table into a DataFrame.

    Supports both canonical table headers (`Substation (...)`) and a lightweight
    `DATA (Substation)` marker style when present.
    """

    table_start_idx = _find_aux_table_header(lines, table_name)
    header_line = lines[table_start_idx]

    header_buffer = [header_line]
    idx = table_start_idx + 1
    while idx < len(lines) and ")" not in " ".join(header_buffer):
        header_buffer.append(lines[idx])
        idx += 1

    header_text = " ".join(header_buffer)
    col_match = re.search(r"\((.*)\)", header_text)
    if col_match is None:
        raise ValueError(f"Could not parse column header for AUX table '{table_name}'")

    columns = [c.strip() for c in col_match.group(1).split(",") if c.strip()]
    if not columns:
        raise ValueError(f"AUX table '{table_name}' has no parsed columns")

    while idx < len(lines) and "{" not in lines[idx]:
        idx += 1
    if idx >= len(lines):
        raise ValueError(f"AUX table '{table_name}' is missing opening '{{'")
    idx += 1

    rows: list[list[str]] = []
    while idx < len(lines):
        line = lines[idx].strip()
        if not line or line.startswith("//"):
            idx += 1
            continue
        if line.startswith("}"):
            break

        tokens = shlex.split(line)
        if len(tokens) < len(columns):
            tokens.extend([""] * (len(columns) - len(tokens)))
        elif len(tokens) > len(columns):
            tokens = tokens[: len(columns)]
        rows.append(tokens)
        idx += 1

    if not rows:
        raise ValueError(f"AUX table '{table_name}' has no data rows")

    return pd.DataFrame(rows, columns=columns)


def _find_aux_table_header(lines: list[str], table_name: str) -> int:
    canonical_regex = re.compile(rf"^\s*{re.escape(table_name)}\s*\(")
    data_regex = re.compile(
        rf"^\s*DATA\s*\(\s*{re.escape(table_name)}\s*\)", re.IGNORECASE
    )

    for idx, line in enumerate(lines):
        if canonical_regex.search(line) or data_regex.search(line):
            return idx
    raise ValueError(
        f"Could not find AUX table '{table_name}' (expected '{table_name} (...)' or 'DATA ({table_name})')"
    )


def split_generator_points(
    gen_df: pd.DataFrame,
    radius: float = 0.003,
    fan_singles: bool = False,
) -> pd.DataFrame:
    """Spread generators on the same bus around the bus center for visibility."""

    out_df = gen_df.copy()
    required_cols = ["GEN_BUS", "lat", "lon"]
    missing_cols = [c for c in required_cols if c not in out_df.columns]
    if missing_cols:
        raise ValueError(
            f"Generator fan-out requires columns: {', '.join(required_cols)}"
        )

    for _, grp in out_df.groupby("GEN_BUS", sort=False):
        row_idx = grp.index.tolist()
        n_gen = len(row_idx)

        center_lat = float(out_df.loc[row_idx[0], "lat"])
        center_lon = float(out_df.loc[row_idx[0], "lon"])
        # single generator: optionally offset east so bus marker underneath is visible
        if n_gen == 1:
            if fan_singles:
                out_df.loc[row_idx[0], "lon"] = center_lon + radius
            continue

        theta = 2 * np.pi / n_gen
        for k, idx in enumerate(row_idx):
            out_df.loc[idx, "lon"] = center_lon + radius * np.cos(k * theta)
            out_df.loc[idx, "lat"] = center_lat + radius * np.sin(k * theta)

    return out_df


def _viridis_colorscale() -> list[tuple[float, str]]:
    """Return viridis colorscale as [(frac, 'rgb(...)'), ...] list."""

    import plotly.colors as pc

    colors = pc.sequential.Viridis  # list of hex/rgb strings
    n = len(colors)
    return [(i / (n - 1), c) for i, c in enumerate(colors)]


def _sample_viridis(colorscale: list[list], frac: float) -> str:
    """Sample a color from a plotly colorscale at position frac in [0, 1]."""

    if frac <= 0:
        return colorscale[0][1]
    if frac >= 1:
        return colorscale[-1][1]
    for i in range(len(colorscale) - 1):
        lo_frac, lo_col = colorscale[i]
        hi_frac, hi_col = colorscale[i + 1]
        if lo_frac <= frac <= hi_frac:
            # just snap to nearest stop for simplicity
            t = (frac - lo_frac) / (hi_frac - lo_frac)
            return lo_col if t < 0.5 else hi_col
    return colorscale[-1][1]


def _discrete_colorscale(bvals: list[float], colors: list[str]) -> list[list]:
    """Build a plotly discrete colorscale from bin boundaries and per-bin colors."""

    n = len(colors)
    nvals = [(v - bvals[0]) / (bvals[-1] - bvals[0]) for v in bvals]
    dcolorscale: list[list] = []
    for k in range(n):
        dcolorscale.extend([[nvals[k], colors[k]], [nvals[k + 1], colors[k]]])
    return dcolorscale


def _branch_lat_array(branch_df: pd.DataFrame, n_points: int) -> list:
    """Interpolate n_points between from_lat and to_lat for each branch, NaN-separated.

    Parallel branches (same F_BUS/T_BUS) are deduplicated — only the first occurrence
    is rendered, matching the grouping in _branch_hover_array.
    """
    deduped = branch_df.drop_duplicates(subset=["F_BUS", "T_BUS"], keep="first")
    arr = np.pad(
        np.linspace(
            deduped["from_lat"].values, deduped["to_lat"].values, n_points, axis=1
        ),
        [(0, 0), (0, 1)],
        constant_values=np.nan,
    )
    return arr.reshape(-1).tolist()


def _branch_lon_array(branch_df: pd.DataFrame, n_points: int) -> list:
    """Interpolate n_points between from_lon and to_lon for each branch, NaN-separated.

    Parallel branches (same F_BUS/T_BUS) are deduplicated — only the first occurrence
    is rendered, matching the grouping in _branch_hover_array.
    """
    deduped = branch_df.drop_duplicates(subset=["F_BUS", "T_BUS"], keep="first")
    arr = np.pad(
        np.linspace(
            deduped["from_lon"].values, deduped["to_lon"].values, n_points, axis=1
        ),
        [(0, 0), (0, 1)],
        constant_values=np.nan,
    )
    return arr.reshape(-1).tolist()


def _branch_hover_detail(
    row: object, has_loading: bool, has_branch_row: bool, has_json_branch: bool = False
) -> str:
    """Build the per-branch detail line (no F_BUS/T_BUS header — used inside groups)."""
    parts = []
    if has_branch_row:
        parts.append(f"branch_row: {row.branch_row}")
    if has_json_branch:
        parts.append(f"json: device_id={row.json_branch_id}")
    if has_loading and not pd.isna(row.loading_pct):
        parts.append(f"Loading: {row.loading_pct:.1f}%")
        pf = getattr(row, "PF", None)
        if pf is not None:
            parts.append(f"PF: {pf:.1f} MW")
        rate_a = getattr(row, "RATE_A", None)
        if rate_a is not None:
            parts.append(f"RATE_A: {rate_a:.1f} MVA")
    return "<br>".join(parts) if parts else ""


def _branch_hover_array(branch_df: pd.DataFrame, n_points: int) -> list:
    """Build hover text repeated n_points times per branch, with empty string separator.

    Parallel branches (same F_BUS/T_BUS) are grouped under one header line that
    states the count, followed by per-branch detail lines separated by ----.
    """
    has_loading = "loading_pct" in branch_df.columns
    has_branch_row = "branch_row" in branch_df.columns
    has_json_branch = "json_branch_id" in branch_df.columns

    groups: dict[tuple, list[str]] = {}
    order: list[tuple] = []
    for row in branch_df.itertuples():
        key = (int(row.F_BUS), int(row.T_BUS))
        detail = _branch_hover_detail(row, has_loading, has_branch_row, has_json_branch)
        if key not in groups:
            groups[key] = [detail]
            order.append(key)
        else:
            groups[key].append(detail)

    result = []
    for key in order:
        f_bus, t_bus = key
        details = groups[key]
        n = len(details)
        count = f" ({n} parallel)" if n > 1 else ""
        header = f"<b>F_BUS: {f_bus} \u2192 T_BUS: {t_bus}{count}</b>"
        body = "<br>---<br>".join(d for d in details if d)
        hover = f"{header}<br>{body}" if body else header
        result.extend([hover] * n_points + [""])
    return result


def _build_branch_line_arrays(
    branch_df: pd.DataFrame,
) -> tuple[list[float | None], list[float | None], list[str]]:
    lat_vals: list[float | None] = []
    lon_vals: list[float | None] = []
    hover_vals: list[str] = []
    has_loading = "loading_pct" in branch_df.columns
    has_branch_row = "branch_row" in branch_df.columns
    has_json_branch = "json_branch_id" in branch_df.columns

    for row in branch_df.itertuples():
        lat_vals.extend([row.from_lat, row.to_lat, None])
        lon_vals.extend([row.from_lon, row.to_lon, None])

        branch_id = f"<br>branch_row: {row.branch_row}" if has_branch_row else ""
        json_bid = (
            f"<br>json: device_id={row.json_branch_id}" if has_json_branch else ""
        )
        if has_loading and not pd.isna(row.loading_pct):
            rate_a = getattr(row, "RATE_A", None)
            pf = getattr(row, "PF", None)
            rate_str = f"<br>RATE_A: {rate_a:.1f} MVA" if rate_a is not None else ""
            pf_str = f"<br>PF: {pf:.1f} MW" if pf is not None else ""
            htxt = (
                f"<b>F_BUS: {row.F_BUS} \u2192 T_BUS: {row.T_BUS}</b>"
                f"{branch_id}{json_bid}<br>"
                f"Loading: {row.loading_pct:.1f}%"
                f"{pf_str}{rate_str}"
            )
        else:
            htxt = f"<b>F_BUS: {row.F_BUS} \u2192 T_BUS: {row.T_BUS}</b>{branch_id}{json_bid}"
        hover_vals.extend([htxt, htxt, ""])

    return lat_vals, lon_vals, hover_vals


def _build_gen_connector_line_arrays(
    gen_df: pd.DataFrame,
) -> tuple[list[float | None], list[float | None], list[str]]:
    required_cols = ["bus_lat", "bus_lon", "lat", "lon", "GEN_BUS"]
    if any(c not in gen_df.columns for c in required_cols):
        return [], [], []

    lat_vals: list[float | None] = []
    lon_vals: list[float | None] = []
    hover_vals: list[str] = []

    for row in gen_df.itertuples():
        if (
            pd.isna(row.bus_lat)
            or pd.isna(row.bus_lon)
            or pd.isna(row.lat)
            or pd.isna(row.lon)
        ):
            continue

        lat_vals.extend([row.bus_lat, row.lat, None])
        lon_vals.extend([row.bus_lon, row.lon, None])

        gen_row = getattr(row, "gen_row", "n/a")
        htxt = f"<b>GEN_BUS: {row.GEN_BUS}</b><br>gen_row: {gen_row}"
        hover_vals.extend([htxt, htxt, ""])

    return lat_vals, lon_vals, hover_vals


def _parse_m_cell_array(m_text: str, attr_name: str) -> list[str]:
    """Extract a MATPOWER curly-brace string cell array from .m file text.

    Handles patterns like::

        mpc.genfuel = {
            'coal';
            'wind';
        };

    Returns a list of strings (quotes and semicolons stripped), or [] if not found.
    """

    pattern = re.compile(
        r"mpc\." + re.escape(attr_name) + r"\s*=\s*\{([^}]*)\}\s*;",
        re.DOTALL,
    )
    m = pattern.search(m_text)
    if m is None:
        return []

    entries = []
    for line in m.group(1).splitlines():
        line = line.strip().rstrip(";").strip().strip("'\"")
        if line and not line.startswith("%"):
            entries.append(line)
    return entries


def _copy_df_attr(case: CaseFrames, attr_name: str) -> pd.DataFrame | None:
    value = getattr(case, attr_name, None)
    if value is None:
        return None
    if not isinstance(value, pd.DataFrame):
        return None
    return value.copy()


def _to_int_series(series: pd.Series, column_name: str) -> pd.Series:
    try:
        return pd.to_numeric(series, errors="raise").astype(int)
    except Exception as exc:
        raise ValueError(
            f"Column '{column_name}' could not be converted to int"
        ) from exc


def _coerce_optional_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _coerce_optional_str(value: Any) -> str | None:
    if value is None:
        return None
    return str(value)
