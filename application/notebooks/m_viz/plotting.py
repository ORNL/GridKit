"""Plotly geo-map rendering for m_viz."""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from .parsing import MatpowerCaseData


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
                        else "<br><b>\u26a0 offline — not in JSON (GEN_STATUS=0)</b>"
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
                                else "<br><b>\u26a0 offline — not in JSON (GEN_STATUS=0)</b>"
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
