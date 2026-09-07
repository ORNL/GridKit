"""Geographic coordinate loading and attachment for m_viz."""

from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
from pathlib import Path
import re
import shlex

import numpy as np
import pandas as pd

from .parsing import MatpowerCaseData


@dataclass
class GeoPrepResult:
    """Result bundle for geo enrichment and split-point preprocessing."""

    case_data: MatpowerCaseData
    bus_geo_df: pd.DataFrame
    split_applied: bool
    n_unique_bus_locations_before: int
    n_unique_bus_locations_after: int


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
