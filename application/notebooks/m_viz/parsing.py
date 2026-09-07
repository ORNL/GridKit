"""MATPOWER .m case loading and normalization for m_viz."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Any

import pandas as pd

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
