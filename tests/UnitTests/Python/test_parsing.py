"""Unit tests for m_viz.parsing (no real MATPOWER .m file required)."""

import pandas as pd
import pytest

from m_viz.parsing import (
    MatpowerCaseData,
    _parse_m_cell_array,
    normalize_bus_df,
    normalize_gen_df,
    normalize_branch_df,
    summarize_case,
    validate_case_tables,
)


def test_normalize_bus_df_sorts_by_bus_i_and_sets_index():
    df = pd.DataFrame({"BUS_I": [3, 1, 2], "PD": [10.0, 20.0, 30.0]})

    out = normalize_bus_df(df)

    assert out.index.tolist() == [1, 2, 3]
    assert out["BUS_I"].tolist() == [1, 2, 3]
    assert out.loc[1, "PD"] == 20.0


def test_normalize_bus_df_requires_bus_i():
    df = pd.DataFrame({"PD": [10.0]})

    with pytest.raises(ValueError, match="BUS_I"):
        normalize_bus_df(df)


def test_normalize_gen_df_renames_row_index_column():
    df = pd.DataFrame({"GEN_BUS": [1, 2], "PG": [10.0, 20.0]})

    out = normalize_gen_df(df)

    assert "gen_row" in out.columns
    assert out["GEN_BUS"].tolist() == [1, 2]
    assert out.index.name == "gen_row"


def test_normalize_gen_df_requires_gen_bus():
    df = pd.DataFrame({"PG": [10.0]})

    with pytest.raises(ValueError, match="GEN_BUS"):
        normalize_gen_df(df)


def test_normalize_branch_df_requires_from_and_to_bus():
    df = pd.DataFrame({"F_BUS": [1, 2]})

    with pytest.raises(ValueError, match="T_BUS"):
        normalize_branch_df(df)


def test_normalize_branch_df_sets_branch_row_index():
    df = pd.DataFrame({"F_BUS": [1], "T_BUS": [2]})

    out = normalize_branch_df(df)

    assert "branch_row" in out.columns
    assert out.index.name == "branch_row"


def test_parse_m_cell_array_extracts_entries():
    m_text = """
    mpc.genfuel = {
        'coal';
        'wind';
    };
    """

    assert _parse_m_cell_array(m_text, "genfuel") == ["coal", "wind"]


def test_parse_m_cell_array_missing_returns_empty():
    assert _parse_m_cell_array("mpc.bus = [1 2 3];", "genfuel") == []


def _make_case_data():
    bus = normalize_bus_df(pd.DataFrame({"BUS_I": [1, 2]}))
    gen = normalize_gen_df(pd.DataFrame({"GEN_BUS": [1, 1]}))
    branch = normalize_branch_df(pd.DataFrame({"F_BUS": [1], "T_BUS": [2]}))
    return MatpowerCaseData(
        case_name="test_case",
        source_path="/tmp/test_case.m",
        bus=bus,
        gen=gen,
        branch=branch,
    )


def test_summarize_case_reports_counts():
    summary = summarize_case(_make_case_data())

    assert summary["n_buses"] == 2
    assert summary["n_generators"] == 2
    assert summary["n_branches"] == 1
    assert summary["n_unique_gen_buses"] == 1


def test_validate_case_tables_flags_missing_gen_bus():
    case_data = _make_case_data()
    case_data.gen["GEN_BUS"] = [1, 99]  # 99 does not exist in bus table

    report = validate_case_tables(case_data)

    assert report["gen_buses_are_valid"] is False
    assert report["missing_gen_bus_ids"] == [99]
