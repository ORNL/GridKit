"""Unit tests for m_viz.geo (no real .gic/.AUX file required)."""

import pandas as pd
import pytest

from m_viz.geo import (
    compute_branch_loading,
    lookup_fault_bus,
    split_generator_points,
    split_points,
)


def test_split_points_spreads_colocated_buses():
    df = pd.DataFrame({"lat": [10.0, 10.0, 20.0], "lon": [-80.0, -80.0, -90.0]})

    out = split_points(df, radius=0.01)

    # unique location is untouched
    assert out.loc[2, "lat"] == 20.0
    assert out.loc[2, "lon"] == -90.0

    # colocated pair is spread apart in longitude, latitude unchanged
    lons = sorted(out.loc[[0, 1], "lon"].tolist())
    assert lons == pytest.approx([-80.0 - 0.01, -80.0 + 0.01])
    assert out.loc[0, "lat"] == pytest.approx(10.0, abs=1e-9)
    assert out.loc[1, "lat"] == pytest.approx(10.0, abs=1e-9)


def test_split_points_missing_lat_lon_raises():
    with pytest.raises(ValueError):
        split_points(pd.DataFrame({"x": [1], "y": [2]}))


def test_split_generator_points_leaves_single_gen_in_place():
    gen_df = pd.DataFrame({"GEN_BUS": [1], "lat": [10.0], "lon": [-80.0]})

    out = split_generator_points(gen_df, radius=0.01, fan_singles=False)

    assert out.loc[0, "lat"] == 10.0
    assert out.loc[0, "lon"] == -80.0


def test_split_generator_points_fans_out_multi_gen_bus():
    gen_df = pd.DataFrame({"GEN_BUS": [1, 1], "lat": [10.0, 10.0], "lon": [-80.0, -80.0]})

    out = split_generator_points(gen_df, radius=0.01, fan_singles=False)

    lons = sorted(out["lon"].tolist())
    assert lons == pytest.approx([-80.0 - 0.01, -80.0 + 0.01])


def test_compute_branch_loading_percent_of_rate_a():
    branch_df = pd.DataFrame(
        {
            "PF": [3.0, 1.0],
            "QF": [4.0, 0.0],
            "PT": [0.0, 0.0],
            "QT": [0.0, 0.0],
            "RATE_A": [10.0, 0.0],
        }
    )

    loading = compute_branch_loading(branch_df)

    assert loading.iloc[0] == pytest.approx(50.0)
    assert pd.isna(loading.iloc[1])  # RATE_A == 0 -> NaN, not divide-by-zero


def test_lookup_fault_bus_by_bus_i():
    bus_df = pd.DataFrame({"BUS_I": [1, 2], "bus_name": ["ALOHA138", "WAENA138"]})

    row = lookup_fault_bus(bus_df, 2)

    assert row["bus_name"].iloc[0] == "WAENA138"


def test_lookup_fault_bus_by_name_case_insensitive_substring():
    bus_df = pd.DataFrame({"BUS_I": [1, 2], "bus_name": ["ALOHA138", "WAENA138"]})

    row = lookup_fault_bus(bus_df, "aloha")

    assert row["BUS_I"].iloc[0] == 1


def test_lookup_fault_bus_no_match_returns_empty():
    bus_df = pd.DataFrame({"BUS_I": [1], "bus_name": ["ALOHA138"]})

    row = lookup_fault_bus(bus_df, "nonexistent")

    assert row.empty
