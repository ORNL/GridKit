"""Utilities for running and post-processing GridKit solve_pf power flow cases."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# grid3bus output parser
# ---------------------------------------------------------------------------


def parse_grid3bus_output(stdout: str) -> pd.DataFrame:
    """Parse theta/V lines from grid3bus stdout into a small DataFrame."""
    rows = []
    current_case = None
    for line in stdout.splitlines():
        ll = line.lower()
        if "monolithic" in ll:
            current_case = "monolithic"
        elif "via a parser" in ll:
            current_case = "parser"
        elif "manually" in ll:
            current_case = "hardwired"
        m = re.match(r"\s*theta(\d+)\s*=\s*([-\d.]+)\s*deg", line)
        if m:
            rows.append(
                {
                    "case": current_case,
                    "var": f"theta{m.group(1)} (deg)",
                    "value": float(m.group(2)),
                }
            )
        m = re.match(r"\s*V(\d+)\s*=\s*([-\d.]+)\s*p\.u\.", line)
        if m:
            rows.append(
                {
                    "case": current_case,
                    "var": f"V{m.group(1)} (p.u.)",
                    "value": float(m.group(2)),
                }
            )
    return pd.DataFrame(rows).pivot(index="var", columns="case", values="value")


# ---------------------------------------------------------------------------
# solve_pf runners
# ---------------------------------------------------------------------------


def _parse_pf_stdout(stdout: str) -> pd.DataFrame:
    """Parse bus lines from solve_pf stdout into a DataFrame."""
    rows = []
    for line in stdout.splitlines():
        m = re.match(
            r"bus\s+(\d+)\s+V=([\d.eE+\-]+)\s+theta_deg=([\d.eE+\-]+)\s+type=(\d+)",
            line.strip(),
        )
        if m:
            rows.append(
                {
                    "bus_i": int(m.group(1)),
                    "V_pu": float(m.group(2)),
                    "theta_deg": float(m.group(3)),
                    "type": int(m.group(4)),
                }
            )
    return pd.DataFrame(rows)


def _extract_nni(stdout: str) -> str | None:
    """Extract KINSOL nonlinear iteration count from KINPrintAllStats stdout."""
    m = re.search(r"Nonlinear iters\s*=\s*(\d+)", stdout)
    return m.group(1) if m else None


def _run(cmd: list[str]) -> tuple[pd.DataFrame, str, int]:
    """Run solve_pf command, parse bus output, inject nni into stderr."""
    r = subprocess.run(cmd, capture_output=True, text=True)
    nni = _extract_nni(r.stdout)
    stderr = r.stderr
    if nni is not None:
        stderr += f"nni={nni}\n"
    return _parse_pf_stdout(r.stdout), stderr, r.returncode


def run_solve_pf(solve_pf_bin: Path, m_path: Path) -> tuple[pd.DataFrame, str, int]:
    """Run solve_pf on a .m file. Returns (bus_df, stderr_text, return_code)."""
    return _run([str(solve_pf_bin), str(m_path)])


def run_solve_pf_flat(
    solve_pf_bin: Path, m_path: Path
) -> tuple[pd.DataFrame, str, int]:
    """Run solve_pf with --flat-start (Vm=1, Va=0). Returns (bus_df, stderr, rc)."""
    return _run([str(solve_pf_bin), str(m_path), "--flat-start"])


def run_solve_pf_out(
    solve_pf_bin: Path, m_path: Path, out_m_path: Path
) -> tuple[pd.DataFrame, str, int]:
    """Run solve_pf with --output-m. Returns (bus_df, stderr, return_code)."""
    return _run([str(solve_pf_bin), str(m_path), "--output-m", str(out_m_path)])


def run_solve_pf_out_flat(
    solve_pf_bin: Path, m_path: Path, out_m_path: Path
) -> tuple[pd.DataFrame, str, int]:
    """Run solve_pf with --flat-start and --output-m. Returns (bus_df, stderr, rc)."""
    return _run(
        [str(solve_pf_bin), str(m_path), "--flat-start", "--output-m", str(out_m_path)]
    )


# ---------------------------------------------------------------------------
# result display helpers
# ---------------------------------------------------------------------------


def pf_summary(case_label: str, df: pd.DataFrame, rc: int, stderr: str) -> None:
    """Print convergence status, ||f||, nni, and voltage range/violations."""
    fnorm_m = re.search(r"Residual 2-norm:\s*([\d.eE+\-]+)", stderr)
    fnorm_str = fnorm_m.group(1) if fnorm_m else "n/a"
    nni_m = re.search(r"nni=(\d+)", stderr)
    nni_str = f"   nni={nni_m.group(1)}" if nni_m else ""
    converged = rc == 0
    print(f"{case_label}")
    print(
        f"  {'CONVERGED' if converged else 'NOT CONVERGED'}"
        f"   rc={rc}   ||f||={fnorm_str}   buses={len(df)}{nni_str}"
    )
    if converged and not df.empty:
        viol = df[(df.V_pu < 0.95) | (df.V_pu > 1.05)]
        print(
            f"  V range: [{df.V_pu.min():.4f}, {df.V_pu.max():.4f}] pu"
            f"   violations (V<0.95 or V>1.05): {len(viol)}"
        )
    print()


def diff_vs_base(
    label: str, perturbed_df: pd.DataFrame, base_df: pd.DataFrame
) -> pd.DataFrame:
    """Compare a perturbed solved solution against the base-case solved solution.

    Prints max/mean |delta V| and |delta theta|, returns a merged DataFrame.
    Warns if both dV and dTheta are near zero (likely stagnation at warm-start).
    """
    cmp = perturbed_df.merge(
        base_df[["bus_i", "V_pu", "theta_deg"]].rename(
            columns={"V_pu": "V_base", "theta_deg": "theta_base"}
        ),
        on="bus_i",
    )
    cmp["dV"] = (cmp.V_pu - cmp.V_base).abs()
    cmp["dTheta"] = (cmp.theta_deg - cmp.theta_base).abs()
    print(f"  vs base-case solution ({label}):")
    print(f"    max |dV|      = {cmp.dV.max():.6f} pu")
    print(f"    mean |dV|     = {cmp.dV.mean():.6f} pu")
    print(f"    max |dTheta|  = {cmp.dTheta.max():.6f} deg")
    print(f"    mean |dTheta| = {cmp.dTheta.mean():.6f} deg")
    if cmp.dV.max() < 1e-4 and cmp.dTheta.max() < 1e-4:
        print("    WARNING: solution is essentially identical to base case.")
        print("             Solver likely stagnated at warm-start point.")
    print()
    return cmp


def parse_raw_m_bus(m_path: Path) -> pd.DataFrame:
    """Read Vm/Va directly from an mpc.bus block (no solve).

    Works on any MATPOWER .m file: the raw/unsolved case (Vm/Va = flat or the
    original TAMU/PowerWorld reference solution embedded by the data provider),
    or an already-"solved" .m (e.g. GridKit's `solve_pf --output-m` / PowerModels'
    `--output-m` output) - both are plain MATPOWER files with the same mpc.bus
    column layout, so this one parser covers both use cases. Returns columns
    matching `_parse_pf_stdout`'s bus_df (bus_i, V_pu, theta_deg, type) so it can
    be used interchangeably with solve_pf/pm_solve.jl stdout output in
    `diff_vs_base`/`pm_summary`/`pf_summary`.
    """
    content = m_path.read_text()
    rows = []
    for row in _parse_block(content, "bus"):
        cols = row.split()
        rows.append(
            {
                "bus_i": int(float(cols[0])),
                "type": int(float(cols[1])),
                "V_pu": float(cols[7]),
                "theta_deg": float(cols[8]),
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# .m file block manipulation helpers
# ---------------------------------------------------------------------------


def _replace_block(content: str, block_name: str, new_block: str) -> str:
    """Replace mpc.<block_name> = [ ... ]; in content with new_block rows."""
    pat = re.compile(
        rf"(mpc\.{block_name}\s*=\s*\[)(.*?)(\];)",
        re.DOTALL,
    )
    return pat.sub(lambda m: m.group(1) + "\n" + new_block + m.group(3), content)


def _parse_block(content: str, block_name: str) -> list[str]:
    """Return data lines from mpc.<block_name> block (comments stripped)."""
    pat = re.compile(rf"mpc\.{block_name}\s*=\s*\[(.*?)\];", re.DOTALL)
    m = pat.search(content)
    if not m:
        raise ValueError(f"block mpc.{block_name} not found")
    lines = []
    for line in m.group(1).splitlines():
        stripped = line.split("%")[0].strip().rstrip(";").strip()
        if stripped:
            lines.append(stripped)
    return lines


# ---------------------------------------------------------------------------
# .m perturbation generators
# ---------------------------------------------------------------------------


def make_perturbed_load_m(
    src: Path, dst: Path, pct: float = 0.05, seed: int = 42
) -> None:
    """Scale Pd (col 3) and Qd (col 4) of every bus by (1 + Uniform(-pct, +pct))."""
    rng = np.random.default_rng(seed)
    content = src.read_text()
    rows = _parse_block(content, "bus")
    new_lines = []
    for row in rows:
        cols = row.split()
        r = rng.uniform(-pct, pct)
        cols[2] = f"{float(cols[2]) * (1 + r):.6f}"  # Pd
        cols[3] = f"{float(cols[3]) * (1 + r):.6f}"  # Qd
        new_lines.append("\t" + "\t".join(cols) + ";")
    content = _replace_block(content, "bus", "\n".join(new_lines) + "\n")
    dst.write_text(content)
    print(f"Written: {dst.name}")


def make_wind_dispatch_m(
    src: Path,
    dst: Path,
    genfuel: list[str] | "pd.Index",
    pct: float = 0.10,
    seed: int = 7,
) -> None:
    """Curtail wind generators: scale Pg by (1 - Uniform(0, pct)).

    All ACTIVSg200 wind generators run at Pg == Pmax in the base case, so upward
    perturbation is not physically possible. Each wind gen is curtailed by a random
    fraction drawn from Uniform(0, pct), independently per generator.

    Parameters
    ----------
    genfuel : list[str] or pd.Index
        Per-generator fuel labels, as returned by ``read_matpower_case(src).genfuel``.
        Generator i is treated as wind when ``genfuel[i].lower() == 'wind'``.
    pct : float
        Maximum curtailment fraction (0..1). Each gen is reduced by Uniform(0, pct).
    """
    rng = np.random.default_rng(seed)
    content = src.read_text()
    gen_rows = _parse_block(content, "gen")
    fuels = [str(f).lower() for f in genfuel]
    if len(fuels) != len(gen_rows):
        raise ValueError(
            f"genfuel length ({len(fuels)}) != mpc.gen rows ({len(gen_rows)})"
        )
    new_gen_lines = []
    wind_count = 0
    for fuel, row in zip(fuels, gen_rows):
        cols = row.split()
        if fuel == "wind":
            r = rng.uniform(0.0, pct)  # curtailment fraction in [0, pct]
            pg_new = max(0.0, float(cols[1]) * (1.0 - r))
            cols[1] = f"{pg_new:.6f}"
            wind_count += 1
        new_gen_lines.append("\t" + "\t".join(cols) + ";")
    content = _replace_block(content, "gen", "\n".join(new_gen_lines) + "\n")
    dst.write_text(content)
    print(
        f"Written: {dst.name}  ({wind_count} wind gens curtailed, up to {pct*100:.0f}%)"
    )


def make_gen_off_m(
    src: Path,
    dst: Path,
    bus_num: int | None = None,
    n_random: int | None = None,
    seed: int = 0,
) -> frozenset[int]:
    """Take one or more generators offline: set GEN_STATUS=0, Pg=0, Qg=0.

    Two modes (mutually exclusive):
      bus_num=<int>          -- turn off all gens at that specific bus
      n_random=<int>, seed=  -- randomly choose n_random online non-slack gens
                                (status=1, Pg>0); reproducible via seed

    All three fields (status, Pg, Qg) are zeroed because GridKit's
    SystemSteadyStateModel does not filter on status; a gen with status=0
    but non-zero Pg would still inject power.
    """
    if (bus_num is None) == (n_random is None):
        raise ValueError("specify exactly one of bus_num or n_random")

    content = src.read_text()
    gen_rows = _parse_block(content, "gen")

    if bus_num is not None:
        off_buses = {bus_num}
        label = f"bus {bus_num}"
    else:
        rng = np.random.default_rng(seed)
        bus_rows = _parse_block(content, "bus")
        slack_buses = {int(r.split()[0]) for r in bus_rows if r.split()[1] == "3"}
        candidates = [
            int(r.split()[0])
            for r in gen_rows
            if int(r.split()[7]) == 1
            and float(r.split()[1]) > 0
            and int(r.split()[0]) not in slack_buses
        ]
        if n_random > len(candidates):
            raise ValueError(
                f"n_random={n_random} exceeds available online non-slack gens"
                f" ({len(candidates)})"
            )
        chosen = list(rng.choice(candidates, size=n_random, replace=False))
        off_buses = set(chosen)
        label = f"{n_random} random gens (seed={seed}): buses {sorted(chosen)}"

    new_gen_lines = []
    off_count = 0
    mw_dropped = 0.0
    for row in gen_rows:
        cols = row.split()
        if int(cols[0]) in off_buses:
            mw_dropped += float(cols[1])
            cols[7] = "0"  # GEN_STATUS
            cols[1] = "0.000000"  # Pg — must zero; GridKit ignores status flag
            cols[2] = "0.000000"  # Qg
            off_count += 1
        new_gen_lines.append("\t" + "\t".join(cols) + ";")
    content = _replace_block(content, "gen", "\n".join(new_gen_lines) + "\n")
    dst.write_text(content)
    print(
        f"Written: {dst.name}"
        f"  ({off_count} gen(s) offline: {label}  |  MW dropped: {mw_dropped:.1f})"
    )
    return frozenset(off_buses)
