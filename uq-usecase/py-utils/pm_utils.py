"""Utilities for running and post-processing PowerModels.jl (pm_solve.jl) power flow cases.

pm_solve.jl's stdout format is identical to GridKit's solve_pf
(`bus <i>  V=<pu>  theta_deg=<deg>  type=<1|2|3>`), so the existing parser and
comparison helpers in pf_utils.py are reused here rather than duplicated.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pandas as pd

from pf_utils import _parse_pf_stdout  # noqa: F401  (re-exported for convenience)

# ---------------------------------------------------------------------------
# pm_solve.jl runners
# ---------------------------------------------------------------------------


def _run(cmd: list[str]) -> tuple[pd.DataFrame, str, int]:
    """Run a pm_solve.jl command, parse bus output from stdout."""
    r = subprocess.run(cmd, capture_output=True, text=True)
    return _parse_pf_stdout(r.stdout), r.stderr, r.returncode


def _tol_args(
    tol: float | None, max_iter: int | None, flat_start: bool = False
) -> list[str]:
    args: list[str] = []
    if tol is not None:
        args += ["--tol", str(tol)]
    if max_iter is not None:
        args += ["--max-iter", str(max_iter)]
    if flat_start:
        args += ["--flat-start"]
    return args


def run_pm_solve(
    julia_bin: Path,
    pm_project_dir: Path,
    pm_solve_jl: Path,
    m_path: Path,
    pf_type: str = "ac",
    solver: str = "ipopt",
    tol: float | None = None,
    max_iter: int | None = None,
    flat_start: bool = False,
) -> tuple[pd.DataFrame, str, int]:
    """Run pm_solve.jl on a .m file. Returns (bus_df, stderr_text, return_code)."""
    return _run(
        [
            str(julia_bin),
            f"--project={pm_project_dir}",
            str(pm_solve_jl),
            str(m_path),
            "--pf-type",
            pf_type,
            "--solver",
            solver,
        ]
        + _tol_args(tol, max_iter, flat_start)
    )


def run_pm_solve_out(
    julia_bin: Path,
    pm_project_dir: Path,
    pm_solve_jl: Path,
    m_path: Path,
    out_m_path: Path,
    pf_type: str = "ac",
    solver: str = "ipopt",
    tol: float | None = None,
    max_iter: int | None = None,
    flat_start: bool = False,
) -> tuple[pd.DataFrame, str, int]:
    """Run pm_solve.jl with --output-m. Returns (bus_df, stderr, return_code)."""
    return _run(
        [
            str(julia_bin),
            f"--project={pm_project_dir}",
            str(pm_solve_jl),
            str(m_path),
            "--pf-type",
            pf_type,
            "--solver",
            solver,
            "--output-m",
            str(out_m_path),
        ]
        + _tol_args(tol, max_iter, flat_start)
    )


def pm_summary(case_label: str, df: pd.DataFrame, rc: int, stderr: str) -> None:
    """Print convergence status, solver params, Ipopt iteration count, and V range."""
    converged = rc == 0
    print(f"{case_label}")
    # solver params line from [pm_solve]
    params_m = re.search(
        r"\[pm_solve\] solver=\S+\s+pf_type=\S+\s+tol=(\S+)\s+max_iter=(\S+)", stderr
    )
    if params_m:
        print(f"  tol={params_m.group(1)}  max_iter={params_m.group(2)}")
    # termination status
    for tag in ("termination_status", "solve_time", "converged"):
        m = re.search(rf"\[pm_solve\] {tag}=(\S+)", stderr)
        if m:
            print(f"  {tag}={m.group(1)}")
    # Ipopt iteration count from print_level=3 output redirected to stderr
    iters_m = re.search(r"Number of Iterations\.+:\s*(\d+)", stderr)
    nlp_m = re.search(r"Overall NLP error\.+:\s*(\S+)", stderr)
    if iters_m:
        iters_str = iters_m.group(1)
        nlp_str = nlp_m.group(1) if nlp_m else "n/a"
        print(f"  ipopt_iters={iters_str}  final_nlp_error={nlp_str}")
    print(
        f"  {'CONVERGED' if converged else 'NOT CONVERGED'}"
        f"   rc={rc}   buses={len(df)}"
    )
    if converged and not df.empty:
        viol = df[(df.V_pu < 0.95) | (df.V_pu > 1.05)]
        print(
            f"  V range: [{df.V_pu.min():.4f}, {df.V_pu.max():.4f}] pu"
            f"   violations (V<0.95 or V>1.05): {len(viol)}"
        )
    print()
