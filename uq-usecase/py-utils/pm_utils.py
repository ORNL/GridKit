"""Utilities for running and post-processing PowerModels.jl (pm_solve.jl) power flow cases.

pm_solve.jl's stdout format is identical to GridKit's solve_pf
(`bus <i>  V=<pu>  theta_deg=<deg>  type=<1|2|3>`), so the existing parser and
comparison helpers in pf_utils.py are reused here rather than duplicated.
"""

from __future__ import annotations

import json
import re
import subprocess
from functools import lru_cache
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


# ---------------------------------------------------------------------------
# backend-agnostic entry point
# ---------------------------------------------------------------------------
#
# The contract, which both backends satisfy by construction (pm_solve.jl was
# written to mirror solve_pf's CLI and stdout):
#
#   in   : a .m with injections set (PD/QD, PG, GEN_STATUS)
#   out  : a .m with solved VM/VA, and PG/QG updated for the slack and
#          PV-bus reactive output
#   also : one line per bus on stdout, "bus <i> V=<pu> theta_deg=<deg> type=<t>"
#   rc   : 0 on convergence
#
# Swapping backends should be a one-constant change. Run
# `compare_backends()` first: a backend that converges but disagrees on the
# operating point is worse than one that fails loudly.

_GRIDKIT = Path.home() / "gridkit"

# Existing-but-broken Julia builds are the worst failure mode here: the
# Sienna-Ops 1.12.5 install hangs forever before executing a single line, so
# `exe.exists()` is True and every downstream call blocks with no output.
# Candidates are therefore probed by *running* them, in preference order.
JULIA_CANDIDATES: tuple[Path, ...] = (
    Path("/nopt/nrel/apps/cpu_stack/software/julia/julia-1.12.1/usr/bin/julia"),
    Path.home() / "bin/julia110",
    Path.home() / "bin/julia",
    Path.home() / "bin/julia112",
)

_JULIA_PROBE_S = 60.0


def julia_probe(exe: str | Path, timeout: float = _JULIA_PROBE_S) -> tuple[bool, str]:
    """Run `println(VERSION)` and report whether this binary actually works."""
    exe = Path(exe)
    if not exe.exists():
        return False, "not found"
    try:
        r = subprocess.run(
            [str(exe), "--startup-file=no", "-e", "println(VERSION)"],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return False, f"hangs (no output in {timeout:g}s) -- broken install"
    if r.returncode != 0:
        return False, f"exited {r.returncode}: {r.stderr.strip()[:120]}"
    return True, r.stdout.strip()


@lru_cache(maxsize=1)
def resolve_julia() -> Path | None:
    """First Julia in JULIA_CANDIDATES that actually runs. None if all fail."""
    for exe in JULIA_CANDIDATES:
        ok, _ = julia_probe(exe)
        if ok:
            return exe
    return None


PF_BACKENDS: dict[str, dict] = {
    "pm": {
        "exe": None,  # resolved by resolve_julia(); see JULIA_CANDIDATES
        "resolve": "julia",
        "project": _GRIDKIT / "uq-usecase/pm-solver",
        "script": _GRIDKIT / "uq-usecase/pm-solver/pm_solve.jl",
        "note": "PowerModels.jl, Ipopt. Validated; current default.",
    },
    "gridkit": {
        "exe": _GRIDKIT / "uq-usecase/pf-solver/build/solve_pf",
        "resolve": None,
        "project": None,
        "script": None,
        "note": (
            "GridKit steady-state solver via uq-usecase/pf-solver. "
            "No PV->PQ switching. Needs a rebuild against GridKit v0.2.0: "
            "the Vector API lost .size()/operator[]."
        ),
    },
}


def backend_exe(name: str) -> Path | None:
    """Resolved executable for a backend, or None if none works."""
    cfg = PF_BACKENDS[name]
    if cfg.get("resolve") == "julia":
        return resolve_julia()
    return Path(cfg["exe"])


def pf_backend_status(name: str) -> dict:
    """Whether a backend is present and runnable, without running a case."""
    cfg = PF_BACKENDS[name]
    exe = backend_exe(name)
    note = cfg["note"]
    if cfg.get("resolve") == "julia":
        if exe is None:
            tried = "\n".join(
                f"      {p}: {julia_probe(p)[1]}" for p in JULIA_CANDIDATES
            )
            note = f"{note}\n    no working Julia found. Tried:\n{tried}"
        else:
            note = f"{note} Using {exe} ({julia_probe(exe)[1]})."
    out = {
        "backend": name,
        "exe": str(exe) if exe is not None else None,
        "exists": exe is not None and exe.exists(),
        "note": note,
    }
    if cfg["script"] is not None:
        out["script_exists"] = Path(cfg["script"]).exists()
    return out


def solve_pf(
    m_path: str | Path,
    out_m_path: str | Path,
    backend: str = "pm",
    pf_type: str = "ac",
    solver: str = "ipopt",
    tol: float | None = None,
    max_iter: int | None = None,
    flat_start: bool = False,
) -> tuple[pd.DataFrame, str, int]:
    """Solve a power flow with the chosen backend. Returns (bus_df, stderr, rc).

    `backend` is a key of `PF_BACKENDS`. Keep this the only place a solver
    is named, so switching is a one-line change.
    """
    if backend not in PF_BACKENDS:
        raise KeyError(f"unknown backend {backend!r}; have {sorted(PF_BACKENDS)}")
    cfg = PF_BACKENDS[backend]
    exe = backend_exe(backend)
    if exe is None or not exe.exists():
        raise FileNotFoundError(
            f"backend {backend!r} has no working executable\n"
            f"{pf_backend_status(backend)['note']}"
        )

    if backend == "pm":
        return run_pm_solve_out(
            exe,
            cfg["project"],
            cfg["script"],
            m_path,
            out_m_path,
            pf_type=pf_type,
            solver=solver,
            tol=tol,
            max_iter=max_iter,
            flat_start=flat_start,
        )

    cmd = [str(exe), str(m_path), "--output-m", str(out_m_path)]
    if flat_start:
        cmd.append("--flat-start")
    return _run(cmd)


def compare_backends(
    m_path: str | Path,
    out_dir: str | Path,
    backends: tuple[str, ...] = ("pm", "gridkit"),
    tol_v: float = 1e-4,
    tol_deg: float = 1e-2,
) -> pd.DataFrame:
    """Solve the same case with each backend and compare bus voltages.

    Conformance test for swapping solvers. Agreement on |V| and angle is the
    minimum bar; a backend that converges to a different operating point will
    silently change every downstream QoI.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    res, rows = {}, []
    for b in backends:
        st = pf_backend_status(b)
        if not st["exists"]:
            rows.append({"backend": b, "status": "missing", "note": st["note"]})
            continue
        try:
            df, err, rc = solve_pf(m_path, out_dir / f"solved_{b}.m", backend=b)
        except Exception as exc:  # a crashing backend must not abort the compare
            rows.append({"backend": b, "status": f"error: {type(exc).__name__}"})
            continue
        (out_dir / f"{b}.log").write_text(err)
        rows.append(
            {"backend": b, "status": "ok" if rc == 0 else f"rc={rc}", "n_bus": len(df)}
        )
        if rc == 0 and len(df):
            res[b] = df.set_index("bus_i")

    summary = pd.DataFrame(rows)
    names = list(res)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            j_ = res[a].join(res[b], lsuffix=f"_{a}", rsuffix=f"_{b}", how="inner")
            dv = (j_[f"V_pu_{a}"] - j_[f"V_pu_{b}"]).abs()
            da = (j_[f"theta_deg_{a}"] - j_[f"theta_deg_{b}"]).abs()
            summary.attrs[f"{a}_vs_{b}"] = {
                "max_dV_pu": float(dv.max()),
                "max_dtheta_deg": float(da.max()),
                "agree": bool(dv.max() <= tol_v and da.max() <= tol_deg),
            }
    return summary


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


def run_pm_solve_batch(
    julia_bin: Path,
    pm_project_dir: Path,
    pm_solve_jl: Path,
    cases: list[dict],
    workdir: Path,
    solver: str = "ipopt",
    tol: float | None = None,
    max_iter: int | None = None,
    flat_start: bool = False,
) -> tuple[list[dict], str, int]:
    """Solve many .m files in ONE Julia process.

    Julia startup plus package load costs ~10 s; a warm solve costs 0.03-3 s.
    Batching therefore turns an N-case sweep from N*(10+solve) into
    10+N*solve.

    ``cases`` is a list of ``{"input": path, "output": path, "pf_type": str}``.

    Returns ``(records, stderr, rc)`` where ``records`` has one entry per case,
    in manifest order, each carrying an ``ok`` flag. A case that fails to
    converge does not abort the batch and does not get an output .m written,
    so a caller must check ``ok`` before advancing that case to a dynamic
    simulation. ``rc != 0`` means the batch itself could not run.
    """
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    manifest = workdir / "pm_batch_manifest.json"
    results = workdir / "pm_batch_results.json"
    payload = [
        {
            "input": str(c["input"]),
            "output": str(c.get("output", "")),
            "pf_type": c.get("pf_type", "ac"),
        }
        for c in cases
    ]
    manifest.write_text(json.dumps(payload, indent=1))
    results.unlink(missing_ok=True)

    cmd = [
        str(julia_bin),
        f"--project={pm_project_dir}",
        str(pm_solve_jl),
        "--batch",
        str(manifest),
        "--results",
        str(results),
        "--solver",
        solver,
    ] + _tol_args(tol, max_iter, flat_start)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not results.exists():
        return [], r.stderr, r.returncode or 1
    return json.loads(results.read_text()), r.stderr, 0


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
