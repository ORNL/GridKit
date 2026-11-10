"""Outer-loop PV->PQ enforcement of generator reactive limits for the AC PF.

PowerModels' AC PF (solve_ac_pf) holds voltage at every PV bus and lets Q float.
This module adds the usual planning-tool loop (MATPOWER enforce_q_lims style):

  1. Solve the PF with the current set of fixed units.
  2. An active PV unit whose Q leaves [qmin, qmax] is fixed at the violated limit.
     Its P and Q move into the bus load (Pd -= Pg, Qd -= Q_fixed), its gen status
     goes to 0, and its bus becomes PQ if no other active unit is attached.
  3. A fixed unit returns to PV when the bus voltage moves back past its setpoint
     (V > Vg for a max-limited unit, V < Vg for a min-limited unit).
  4. Repeat until no unit changes state. A unit that has returned to PV twice is
     kept fixed (anti-cycling).

Slack units are never switched.
"""

import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
import pm_utils  # noqa: E402

PROJ = Path("/home/isatkaus/gridkit/uq-usecase/pm-solver")
PM_SOLVE = PROJ / "pm_solve.jl"
TOL_Q = 0.1     # MVAr
TOL_V = 1e-4    # pu
MAX_ITER = 30
MAX_RETURNS = 2


@dataclass
class Unit:
    row: int
    bus: int
    pg: float
    qmax: float
    qmin: float
    vg: float
    slack: bool


def _rows(txt, name):
    i = txt.index(name)
    j = txt.index("[", i) + 1
    k = txt.index("];", j)
    out = []
    for line in txt[j:k].splitlines():
        body = line.split("%")[0].strip().rstrip(";").strip()
        if body:
            out.append(body.split())
    return out


def _replace_block(txt, name, rows):
    i = txt.index(name)
    j = txt.index("[", i) + 1
    k = txt.index("];", j)
    body = "\n" + "\n".join("\t".join(r) + ";" for r in rows) + "\n"
    return txt[:j] + body + txt[k:]


def _units(gen_rows, bus_rows):
    slack_buses = {int(r[0]) for r in bus_rows if r[1] == "3"}
    units = {}
    for idx, r in enumerate(gen_rows):
        if r[7] == "0":
            continue
        bus = int(r[0])
        units[idx] = Unit(row=idx, bus=bus, pg=float(r[1]), qmax=float(r[3]),
                          qmin=float(r[4]), vg=float(r[5]), slack=bus in slack_buses)
    return units


def write_case(base_txt, bus_rows, gen_rows, units, state, vm, out_path):
    """Write a .m with fixed units moved into the bus load and gen status set to 0."""
    new_bus = [list(r) for r in bus_rows]
    new_gen = [list(r) for r in gen_rows]
    bus_index = {int(r[0]): n for n, r in enumerate(new_bus)}
    pd_shift = {}
    qd_shift = {}
    for i, (kind, qfix) in state.items():
        u = units[i]
        pd_shift[u.bus] = pd_shift.get(u.bus, 0.0) + u.pg
        qd_shift[u.bus] = qd_shift.get(u.bus, 0.0) + qfix
        new_gen[i][2] = str(qfix)  # status-0 rows keep their Qg; store the limit so it is not stale
        new_gen[i][7] = "0"
    active_at = {}
    for i, u in units.items():
        if i not in state:
            active_at.setdefault(u.bus, []).append(i)
    for bus, n in bus_index.items():
        r = new_bus[n]
        pd0, qd0 = float(bus_rows[n][2]), float(bus_rows[n][3])
        r[2] = f"{pd0 - pd_shift.get(bus, 0.0):.10g}"
        r[3] = f"{qd0 - qd_shift.get(bus, 0.0):.10g}"
        if vm and bus in vm:
            r[7] = f"{vm[bus]:.10g}"
        if r[1] == "3":
            continue
        r[1] = "2" if active_at.get(bus) else "1"
    txt = _replace_block(base_txt, "mpc.bus =", new_bus)
    txt = _replace_block(txt, "mpc.gen =", new_gen)
    out_path.write_text(txt)


def read_solved(path):
    txt = path.read_text()
    bus = {int(r[0]): float(r[7]) for r in _rows(txt, "mpc.bus =")}
    gen = [r for r in _rows(txt, "mpc.gen =")]
    return bus, [float(r[2]) for r in gen], [int(r[0]) for r in gen]


def run(base_m, out_dir, log=print):
    out_dir.mkdir(parents=True, exist_ok=True)
    base_txt = base_m.read_text()
    bus_rows = _rows(base_txt, "mpc.bus =")
    gen_rows = _rows(base_txt, "mpc.gen =")
    units = _units(gen_rows, bus_rows)
    jb = pm_utils.resolve_julia()
    state, returns, vm, history = {}, {}, None, []
    converged = False
    for it in range(MAX_ITER):
        m_in = out_dir / f"it{it:02d}.m"
        m_out = out_dir / f"it{it:02d}.solved.m"
        write_case(base_txt, bus_rows, gen_rows, units, state, vm, m_in)
        _, err, rc = pm_utils.run_pm_solve_out(jb, PROJ, PM_SOLVE, m_in, m_out)
        if rc != 0 or not m_out.exists():
            raise RuntimeError(f"PF failed at iteration {it}, rc={rc}\n{err[-1500:]}")
        vm, qg, gbus = read_solved(m_out)
        changes = []
        for i, u in units.items():
            if u.slack:
                continue
            if i in state:
                kind = state[i][0]
                if returns.get(i, 0) >= MAX_RETURNS:
                    continue
                if (kind == "max" and vm[u.bus] > u.vg + TOL_V) or \
                   (kind == "min" and vm[u.bus] < u.vg - TOL_V):
                    del state[i]
                    returns[i] = returns.get(i, 0) + 1
                    changes.append(f"return unit{i} bus{u.bus} V={vm[u.bus]:.5f} Vg={u.vg:.5f}")
            else:
                q = qg[i]
                if q > u.qmax + TOL_Q:
                    state[i] = ("max", u.qmax)
                    changes.append(f"fix max unit{i} bus{u.bus} Q={q:.2f} qmax={u.qmax:.2f}")
                elif q < u.qmin - TOL_Q:
                    state[i] = ("min", u.qmin)
                    changes.append(f"fix min unit{i} bus{u.bus} Q={q:.2f} qmin={u.qmin:.2f}")
        history.append({"iter": it, "n_fixed": len(state), "changes": changes})
        log(f"iter {it}: fixed={len(state)} changes={len(changes)}")
        for c in changes:
            log(f"    {c}")
        if not changes:
            converged = True
            break
    return {"converged": converged, "state": state, "units": units,
            "history": history, "final_solved": m_out, "vm": vm, "qg": qg}


def main():
    base = Path("/scratch/isatkaus/gridkit-campaigns/ida_zoom/mw_balanced/pf/warm/mw200.m")
    out = Path("/scratch/isatkaus/gridkit-campaigns/ida_zoom/q_limit/mw200")
    res = run(base, out)
    print("converged:", res["converged"], "fixed units:", len(res["state"]))
    for i, (kind, q) in sorted(res["state"].items()):
        u = res["units"][i]
        print(f"  unit{i} bus{u.bus} {kind} Qfix={q:.2f} MVAr")


if __name__ == "__main__":
    main()
