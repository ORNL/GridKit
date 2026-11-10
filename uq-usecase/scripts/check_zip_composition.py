"""Does GridKit compose two LoadZIP devices on one bus correctly?

Three variants of the SAME total load at one bus, same dispatch, same fault:

  A  one device,  Pnom=P,     alphaZ=1          <- how the case ships
  B  two devices, Pnom=P/2 each, alphaZ=1       <- must be IDENTICAL to A
  C  two devices, host alphaZ=1 + second alphaP=1
  D  one device,  Pnom=P,     alphaP=1          <- known to work on existing loads

A vs B isolates "can two LoadZIPs share a bus at all" from any ZIP question.
C vs D isolates "is a second device different from changing the first".
"""

import copy
import sys

import pandas as pd

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
import campaign as cp  # noqa: E402

DATA = "/kfs2/projects/scidac/scidac-data/ACTIVSg200/grid-v1"
HOUR = 6164
BUS = 129  # 87.1 MW, the largest single load; runs at alphaP=1 as one device


def qoi(csv):
    df = pd.read_csv(csv)
    df.columns = [c.strip() for c in df.columns]
    df = df.drop_duplicates(subset="t", keep="last")
    vr = df[[c for c in df.columns if c.endswith("_Vr")]].to_numpy()
    vi = df[[c for c in df.columns if c.endswith("_Vi")]].to_numpy()
    w = df[[c for c in df.columns if c.endswith("_omega")]].to_numpy()
    return (vr**2 + vi**2) ** 0.5, w


def main():
    cfg = cp.CampaignConfig(name="ziptest", scratch="/tmp/ziptest")
    dem = pd.read_csv(f"{DATA}/demand.csv")
    wnd = pd.read_csv(f"{DATA}/wind.csv")
    demc = [c for c in dem.columns if c != "DateTime"]
    wndc = [c for c in wnd.columns if c != "DateTime"]
    wg = {c: cfg.gid[int(c.replace("gen", "")) - 1] for c in wndc}
    d = {int(x.replace("bus", "")): float(dem.loc[HOUR, x]) for x in demc}
    w = {wg[x]: float(wnd.loc[HOUR, x]) for x in wndc}
    net = sum(d.values()) - sum(w.values())
    lo, hi = cp.commitment_band(cfg, net, cp.merit_order(cfg))
    n = int(round(lo + 0.5 * (hi - lo)))

    ((case, rec),) = cp.solve_points_batch(
        cfg, [{"demand": d, "wind": w, "n_commit": n, "workdir": cfg.scratch / "p"}]
    )
    host = next(
        x
        for x in case["devices"]
        if x["class"] == "LoadZIP"
        and x["ports"]["bus"] == BUS
        and not x["id"].endswith("_ibr_stand_in")
    )
    P, Q = host["params"]["Pnom"], host["params"]["Qnom"]
    print(
        f"dispatch ok. bus {BUS} host load {P*cp.BASEMVA:.1f} MW "
        f"({len([x for x in case['devices'] if x['class']=='LoadZIP' and x['ports']['bus']==BUS])} device)"
    )

    def build(variant):
        c = copy.deepcopy(case)
        h = next(x for x in c["devices"] if x["id"] == host["id"])
        if variant == "A":
            pass
        elif variant == "D":
            h["params"]["alphaP"] = 1.0
            h["params"]["alphaI"] = 0.0
        else:
            h["params"]["Pnom"], h["params"]["Qnom"] = P / 2, Q / 2
            second = {
                "class": "LoadZIP",
                "ports": {"bus": BUS},
                "id": f"second_{BUS}",
                "params": {
                    "Pnom": P / 2,
                    "Qnom": Q / 2,
                    "alphaI": 0.0,
                    "alphaP": 1.0 if variant == "C" else 0.0,
                },
            }
            c["devices"].append(second)
        return c

    dist = cp.Disturbance(fault_bus=66, duration=0.1)
    out = {}
    for v in ("A", "B", "C", "D"):
        wd = cfg.scratch / v
        r = cp.run_disturbance(cfg, build(v), dist, wd, with_cct=False)
        csv = next(iter(wd.glob("*.csv")), None)
        out[v] = (r["status"], csv)
        ndev = len(
            [
                x
                for x in build(v)["devices"]
                if x["class"] == "LoadZIP" and x["ports"]["bus"] == BUS
            ]
        )
        print(f"  {v}: {ndev} device(s) at bus {BUS}  -> {r['status']}")

    if out["A"][1] and out["B"][1]:
        va, wa = qoi(out["A"][1])
        vb, wb = qoi(out["B"][1])
        if va.shape == vb.shape:
            dv = abs(va - vb).max()
            dw = abs(wa - wb).max()
            print(
                f"\nA vs B (must be identical): max|dV| = {dv:.3e} pu, "
                f"max|domega| = {dw:.3e} pu"
            )
            # Integrator noise floor is ~1e-8 pu (a flat run measures
            # max|omega| = 7.2e-09). Splitting one device into two changes the
            # residual structure, so IDA takes slightly different steps.
            # Anything at or below ~1e-6 pu is the same trajectory.
            tol = 1e-6
            print(f"  (tolerance {tol:.0e} pu; integrator noise floor ~1e-8)")
            print(
                "  VERDICT:",
                (
                    "two LoadZIPs compose correctly"
                    if dv < tol and dw < tol
                    else "*** MISMATCH - composition bug ***"
                ),
            )
        else:
            print(f"\nA vs B shapes differ: {va.shape} vs {vb.shape}")


if __name__ == "__main__":
    main()
