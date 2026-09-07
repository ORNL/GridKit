"""Parameter sampling for UQ workflows."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import norm as sp_norm, qmc


def generate_samples(param_specs, N, seed=None, method="lhs"):
    """
    Sample parameters according to per-spec distributions.

    Parameters
    ----------
    param_specs : list of dict, each with keys:
        "class"  - device class string (e.g. "Genrou")
        "id"     - device id string    (e.g. "genrou_2_1")
        "param"  - parameter name      (e.g. "H")
        "dist"   - one of "uniform" (default) or "normal"

        For dist="uniform":
            Option A — percent:  "nominal" + "pct"
                lo = nominal * (1 - pct),  hi = nominal * (1 + pct)
            Option B — fixed:    "lo" + "hi"

        For dist="normal":
            "mean" and "std"

    N      : int, number of samples
    seed   : int or None, for reproducibility
    method : "lhs"    — Latin Hypercube (better space-filling, recommended)
             "random" — independent random draws per parameter

    Returns
    -------
    pd.DataFrame shape (N, len(param_specs))
        columns named "{class}_{id}_{param}", e.g. "Genrou_genrou_2_1_H"

    Notes
    -----
    LHS works by drawing N uniform samples stratified across [0,1] per
    dimension, then mapping each through the inverse CDF (ppf) of the
    requested distribution. This gives correct marginals for both uniform
    and normal distributions while maximising coverage of the input space.
    """
    cols = [f"{s['class']}_{s['id']}_{s['param']}" for s in param_specs]
    d = len(param_specs)

    if method == "lhs":
        sampler = qmc.LatinHypercube(d=d, seed=seed)
        unit = sampler.random(N)  # N × d uniform in (0, 1)
    else:
        rng = np.random.default_rng(seed)
        unit = rng.uniform(0, 1, size=(N, d))

    data = {}
    for j, spec in enumerate(param_specs):
        dist = spec.get("dist", "uniform")
        u = unit[:, j]

        if dist == "uniform":
            if "lo" in spec and "hi" in spec:
                lo, hi = spec["lo"], spec["hi"]
            else:
                lo = spec["nominal"] * (1 - spec["pct"])
                hi = spec["nominal"] * (1 + spec["pct"])
            # uniform ppf: lo + u*(hi-lo)
            data[cols[j]] = lo + u * (hi - lo)

        elif dist == "normal":
            mean, std = spec["mean"], spec["std"]
            # normal ppf via scipy
            data[cols[j]] = sp_norm.ppf(u, loc=mean, scale=std)

        else:
            raise ValueError(f"Unknown dist '{dist}' in param_spec for {cols[j]}")

    return pd.DataFrame(data)
