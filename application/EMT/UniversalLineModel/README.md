# Universal Line Model

`UniversalLineModel` sweeps an overhead line description over frequency and
fits the time-domain line model coefficients: the characteristic admittance
$\mathbf{Y}_c(s)$ to `yc.model.json`, and the propagation function
$\mathbf{H}(s)$ to `propagation.model.json`.

## Propagation

The propagation function is a sum of rational matrices, one per mode, each
behind its own transport delay:

```math
\mathbf{H}(s) = \sum_m \mathbf{H}^{\min}_m(s)\,e^{-s\tau_m}
```

Each fitting target is the mode's rank-one dyad unwound by its own delay,

```math
\mathbf{H}^{\min}_m(s) = \mathrm{conj}(\mathbf{t}_{i,m})\,
  h_m(s)\,e^{+s\tau_m}\,\mathbf{t}_{v,m}^{\mathsf{T}},
```

where $\tau_m$ is the smallest sampled delay over the frequencies at which
$|h_m|$ still stands above the fitter's magnitude floor, so samples the fit
treats as zeros never decide a delay. The eigenvector normalization cancels
inside each dyad, so the targets carry no per-frequency phase convention.

Run `UniversalLineModel <line-json-file> --help` for the option list.
