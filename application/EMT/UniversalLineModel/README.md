# Universal Line Model

`UniversalLineModel` sweeps an overhead line description over frequency and
fits the time-domain line model coefficients: the characteristic admittance
$\mathbf{Y}_c(s)$ to `yc.model.json`, and the propagation function
$\mathbf{H}(s)$ to `propagation.model.json` under one of two treatments
selected by `--h-domain` (default `phase`). The artifact names its
treatment in a top-level `domain` field.

## Propagation Treatments

`--h-domain modal` extracts one minimum delay per mode and fits the
backwound factor pair:

```math
\mathbf{H}(s) = \mathbf{G}^{\mathrm{out}}(s)\,
  \mathrm{diag}\!\left(e^{-s\tau_1},\ldots,e^{-s\tau_M}\right)
  \mathbf{G}^{\mathrm{in}}(s)
```

`--h-domain phase` extracts the single minimum delay among all modes and
fits one phase-domain matrix:

```math
\mathbf{H}(s) = \hat{\mathbf{H}}(s)\,e^{-s\tau},
\qquad \tau = \min_m \tau_m
```

The phase-domain target $\hat{\mathbf{H}}$ is assembled per sample as
$\mathrm{conj}(\mathbf{T}_i)\,\mathrm{diag}(\mathbf{H})\,\mathbf{T}_v^{\mathsf{T}}$
before backwinding, so the eigenvector gauge cancels identically and the
fit carries no per-frequency phase convention. The modal treatment keeps
per-mode delays and is preferable when the modal delays separate into
distinct groups; the phase treatment fits one function and is preferable
when they cluster. The single delay carries the inter-mode spread as
residual winding proportional to frequency, so wider sweep bands shift
the balance toward the modal treatment.

## Files

| File | Contents |
| ---- | -------- |
| `UniversalLineModel.cpp` | Sweep, delay extraction, fits, artifacts |

Run `UniversalLineModel <line-json-file> --help` for the full option list.
