# Universal Line Model

`UniversalLineModel` sweeps an overhead line description over frequency and
fits the time-domain line model coefficients: the characteristic admittance
$\mathbf{Y}_c(s)$ to `yc.model.json`, and the propagation function
$\mathbf{H}(s)$ to `propagation.model.json`.

## Propagation

The propagation function is one rational matrix behind one shared delay,
the smallest transport delay the line carries:

```math
\mathbf{H}(s) = \mathbf{H}^{\min}(s)\,e^{-s\tau},
\qquad \tau = \min_{m,\omega}\tau_m(\omega)
```

The fitting target $\mathbf{H}^{\min}$ is the modal propagation unwound by
that delay and carried back to phase coordinates:

```math
\mathbf{H}^{\min}(s) = \mathrm{conj}(\mathbf{T}_i)\,
  \mathrm{diag}\!\left(h_1(s)e^{+s\tau},\ldots,h_M(s)e^{+s\tau}\right)
  \mathbf{T}_v^{\mathsf{T}}
```

The eigenvector normalization cancels identically in that product, so the
target carries no per-frequency phase convention and no structural zeros,
and the inter-mode cancellation survives into the fit. What the single
delay leaves behind is the spread of the modal delays, which the fit must
absorb as a residual winding growing with frequency; the application
reports that spread, and a line whose modal delays separate widely needs
either more poles or a narrower band.

## Files

| File | Contents |
| ---- | -------- |
| `UniversalLineModel.cpp` | Sweep, delay extraction, fits, artifacts |

Run `UniversalLineModel <line-json-file> --help` for the full option list.
