# Rational Approximation

This directory holds the model space for rational approximation of sampled
frequency responses: the sampled-data carrier, the rational model itself, and
the operations performed on a model or its samples — delay preprocessing,
passivity assessment, and (planned) exact-rank state-space realization.
Estimators live
beside this directory; [vector fitting](../VectorFitting/README.md) is the
first. In short: everything one does with a rational model except fit it.

## Rational Model

`RationalModel` stores the coefficients of

```math
\hat{\mathbf{H}}(s) = \mathbf{D} + s\mathbf{E}
  + \sum_{q=1}^{Q} \dfrac{\mathbf{R}_q}{s - p_q}
```

with output dimension $N$, input dimension $K$, and pole count $Q$. The
coefficient conventions, units, and conjugate-pair ordering are those of the
[`VectorFit` model](../../../Model/EMT/Operators/Rational/VectorFit/README.md):
real poles carry real residues, and every nonreal pole and residue is
immediately followed by its conjugate. Estimators in this family produce
models satisfying that contract by construction, so the coefficients feed the
EMT rational operators without translation. Stability is always computed from
the poles, never cached.

## Sampled Response

`SampledResponse` carries $\mathbf{H}_m = \mathbf{H}(\mathrm{j}\omega_m) \in
\mathbb{C}^{N \times K}$ on strictly increasing, positive angular frequencies
$\omega_1 < \dots < \omega_M$, stored sample-major and row-major within each
sample, with an element accessor so callers never index the flat storage
directly.

## Delay Preprocessing

A propagation function with transport delay is fitted after minimum-phase
shifting. With the delay $\tau = \min_\omega \tau(\omega)$ taken over the
sweep, the shifted response

```math
\mathbf{H}^{\mathrm{mps}}(\mathrm{j}\omega)
  = \mathbf{H}(\mathrm{j}\omega)\,e^{\mathrm{j}\omega\tau}
```

is smooth enough for a low-order strictly proper fit [2]. The extracted
delays accompany the fitted model for the `Delay` and `Propagation`
operators.

## State-Space Realization

`StateSpaceRealization.hpp` fixes the carrier for the factorized form
$\hat{\mathbf{H}}(s) = \mathbf{D} + s\mathbf{E}
+ \mathbf{C}(s\mathbf{I} - \mathbf{P})^{-1}\mathbf{B}$, matching the
[`StateSpace` model](../../../Model/EMT/Operators/Rational/StateSpace/README.md)
coefficient contract (`D`, `E`, `poles`, `C`, `B`). The factorization of
a rational model into this carrier -- each residue matrix at its
numerical rank, so the realization is exact and reduces to the rank-one
form for scalar and vector responses -- is planned and not implemented
yet; nothing exports it today.

## Passivity

A fitted admittance intended for time-domain use is first gated on
stability, then screened by sweeping the smallest eigenvalue of the
Hermitian part of $\mathbf{Y}(\mathrm{j}\omega)$ at DC, over a
logarithmic grid spanning the union of the fitted band and the pole
magnitudes, and at the asymptotic limit through the constant term.
Violation band edges inside the grid are refined by bisection and
reported as frequency bands; a DC violation starts at zero exactly, and
an indefinite constant term extends its band to infinity. The exact
Hamiltonian eigenvalue certification [1] and enforcement by minimal
coefficient perturbation are planned follow-ons; the latter reuses the
constrained coefficient identification of the estimator.

## Files

| File | Contents |
| ---- | -------- |
| `SampledResponse.hpp` | Sampled-response carrier |
| `RationalModel.hpp` | Rational model, evaluation, stability check |
| `MinimumPhase.hpp` | Delay extraction and minimum-phase shifting |
| `StateSpaceRealization.hpp` | Exact-rank factorized realization |
| `Passivity.hpp` | Passivity report and assessment |

## References

1. B. Gustavsen and A. Semlyen, "Enforcing passivity for admittance matrices
   approximated by rational functions," IEEE Transactions on Power Systems,
   vol. 16, no. 1, pp. 97-104, 2001.
2. A. Morched, B. Gustavsen, and M. Tartibi, "A universal model for accurate
   calculation of electromagnetic transients on overhead lines and underground
   cables," IEEE Transactions on Power Delivery, vol. 14, no. 3,
   pp. 1032-1038, 1999.
