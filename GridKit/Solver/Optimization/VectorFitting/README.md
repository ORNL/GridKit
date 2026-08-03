# Vector Fitting

`VectorFitting` computes a rational approximation of a sampled frequency
response. Its two linear stages are convex least squares problems solved in
closed form by orthogonal block elimination, with one small eigenvalue
computation per pole-relocation pass; only the genuinely nonlinear
refinement stage is posed to the `Optimization` family's Ipopt backend. The
estimator consumes and produces the
[rational-approximation model space](../Rational/README.md): a
`SampledResponse` in, a `RationalModel` out.

## Problem Statement

Given frequency samples $\omega_1 < \omega_2 < \dots < \omega_M$ with
responses $\mathbf{H}_m = \mathbf{H}(\mathrm{j}\omega_m) \in
\mathbb{C}^{N \times K}$ and weights $w_m > 0$, find the lowest pole count
$Q$ and coefficients such that

```math
\hat{\mathbf{H}}(s) = \mathbf{D} + s\mathbf{E}
  + \sum_{q=1}^{Q} \dfrac{\mathbf{R}_q}{s - p_q}
```

solves

```math
\min \; \sum_{m=1}^{M} w_m
  \left\lVert \mathbf{H}_m - \hat{\mathbf{H}}(\mathrm{j}\omega_m)
  \right\rVert_F^2
\quad \text{subject to} \quad
\mathrm{Re}(p_q) \le -\sigma_{\min}
```

within a target relative error. All matrix elements share the poles. The
fitted `RationalModel` satisfies the conjugate-pair contract of the model
space by construction.

## Real Parametrization

The optimization variables are real. A pole is either real, $p_q = a_q$, or
one of a conjugate pair stored as $(a_q, \omega_q)$ with
$p_q = a_q + \mathrm{j}\omega_q$. Residues follow the same split:
$\mathbf{R}_q = \mathbf{A}_q$ for a real pole and
$\mathbf{R}_q = \mathbf{A}_q + \mathrm{j}\mathbf{B}_q$ for a pair. Each real
pole contributes the basis function

```math
\varphi_q(s) = \dfrac{1}{s - p_q}
```

and each conjugate pair contributes the two real-coefficient functions

```math
\varphi_q(s) = \dfrac{1}{s - p_q} + \dfrac{1}{s - p_q^{\ast}},
\qquad
\psi_q(s) = \dfrac{\mathrm{j}}{s - p_q} - \dfrac{\mathrm{j}}{s - p_q^{\ast}}.
```

Real and imaginary parts of every sample equation are stacked as separate
real rows. Conjugate symmetry of the fitted model is therefore exact by
construction and never repaired after the fact.

## Algorithm

### Stage A: Pole Relocation

Relaxed vector fitting [2] with the current poles $p_1, \dots, p_Q$
introduces the scaling function

```math
\sigma(s) = \tilde{d} + \sum_{q=1}^{Q} \dfrac{\tilde{c}_q}{s - p_q}
```

and requires $\sigma(s)\mathbf{H}(s)$ to be rational with the same poles. At
the samples this relation is linear in the unknowns
$(\tilde{\mathbf{c}}, \tilde{d})$ and the surrogate coefficients
$(\tilde{\mathbf{D}}, \tilde{\mathbf{E}}, \tilde{\mathbf{R}}_q)$. One
relocation pass is the equality-constrained convex quadratic program

```math
\min \; \sum_{m=1}^{M} w_m
  \left\lVert \sigma(\mathrm{j}\omega_m)\mathbf{H}_m
  - \tilde{\mathbf{D}} - \mathrm{j}\omega_m\tilde{\mathbf{E}}
  - \sum_{q=1}^{Q}
    \dfrac{\tilde{\mathbf{R}}_q}{\mathrm{j}\omega_m - p_q}
  \right\rVert_F^2
\quad \text{subject to} \quad
\sum_{m=1}^{M} \mathrm{Re}\,\sigma(\mathrm{j}\omega_m) = M.
```

Classical implementations append the relaxation condition as a weighted
extra row of an augmented least-squares system; here it is an honest linear
equality, honored exactly. The coupled system is block arrow — each matrix
element's rows touch only that element's surrogate coefficients plus the
shared $(\tilde{\mathbf{c}}, \tilde{d})$ — so the surrogate blocks are
eliminated by one unpivoted thin QR per element, the fast formulation of
[3]. The surviving triangles stack into a small least-squares problem over
$(\tilde{\mathbf{c}}, \tilde{d})$ alone, and the relaxation equality is
absorbed by a Householder change of variables. Cost grows linearly with the
element count, no normal equations are formed, and the pass yields the
exact minimizer of the surrogate program.

The relocated poles are the zeros of $\sigma$, obtained as the eigenvalues
of the companion update assembled in its real block form (a 2-by-2 rotation
block per conjugate pair), so eigenvalues emerge in exact conjugate pairs.
Poles in the right half plane are reflected during iteration. Passes repeat
until the largest relative pole displacement falls below the tolerance.

### Stage B: Coefficient Identification

With poles frozen, identifying $\mathbf{D}$, $\mathbf{E}$, and the residues
is linear: every matrix element regresses on the same weighted basis, so
one rank-revealing factorization serves all elements with a triangular
solve each.

### Stage C: Refinement (optional)

Starting from the Stage A/B solution, the full nonlinear program over pole
parameters and coefficients minimizes the true weighted misfit subject to
$\mathrm{Re}(p_q) \le -\sigma_{\min}$ as inequality constraints. This
replaces the reflection heuristic as the final word on stability and yields
a local optimum of the original problem rather than of a surrogate.

Passivity assessment of the fitted model is a
[model-space operation](../Rational/README.md) invoked by the caller, not a
stage of the estimator.

## Order Search

When enabled, the pole count starts at a minimum and increases (in conjugate
pairs) until the target relative error is met, returning the lowest-order
model that satisfies it. This realizes the lowest-pole-count contract used
by the EMT line-fitting application.

An optional plateau stop ends the search once the best error has failed to
improve by a configured fraction for a configured number of consecutive
pole counts, so a structurally unreachable target does not force the full
ladder. An early stop carries the same verdict as an exhausted search.

## Restarts

When a fit stalls above the restart threshold, the current best poles are
perturbed by a deterministic bounded scaling and the iteration re-enters
Stage A. The perturbation sequence is reproducible; no random state exists.

## Initial Poles

Starting poles are logarithmically spaced over the sampled band as weakly
damped conjugate pairs, $p = -\alpha\beta \pm \mathrm{j}\beta$ with
$\alpha = 0.01$, following [1]. Real-only seeding and user-supplied starting
poles are supported; option combinations that conflict (custom weights
without the custom scheme, user poles without user seeding) are hard errors,
never silent precedence.

## Weighting

Weights multiply the squared residual of each sample block, as written in
the problem statement. Named schemes:

| Scheme | $w_m$ | Emphasis |
| ------ | ----- | -------- |
| `UNIFORM` | $1$ | absolute error |
| `INVERSE_MAGNITUDE` | $1 / \max(\lVert\mathbf{H}_m\rVert_F, \epsilon)$ | balanced |
| `INVERSE_SQUARED_MAGNITUDE` | $1 / \max(\lVert\mathbf{H}_m\rVert_F^2, \epsilon)$ | relative error |
| `CUSTOM` | user supplied | caller defined |

## Error Metrics

The fit statistics report absolute and relative root-mean-square error over
all samples and channels, the same metrics per channel, and the per-pass
pole displacement history. Failures of the eigenvalue computation or the
optimizer are reported as errors, never as silent convergence.

## Configuration

The module requires `GRIDKIT_ENABLE_EIGEN` (dense factorizations and
eigenvalue computations) and `GRIDKIT_ENABLE_IPOPT` (the nonlinear
refinement stage). Both dependencies are linked privately; public headers
expose only the standard library and the model space.

## Files

| File | Contents |
| ---- | -------- |
| `VectorFitting.hpp` | The fitting algorithm: parameters, statistics, `fit` |

## References

1. B. Gustavsen and A. Semlyen, "Rational approximation of frequency domain
   responses by vector fitting," IEEE Transactions on Power Delivery,
   vol. 14, no. 3, pp. 1052-1061, 1999.
2. B. Gustavsen, "Improving the pole relocating properties of vector
   fitting," IEEE Transactions on Power Delivery, vol. 21, no. 3,
   pp. 1587-1592, 2006.
3. D. Deschrijver, M. Mrozowski, T. Dhaene, and D. De Zutter,
   "Macromodeling of multiport systems using a fast implementation of the
   vector fitting method," IEEE Microwave and Wireless Components Letters,
   vol. 18, no. 6, pp. 383-385, 2008.
