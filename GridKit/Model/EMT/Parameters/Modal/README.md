# Modal Decomposition

`ModalDecomposition` turns the propagation matrix $\Gamma(\omega)$, integrated
by the frequency sweep, into per-mode quantities: eigenvalues
$\lambda_m = \alpha_m + \mathrm{j}\beta_m$, the voltage and current
transformations $T_v$ and $T_i$, the propagation factors
$h_m = e^{-\lambda_m \ell}$, and the phase delays
$\tau_m = \ell\beta_m/\omega$.

## Why observation, not integration

The sweep DAE integrates only quantities analytic in $\omega$: the
per-unit-length matrices, $\Gamma$ through $\Gamma^2 = ZY$, and the
characteristic admittance. Individual eigenpairs are not analytic: their
identity is undefined at eigenvalue degeneracy, and symmetric tower
geometries make degeneracy exact and persistent. Continuing eigenvectors as
DAE states is therefore singular precisely where symmetric lines live.
Observing them pointwise at each emitted sample is exact everywhere, and
identity across samples becomes an explicit, documented rule instead of an
integrator's implicit one.

## The identity rule

Each `decompose()` call factors the current $\Gamma$ exactly, then:

1. **Assignment.** Identity follows eigenvectors, not eigenvalues: the
   previous sample's duals score every raw vector through the biorthogonal
   overlap $|T_i^H v|$, and each mode label claims its best match. Through a
   transversal eigenvalue crossing the vectors vary continuously while the
   value curves intersect, so labels never swap. On a grid coarser than a
   veering region this yields the diabatic labeling — vectors persist and
   the value traces cross cleanly — which is a convention, not a theorem,
   and the smoothest choice for rational fitting.

2. **Clusters.** Modes whose eigenvalues sit within the relative
   `cluster_gap` of one another form a cluster: only their invariant
   subspace is defined. The reported basis is the orthonormal span of the
   cluster (accurate even where the member vectors are not individually
   meaningful), rotated by the polar factor of its overlap with the
   previous frame — the closest admissible continuation. A simple mode is
   the one-by-one case of the same rule, its polar factor a phase, so a
   single alignment law covers every mode. Eigenvalues keep their own
   trajectories through a cluster by nearest-predecessor matching.

3. **Duals.** $T_i = T_v^{-H}$ is computed at every sample, never
   transported, so $T_i^H T_v = I$ holds to machine precision regardless of
   sweep tolerance.

The first sample fixes the canonical frame: modes ascend by
$(\mathrm{Im}\,\lambda, \mathrm{Re}\,\lambda)$ and each column's largest
entry is made real positive.

## Error contract

Presenting a cluster through its own eigenvalues bounds the defect of
$T_v\,\mathrm{diag}(e^{-\lambda_m\ell})\,T_i^H$ against $e^{-\Gamma\ell}$ by
$O(\ell\,|\gamma|\,\texttt{cluster\_gap})$; at exact degeneracy the defect
vanishes, which is exactly the configuration that deadlocks eigenvector
continuation. The default gap of $10^{-8}$ keeps the committed error orders
below the fitting targets downstream.

Known limitation: a cluster that is nearly defective (a genuine Jordan
block approached through non-normal veering) has an eigenvector span that
degrades before the invariant subspace does. The observed veering gaps of
overhead-line matrices sit far above the cluster threshold, so the cluster
path engages only for the semisimple degeneracies of symmetric geometry,
where the span is exact.

## Usage

Results live in one flat buffer described by `Signal` views, so a variable
monitor observes them exactly as it observes element states:

```cpp
ModalDecomposition<double, size_t> modes(k, length);
modes.decompose(omega, alpha, beta);   // per emitted sample, in order
monitor->setComplexMatrix(Variable::Tv, modes.tvReal(), modes.tvImag(), modes.data());
```

`reset()` forgets the tracked frame so the same object can serve a fresh
sweep.

## Files

| File | Contents |
| ---- | -------- |
| `ModalDecomposition.hpp` | The observation engine: parameters, `decompose`, signal views |
