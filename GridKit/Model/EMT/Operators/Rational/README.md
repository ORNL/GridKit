# Rational Operators

Rational operators represent fitted transfer functions with real states and
rectangular input/output maps.

## Models

- [VectorFit](VectorFit/README.md): general residue matrices, $KQ$ real states.
- [StateSpace](StateSpace/README.md): factorized residues, $Q$ real states.

An input derivative is required only when its column of $\mathbf{E}$ is nonzero.
Algebraic inputs require a zero column. Singular $\mathbf{E}$ is supported without
inversion. Consumers impose their own physical constraints on the coefficients.

An output may be read through `output()` and `appendOutputGradient()` without
a residual destination. A bound output destination receives the operator
contribution in its equation. Computed current outputs require zero derivative
feedthrough; derivative-dependent outputs contribute to residual equations.
