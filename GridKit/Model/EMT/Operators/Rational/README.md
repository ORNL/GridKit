# Rational Operators

Rational operators represent fitted transfer functions with real runtime
states. The shared `Rational` implementation evaluates residuals, initializes
memory states, and assembles the exact linear Jacobian, including gradients of
computed algebraic input signals. It supports rectangular input/output maps;
three-phase components use the same implementation with three inputs and
three outputs.

## Models

- [VectorFit](VectorFit/README.md): general residue matrices, $KQ$ real states.
- [StateSpace](StateSpace/README.md): factorized residues, $Q$ real states.

An input derivative is required only when its column of $E$ is nonzero.
Algebraic inputs require a zero column. Singular $E$ is supported without
inversion. Consumers impose their own physical constraints on the coefficients.
