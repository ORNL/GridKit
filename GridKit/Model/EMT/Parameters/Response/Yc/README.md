# Yc Model

`Yc` computes the full-conductor characteristic admittance from the
per-unit-length series impedance and shunt admittance matrices.

Notes:
- For multiconductor lines, the V-definition and I-definition differ by the
  ordering of the matrix products.
- This model uses the V-definition: the series matrix pair multiplies
  $\mathbf{Y}_c^2$ on the left.
- The I-definition places the series matrix pair on the right.
- The real and imaginary parts of $\mathbf{Y}_c$ are $\mathbf{G}_c$ and
  $\mathbf{B}_c$.

The characteristic admittance is decomposed as
$\mathbf{Y}_c=\mathbf{G}_c+j\mathbf{B}_c$.

The complex form typically given for this convention is:

```math
\mathbf{Y}_c =
\sqrt{
(\mathbf{R}' + j\mathbf{X}')^{-1}
(\mathbf{G}' + j\mathbf{B}')
}
```

## Model Parameters

None. `Yc` combines full-conductor parameter outputs for $K$ physical
conductors.

### Parameter Validation

The full-conductor inputs must be square and have the same dimension.

### Model Derived Parameters

None.

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{G}_c$ | [S] | Characteristic conductance | $\mathbb{R}^{K\times K}$
$\mathbf{B}_c$ | [S] | Characteristic susceptance | $\mathbb{R}^{K\times K}$

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{R}'$ | [$\Omega$/m] | Full-conductor resistance | from [`SeriesImpedance`](../../Effects/SeriesImpedance/README.md)
$\mathbf{X}'$ | [$\Omega$/m] | Full-conductor reactance | $\omega\mathbf{L}'$ from [`SeriesImpedance`](../../Effects/SeriesImpedance/README.md)
$\mathbf{G}'$ | [S/m] | Full-conductor shunt conductance | from [`ShuntAdmittance`](../../Effects/ShuntAdmittance/README.md)
$\mathbf{B}'$ | [S/m] | Full-conductor shunt susceptance | $\omega\mathbf{C}'$ from [`ShuntAdmittance`](../../Effects/ShuntAdmittance/README.md)

## Model Equations

### Differential Equations

None.

### Algebraic Equations

```math
\begin{aligned}
\mathbf{0} &= -\mathbf{G}'
  + \mathbf{R}'\left(\mathbf{G}_c\mathbf{G}_c
  - \mathbf{B}_c\mathbf{B}_c\right)
  - \mathbf{X}'\left(\mathbf{G}_c\mathbf{B}_c
  + \mathbf{B}_c\mathbf{G}_c\right) \\
\mathbf{0} &= -\mathbf{B}'
  + \mathbf{R}'\left(\mathbf{G}_c\mathbf{B}_c
  + \mathbf{B}_c\mathbf{G}_c\right)
  + \mathbf{X}'\left(\mathbf{G}_c\mathbf{G}_c
  - \mathbf{B}_c\mathbf{B}_c\right)
\end{aligned}
```

## Initialization

Validate the full-conductor matrix dimensions. Initialize
$\mathbf{G}_c$ and $\mathbf{B}_c$ from the current $\mathbf{R}'$,
$\mathbf{X}'$, $\mathbf{G}'$, and $\mathbf{B}'$ inputs.

## Model Outputs

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{G}_c$ | [S] | Characteristic conductance | $\mathbb{R}^{K\times K}$
$\mathbf{B}_c$ | [S] | Characteristic susceptance | $\mathbb{R}^{K\times K}$
