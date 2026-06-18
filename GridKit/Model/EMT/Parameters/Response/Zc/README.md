# Zc Model

`Zc` computes the full-conductor characteristic impedance from the
per-unit-length series impedance and shunt admittance matrices.

Notes:
- For multiconductor lines, the V-definition and I-definition differ by the
  ordering of the matrix products.
- This model uses the V-definition: the shunt matrix pair multiplies
  $\mathbf{Z}_c^2$ on the left.
- The I-definition places the shunt matrix pair on the right.
- The real and imaginary parts of $\mathbf{Z}_c$ are $\mathbf{R}_c$ and
  $\mathbf{X}_c$.

The characteristic impedance is decomposed as
$\mathbf{Z}_c=\mathbf{R}_c+j\mathbf{X}_c$.

The complex form typically given for this convention is:

```math
\mathbf{Z}_c =
\sqrt{
(\mathbf{G}' + j\mathbf{B}')^{-1}
(\mathbf{R}' + j\mathbf{X}')
}
```

## Model Parameters

None. `Zc` combines full-conductor parameter outputs for $K$ physical
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
$\mathbf{R}_c$ | [$\Omega$] | Characteristic resistance | $\mathbb{R}^{K\times K}$
$\mathbf{X}_c$ | [$\Omega$] | Characteristic reactance | $\mathbb{R}^{K\times K}$

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
\mathbf{0} &= -\mathbf{R}'
  + \mathbf{G}'\left(\mathbf{R}_c\mathbf{R}_c
  - \mathbf{X}_c\mathbf{X}_c\right)
  - \mathbf{B}'\left(\mathbf{R}_c\mathbf{X}_c
  + \mathbf{X}_c\mathbf{R}_c\right) \\
\mathbf{0} &= -\mathbf{X}'
  + \mathbf{G}'\left(\mathbf{R}_c\mathbf{X}_c
  + \mathbf{X}_c\mathbf{R}_c\right)
  + \mathbf{B}'\left(\mathbf{R}_c\mathbf{R}_c
  - \mathbf{X}_c\mathbf{X}_c\right)
\end{aligned}
```

## Initialization

Validate the full-conductor matrix dimensions. Initialize
$\mathbf{R}_c$ and $\mathbf{X}_c$ from the current $\mathbf{R}'$,
$\mathbf{X}'$, $\mathbf{G}'$, and $\mathbf{B}'$ inputs.

## Model Outputs

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{R}_c$ | [$\Omega$] | Characteristic resistance | $\mathbb{R}^{K\times K}$
$\mathbf{X}_c$ | [$\Omega$] | Characteristic reactance | $\mathbb{R}^{K\times K}$
