# Zc Model

`Zc` computes the full-conductor characteristic impedance from the
per-unit-length series impedance and shunt admittance matrices.

Notes:
- With $\mathbf{Z}'=\mathbf{R}'+j\mathbf{X}'$ and
  $\mathbf{Y}'=\mathbf{G}'+j\mathbf{B}'$, the characteristic impedance
  of a coupled line satisfies the sandwich equation
  $\mathbf{Z}_c\mathbf{Y}'\mathbf{Z}_c=\mathbf{Z}'$. The commuting forms
  $\mathbf{Y}'\mathbf{Z}_c^2=\mathbf{Z}'$ and
  $\mathbf{Z}_c^2\mathbf{Y}'=\mathbf{Z}'$ agree with it only when
  $\mathbf{Z}'$ and $\mathbf{Y}'$ commute, as in the scalar and ideally
  transposed cases.
- Because $\mathbf{Z}'$ and $\mathbf{Y}'$ are symmetric for reciprocal
  lines, $\mathbf{Z}_c$ is exactly symmetric. The voltage- and
  current-wave conventions differ only in the propagation matrix,
  $\mathbf{Z}'\mathbf{Y}'$ versus
  $\mathbf{Y}'\mathbf{Z}'=(\mathbf{Z}'\mathbf{Y}')^\top$; the
  characteristic impedance is the same in both.
- The real and imaginary parts of $\mathbf{Z}_c$ are $\mathbf{R}_c$ and
  $\mathbf{X}_c$.

The characteristic impedance is decomposed as
$\mathbf{Z}_c=\mathbf{R}_c+j\mathbf{X}_c$.

The closed form is:

```math
\mathbf{Z}_c =
(\mathbf{G}' + j\mathbf{B}')^{-1}
\sqrt{
(\mathbf{G}' + j\mathbf{B}')
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

The real and imaginary parts of
$\mathbf{0}=-\mathbf{Z}'+\mathbf{Z}_c\mathbf{Y}'\mathbf{Z}_c$, split
over the product $\mathbf{P}=\mathbf{Y}'\mathbf{Z}_c$:

```math
\begin{aligned}
\mathbf{P}_r &= \mathbf{G}'\mathbf{R}_c - \mathbf{B}'\mathbf{X}_c \\
\mathbf{P}_i &= \mathbf{G}'\mathbf{X}_c + \mathbf{B}'\mathbf{R}_c \\
\mathbf{0} &= -\mathbf{R}'
  + \mathbf{R}_c\mathbf{P}_r
  - \mathbf{X}_c\mathbf{P}_i \\
\mathbf{0} &= -\mathbf{X}'
  + \mathbf{R}_c\mathbf{P}_i
  + \mathbf{X}_c\mathbf{P}_r
\end{aligned}
```

## Initialization

Validate the full-conductor matrix dimensions. Initialize from the
closed form
$\mathbf{Z}_c=\mathbf{Y}'^{-1}\sqrt{\mathbf{Y}'\mathbf{Z}'}$ at the
current inputs, which solves the algebraic equations exactly.

## Model Outputs

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{R}_c$ | [$\Omega$] | Characteristic resistance | $\mathbb{R}^{K\times K}$
$\mathbf{X}_c$ | [$\Omega$] | Characteristic reactance | $\mathbb{R}^{K\times K}$
