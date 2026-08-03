# Yc Model

`Yc` computes the full-conductor characteristic admittance from the
per-unit-length series impedance and shunt admittance matrices.

Notes:
- With $\mathbf{Z}'=\mathbf{R}'+j\mathbf{X}'$ and
  $\mathbf{Y}'=\mathbf{G}'+j\mathbf{B}'$, the characteristic admittance
  of a coupled line satisfies the sandwich equation
  $\mathbf{Y}_c\mathbf{Z}'\mathbf{Y}_c=\mathbf{Y}'$. The commuting forms
  $\mathbf{Z}'\mathbf{Y}_c^2=\mathbf{Y}'$ and
  $\mathbf{Y}_c^2\mathbf{Z}'=\mathbf{Y}'$ agree with it only when
  $\mathbf{Z}'$ and $\mathbf{Y}'$ commute, as in the scalar and ideally
  transposed cases.
- Because $\mathbf{Z}'$ and $\mathbf{Y}'$ are symmetric for reciprocal
  lines, $\mathbf{Y}_c$ is exactly symmetric. The voltage- and
  current-wave conventions differ only in the propagation matrix,
  $\mathbf{Z}'\mathbf{Y}'$ versus
  $\mathbf{Y}'\mathbf{Z}'=(\mathbf{Z}'\mathbf{Y}')^\top$; the
  characteristic admittance is the same in both.
- The real and imaginary parts of $\mathbf{Y}_c$ are $\mathbf{G}_c$ and
  $\mathbf{B}_c$.

The characteristic admittance is decomposed as
$\mathbf{Y}_c=\mathbf{G}_c+j\mathbf{B}_c$.

The closed form is:

```math
\mathbf{Y}_c =
(\mathbf{R}' + j\mathbf{X}')^{-1}
\sqrt{
(\mathbf{R}' + j\mathbf{X}')
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

The real and imaginary parts of
$\mathbf{0}=-\mathbf{Y}'+\mathbf{Y}_c\mathbf{Z}'\mathbf{Y}_c$, split
over the product $\mathbf{P}=\mathbf{Z}'\mathbf{Y}_c$:

```math
\begin{aligned}
\mathbf{P}_r &= \mathbf{R}'\mathbf{G}_c - \mathbf{X}'\mathbf{B}_c \\
\mathbf{P}_i &= \mathbf{R}'\mathbf{B}_c + \mathbf{X}'\mathbf{G}_c \\
\mathbf{0} &= -\mathbf{G}'
  + \mathbf{G}_c\mathbf{P}_r
  - \mathbf{B}_c\mathbf{P}_i \\
\mathbf{0} &= -\mathbf{B}'
  + \mathbf{G}_c\mathbf{P}_i
  + \mathbf{B}_c\mathbf{P}_r
\end{aligned}
```

## Initialization

Validate the full-conductor matrix dimensions. Initialize from the
closed form
$\mathbf{Y}_c=\mathbf{Z}'^{-1}\sqrt{\mathbf{Z}'\mathbf{Y}'}$ at the
current inputs, which solves the algebraic equations exactly.

## Model Outputs

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{G}_c$ | [S] | Characteristic conductance | $\mathbb{R}^{K\times K}$
$\mathbf{B}_c$ | [S] | Characteristic susceptance | $\mathbb{R}^{K\times K}$
