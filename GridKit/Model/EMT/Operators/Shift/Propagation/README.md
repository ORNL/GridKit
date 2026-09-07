# Propagation Model

For input units $[u]$, `Propagation` is the $K$-channel current-form propagation
operator used by `LineDistributed`. Each mode applies a proper rational matrix
and a scalar transport delay. The mode outputs are summed, preserving the
input units.

```math
\begin{aligned}
\mathbf{H}(s)
  &= \sum_{m=1}^M \mathbf{H}^\mathrm{mps}_m(s) e^{-s\tau_m}
\end{aligned}
```

## Block Diagram

![Propagation operator block diagram](../../../../../../docs/Figures/EMT/Propagation/diagram.png)

Figure 1: Propagation model

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$K$ | [-] | `K` | Signal dimension | Required, positive integer

### Parameter Validation

```math
K \in \mathbb{Z}_{>0}
```

### Derived Parameters

Each entry of `modes` supplies `tau` and a `K`-by-`K` VectorFit coefficient
set `H`. Modes with the same delay may be grouped into one matrix; consequently
the number of fitted delay groups need not equal the channel count.

```math
\begin{aligned}
M &= \operatorname{size}(\texttt{modes}) \\
\boldsymbol{\tau} &= [\tau_1,\ldots,\tau_M]^\mathsf T
\end{aligned}
```

The modal delays are produced by the offline propagation fitting. Each scalar
delay is applied to all $K$ outputs of its corresponding matrix.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{u}$ | `input` | Input | $[u]$ | Input vector port | $\mathbf{u} \in \mathbb{R}^K$
$\mathbf{y}$ | `out` | Output | $[u]$ | Output vector port | $\mathbf{y} \in \mathbb{R}^K$

## Submodels

Symbol | Description | Type | Order | JSON | Inputs | Outputs
------ | ----------- | ---- | ----- | ---- | ------ | -------
$\mathbf{h}^{\mathrm{mps}}_m$ | Rational part of mode $m$ | [VectorFit](../../Rational/VectorFit/README.md) | $KQ_m$ | `modes[m].H` | $\mathbb{R}^K$ | $\mathbb{R}^K$
$\mathbf{d}_m$ | Transport delay of mode $m$ | [Delay](../Delay/README.md) | $K$ algebraic rows and history | `modes[m].tau` | $\mathbb{R}^K$ | $\mathbb{R}^K$


### Submodel Validation

Every rational matrix must have stable poles and no term linear in $s$.
There must be at least one mode, and all delays must be finite and positive.

```math
\mathbf{E}_m=\mathbf{0},\qquad
\operatorname{Re}(p_{mq})<0,\qquad \tau_m>0
```

### Submodel Wiring

```math
\begin{aligned}
\mathbf{w}_m &\leftarrow \mathbf{h}^{\mathrm{mps}}_m[\mathbf{u}] \\
\mathbf{z}_m &\leftarrow \mathbf{d}_m[\mathbf{w}_m]
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

None.

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{u}$ | $[u]$ | Input vector | $\mathbf{u} \in \mathbb{R}^K$

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

None.

### External Equations

```math
\mathbf{y} \leftarrow \sum_{m=1}^M \mathbf{z}_m
```

## Initialization

For a supplied constant or harmonic input prehistory, initialize each rational
matrix with the same input value, derivative, and angular frequency. Its
frequency response gives the filtered prehistory supplied to the delay bank.
The two propagation directions of a line use independent instances.

For $\omega>0$, write
$\widehat{\mathbf{u}}=\mathbf{u}(t_0)-j\mathbf{u}'(t_0)/\omega$. Then

```math
\begin{aligned}
\widehat{\mathbf{w}}_m
  &\leftarrow \mathbf{H}^{\mathrm{mps}}_m(j\omega)\widehat{\mathbf{u}} \\
\mathbf{y}(t_0)
  &\leftarrow \operatorname{Re}\!\left(
    \sum_{m=1}^M e^{-j\omega\tau_m}\widehat{\mathbf{w}}_m\right).
\end{aligned}
```

At $\omega=0$, the supplied derivative must vanish and the prehistory is
constant. A zero prehistory represents an initially unenergized line.

## Monitors

None.
