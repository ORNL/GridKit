# Propagation Model

For input units $[u]$, `Propagation` is the $K$-channel current-form propagation
operator used by `LineDistributed`. It applies a fitted input factor, one scalar
delay per mode, and a fitted output factor while preserving the input units.

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

The full modal basis has one mode per channel:

```math
\begin{aligned}
M &= K \\
\boldsymbol{\tau} &= [\tau_1,\ldots,\tau_M]^\mathsf T
\end{aligned}
```

The modal delays are produced by the offline propagation fitting and enter
through the delay bank's `delays` coefficient set.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{u}$ | `input` | Input | $[u]$ | Input vector port | $\mathbf{u} \in \mathbb{R}^K$
$\mathbf{y}$ | `out` | Output | $[u]$ | Output vector port | $\mathbf{y} \in \mathbb{R}^K$

## Submodels

Symbol | Description | Type | Order | JSON | Inputs | Outputs
------ | ----------- | ---- | ----- | ---- | ------ | -------
$\mathbf{g}_\mathrm{in}$ | Input factor | [VectorFit](../../Rational/VectorFit/README.md) | $KQ_{\mathbf{g}_\mathrm{in}}$ | `input` | $\mathbb{R}^K$ | $\mathbb{R}^M$
$\mathbf{d}$ | Modal delay bank | [Delay](../Delay/README.md) | History | `delays` | $\mathbb{R}^M$ | $\mathbb{R}^M$
$\mathbf{g}_\mathrm{out}$ | Output factor | [VectorFit](../../Rational/VectorFit/README.md) | $KQ_{\mathbf{g}_\mathrm{out}}$ | `output` | $\mathbb{R}^M$ | $\mathbb{R}^K$


### Submodel Validation

Both rational factors must have stable poles and no term linear in $s$. The
input factor has algebraic input; the output factor has differential input
from the modal delay bank.

```math
\mathbf{E}^\mathrm{in}=\mathbf{E}^\mathrm{out}=\mathbf{0}
```

### Submodel Wiring

```math
\begin{aligned}
\mathbf{w} &\leftarrow \mathbf{g}_\mathrm{in}[\mathbf{u}] \\
\mathbf{z} &\leftarrow \mathbf{d}[\mathbf{w}]
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
\mathbf{y} \leftarrow \mathbf{g}_\mathrm{out}[\mathbf{z}]
```

## Initialization

TBD

## Monitors

None.
