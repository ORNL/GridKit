# Propagation Model

For input units $[u]$, `Propagation` is the $K$-channel current-form propagation
operator used by `LineDistributed`. It applies a fitted input factor, one scalar
delay per mode, and a fitted output factor while preserving the input units.

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
\mathcal{M} &= \{1,\ldots,M\} \\
\boldsymbol{\tau} &= [\tau_1,\ldots,\tau_M]^\mathsf T
\end{aligned}
```

## Submodels

Symbol | Description | Type | Order | JSON | Inputs | Outputs
------ | ----------- | ---- | ----- | ---- | ------ | -------
$\mathbf{g}_\mathrm{in}$ | Input factor | [VectorFit](../../Rational/VectorFit/README.md) | $KQ_{\mathbf{g}_\mathrm{in}}$ | `input` | $\mathbb{R}^K$ | $\mathbb{R}^M$
$d_m$, $m \in \mathcal{M}$ | Modal delay | [Delay](../Delay/README.md) | History | `delays` | $\mathbb{R}$ | $\mathbb{R}$
$\mathbf{g}_\mathrm{out}$ | Output factor | [VectorFit](../../Rational/VectorFit/README.md) | $KQ_{\mathbf{g}_\mathrm{out}}$ | `output` | $\mathbb{R}^M$ | $\mathbb{R}^K$

The offline fitting targets and propagation factorization are

```math
\begin{aligned}
\mathbf{G}^\mathrm{in}(s)
  &\approx \mathbf{H}^\mathrm{mps}(s)\mathbf{T}_i^{-1}(s) \\
\mathbf{G}^\mathrm{out}(s) &\approx \mathbf{T}_i(s) \\
\mathbf{H}^\mathrm{mps}(s)
  &= \mathrm{diag}(h_1^\mathrm{mps}(s),\ldots,h_M^\mathrm{mps}(s)) \\
\mathbf{D}_{\boldsymbol{\tau}}(s)
  &= \mathrm{diag}(e^{-s\tau_1},\ldots,e^{-s\tau_M}) \\
\mathbf{H}(s)
  &= \mathbf{T}_i(s)\mathbf{D}_{\boldsymbol{\tau}}(s)
     \mathbf{H}^\mathrm{mps}(s)\mathbf{T}_i^{-1}(s) \\
  &\approx \mathbf{G}^\mathrm{out}(s)\mathbf{D}_{\boldsymbol{\tau}}(s)
     \mathbf{G}^\mathrm{in}(s)
\end{aligned}
```

$\mathbf{H}^\mathrm{mps}$ is the diagonal modal minimum-phase-shift propagation
function with the modal delays removed. The current modal transformation
$\mathbf{T}_i$ maps modal currents to phase coordinates, and
$\mathbf{T}_i^{-1}$ maps phase currents to modal coordinates.

### Submodel Validation

Both rational factors must have stable poles and no term linear in $s$. The
input factor has algebraic input; the output factor has differential input
from the modal delays.

```math
\mathbf{E}^\mathrm{in}=\mathbf{E}^\mathrm{out}=\mathbf{0}
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

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{u}$ | `input` | Input | $[u]$ | Input vector port | $\mathbf{u} \in \mathbb{R}^K$
$\mathbf{y}$ | `out` | Output | $[u]$ | Output vector port | $\mathbf{y} \in \mathbb{R}^K$

## Model Equations

### Differential Equations

None.

### Algebraic Equations

None.

### Wiring

```math
\begin{aligned}
\mathbf{w} &\leftarrow \mathbf{g}_\mathrm{in}[\mathbf{u}] \\
z_m &\leftarrow d_m[w_m] \\
\mathbf{z} &\leftarrow [z_1,\ldots,z_M]^\mathsf T \\
\mathbf{y} &\leftarrow \mathbf{g}_\mathrm{out}[\mathbf{z}]
\end{aligned}
```

## Initialization

### Input Initialization

```math
\begin{aligned}
\widehat{\mathbf{u}}
  &\leftarrow \text{RMS input phasor} \\
\mathbf{u}
  &\leftarrow \sqrt{2}\,\mathrm{Re}(\widehat{\mathbf{u}}) \\
\dfrac{\mathrm{d}\mathbf{u}}{\mathrm{d}t}
  &\leftarrow \sqrt{2}\,\mathrm{Re}(s_0\widehat{\mathbf{u}}).
\end{aligned}
```

### Internal Initialization

The input-factor output phasor is

```math
\widehat{\mathbf{w}}
  = \mathbf{G}^\mathrm{in}(s_0)\widehat{\mathbf{u}}.
```

The input factor initializes from $\widehat{\mathbf{u}}$. Each modal delay
initializes from $\widehat{w}_m$ and provides $\widehat{z}_m$; the output
factor initializes from $\widehat{\mathbf{z}}$.

### Output Initialization

```math
\begin{aligned}
\widehat{\mathbf{y}}
  &= \mathbf{G}^\mathrm{out}(s_0)\widehat{\mathbf{z}} \\
\mathbf{y}
  &\leftarrow \sqrt{2}\,\mathrm{Re}(\widehat{\mathbf{y}}) \\
\dfrac{\mathrm{d}\mathbf{y}}{\mathrm{d}t}
  &\leftarrow \sqrt{2}\,\mathrm{Re}(s_0\widehat{\mathbf{y}}).
\end{aligned}
```

## Monitors

None.
