# Delay Model

For input units $[u]$, `Delay` maps scalar input $u$ to delayed output
$y_{\mathrm{out}}$ through a chain of first-order lag sections.

> [!WARNING]
> The lag chain is an exact sampled $J$-step delay under forward Euler when
> $h = T$.
> Otherwise it is an approximation. The `fmax` parameter controls section
> density. It does not guarantee accuracy over a signal bandwidth.

## Block Diagram

![Delay operator block diagram](../../../../../../docs/Figures/EMT/Delay/diagram.png)

Figure 1: Delay model

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$\tau$ | [s] | `tau` | Total delay | Required, positive
$f_{\max}$ | [Hz] | `fmax` | Lag-chain section rate | Required, positive

### Parameter Validation

```math
\begin{aligned}
\tau &> 0 \\
f_{\max} &> 0
\end{aligned}
```

### Derived Parameters

```math
\begin{aligned}
J &= \lceil f_{\max}\tau \rceil \\
T &= \dfrac{\tau}{J} \\
\mathcal{J} &= \{1,\ldots,J\}
\end{aligned}
```

## Submodels

None.

### Submodel Validation

None.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{y}$ | $[u]$ | Section differential states | $\mathbf{y} \in \mathbb{R}^J$

#### Algebraic

None.

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$u$ | $[u]$ | Input signal | $u \in \mathbb{R}$

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$u$ | `input` | Input | $[u]$ | Input signal port | $u \in \mathbb{R}$
$y_{\mathrm{out}}$ | `out` | Output | $[u]$ | Delayed output port | $y_J$

## Model Equations

### Differential Equations

```math
\begin{aligned}
0 &= -T\dfrac{\mathrm{d}y_1}{\mathrm{d}t} - y_1 + u \\
0 &= -T\dfrac{\mathrm{d}y_\ell}{\mathrm{d}t} - y_\ell + y_{\ell-1},
     \quad \ell \in \mathcal{J} \setminus \{1\}
\end{aligned}
```

### Algebraic Equations

None.

### Wiring

```math
y_{\mathrm{out}} \leftarrow y_J
```

## Initialization

The uppercase symbols $U$, $Y_\ell$, and $Y_\mathrm{out}$ denote RMS phasors
of the corresponding lowercase variables.

### Input Initialization

```math
\begin{aligned}
U
  &\leftarrow \text{RMS input phasor} \\
u
  &\leftarrow \sqrt{2}\,\mathrm{Re}(U) \\
\dfrac{\mathrm{d}u}{\mathrm{d}t}
  &\leftarrow \sqrt{2}\,\mathrm{Re}(\mathrm{j}\omega_0U).
\end{aligned}
```

### Internal Initialization

```math
\begin{aligned}
Y_\ell
  &= (1+\mathrm{j}\omega_0T)^{-\ell}U,
  \quad \ell \in \mathcal{J}.
\end{aligned}
```

At $t=0$,

```math
\begin{aligned}
y_\ell
  &\leftarrow \sqrt{2}\,\mathrm{Re}(Y_\ell) \\
\dfrac{\mathrm{d}y_\ell}{\mathrm{d}t}
  &\leftarrow \sqrt{2}\,\mathrm{Re}(\mathrm{j}\omega_0Y_\ell),
  \quad \ell \in \mathcal{J}.
\end{aligned}
```

### Output Initialization

```math
\begin{aligned}
Y_\mathrm{out}
  &= Y_J \\
y_{\mathrm{out}}
  &\leftarrow \sqrt{2}\,\mathrm{Re}(Y_\mathrm{out}) \\
\dfrac{\mathrm{d}y_{\mathrm{out}}}{\mathrm{d}t}
  &\leftarrow \sqrt{2}\,\mathrm{Re}(\mathrm{j}\omega_0Y_\mathrm{out}).
\end{aligned}
```

## Monitors

None.
