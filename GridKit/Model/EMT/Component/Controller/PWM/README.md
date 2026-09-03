# PWM Model

`PWM` produces a three-phase sinusoidal PWM switching signal. The model adds no
DAE variables or residual rows.

## Block Diagram

![PWM model switching signal](../../../../../../docs/Figures/EMT/Controller/PWM/diagram.png)

Figure 1: PWM switching signal for $M=0.8$, $f_{\mathrm{m}}=60\,\mathrm{Hz}$, and $f_{\mathrm{c}}=900\,\mathrm{Hz}$ at $\mu=4f_{\mathrm{c}}f_{\mathrm{m}}$ and $\mu=f_{\mathrm{c}}$

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$M$ | [-] | `M` | Modulation index | Required, $M \in [0,1]$
$f_{\mathrm{m}}$ | [Hz] | `fm` | Modulation frequency | Required, positive
$f_{\mathrm{c}}$ | [Hz] | `fc` | Carrier frequency | Required, $f_{\mathrm{c}}>f_{\mathrm{m}}$

### Parameter Validation

```math
\begin{aligned}
0 &\le M \le 1 \\
f_{\mathrm{c}} &> f_{\mathrm{m}} > 0
\end{aligned}
```

### Derived Parameters

*Important!*: We can and should smooth approximate sawtooth in a cheap way.

```math
\begin{aligned}
\operatorname{st}(x) &:= x-\text{floor}(x) \\
\tau_m &:= t-\dfrac{\operatorname{st}(f_{\mathrm{c}}t)-m}{f_{\mathrm{c}}}
\qquad m\in\{-1,0,1\}
\end{aligned}
```

The three-phase pulse-width vector is

```math
\mathbf{d}(t)
=
\dfrac{1}{2f_{\mathrm{c}}}
\begin{bmatrix}
1 + M\sin\left(2\pi f_{\mathrm{m}}t\right) \\
1 + M\sin\left(2\pi f_{\mathrm{m}}t-\dfrac{2\pi}{3}\right) \\
1 + M\sin\left(2\pi f_{\mathrm{m}}t+\dfrac{2\pi}{3}\right)
\end{bmatrix}
```

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{s}$ | `s` | Output | [-] | Three-phase switching function | $\mathbf{s} \in [0,1]^3$

## Submodels

None.

### Submodel Validation

None.

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

None.

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

None.

### External Equations

$\sigma$ denotes the GridKit [`sigmoid`](../../../../../CommonMath.md#primitives).

```math
\begin{aligned}
s_a &\leftarrow \sum_{m=-1}^{1}
  \left[
    \sigma\left(t-\tau_m\right)
    -\sigma\left(t-\tau_m-d_a\left(\tau_m\right)\right)
  \right] \\
s_b &\leftarrow \sum_{m=-1}^{1}
  \left[
    \sigma\left(t-\tau_m\right)
    -\sigma\left(t-\tau_m-d_b\left(\tau_m\right)\right)
  \right] \\
s_c &\leftarrow \sum_{m=-1}^{1}
  \left[
    \sigma\left(t-\tau_m\right)
    -\sigma\left(t-\tau_m-d_c\left(\tau_m\right)\right)
  \right]
\end{aligned}
```

## Initialization

None beyond the EMT initialization contract.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`s` | [-] | Three-phase switching function | $\mathbf{s} \in [0,1]^3$
