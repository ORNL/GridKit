# PWM Model

`PWM` produces a three-phase sinusoidal PWM switching signal. The model adds no
DAE variables or residual rows.

## Block Diagram

![PWM model switching signal](../../../../../../docs/Figures/EMT/Controller/PWM/diagram.png)

Figure 1: PWM switching signal for $M=0.8$, $f_{\mathrm{m}}=60\,\mathrm{Hz}$, and $f_{\mathrm{c}}=900\,\mathrm{Hz}$ at $\mu=240$ and $\mu=1$

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
f_{\mathrm{c}} &> f_{\mathrm{m}} > 0 \\
\dfrac{f_{\mathrm{c}}}{f_{\mathrm{m}}} &\in 3\mathbb{N}
\end{aligned}
```

### Derived Parameters

```math
\begin{aligned}
T_{\mathrm{m}} &:= \dfrac{1}{f_{\mathrm{m}}} \\
m_f &:= \dfrac{f_{\mathrm{c}}}{2f_{\mathrm{m}}}
\end{aligned}
```

The normalized pulse half-width is

```math
D_m := \dfrac{1}{4}
\left[
  1+M\sin\left(\dfrac{\pi m}{m_f}\right)
\right]
\qquad
m\in\{0,\ldots,2m_f-1\}
```

The pulse-aligned carrier phase is

```math
\theta_m(x)
:=
\operatorname{mod}\left(f_{\mathrm{c}}x-m-D_m+m_f,2m_f\right)
-m_f
```

The switching function uses the GridKit
[`sigmoid`](../../../../../CommonMath.md#primitives)

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

```math
\begin{aligned}
s_a(t) &\leftarrow \sum_{m=0}^{2m_f-1}
  \left[
    \sigma\left(\theta_m(t)+D_m\right)
    -\sigma\left(\theta_m(t)-D_m\right)
  \right] \\
s_b\left(t+\dfrac{T_{\mathrm{m}}}{3}\right) &\leftarrow \sum_{m=0}^{2m_f-1}
  \left[
    \sigma\left(\theta_m(t)+D_m\right)
    -\sigma\left(\theta_m(t)-D_m\right)
  \right] \\
s_c\left(t-\dfrac{T_{\mathrm{m}}}{3}\right) &\leftarrow \sum_{m=0}^{2m_f-1}
  \left[
    \sigma\left(\theta_m(t)+D_m\right)
    -\sigma\left(\theta_m(t)-D_m\right)
  \right]
\end{aligned}
```

## Initialization

None beyond the EMT initialization contract.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`s` | [-] | Three-phase switching function | $\mathbf{s} \in [0,1]^3$
