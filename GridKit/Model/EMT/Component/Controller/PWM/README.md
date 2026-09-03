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
\dfrac{f_{\mathrm{c}}}{f_{\mathrm{m}}} &\in \mathbb{N}
\end{aligned}
```

### Derived Parameters

```math
m_f := \dfrac{f_{\mathrm{c}}}{f_{\mathrm{m}}}
```

The phase-offset vector is

```math
\boldsymbol{\phi}
:=
\begin{bmatrix}
0 & -\dfrac{2\pi}{3} & \dfrac{2\pi}{3}
\end{bmatrix}^{\mathsf{T}}
```

The normalized pulse half-widths are

```math
D_{k,m} := \dfrac{1}{4}
\left[
  1+M\sin\left(\dfrac{2\pi m}{m_f}+\phi_k\right)
\right]
\qquad
k\in\{a,b,c\},\quad m\in\{0,\ldots,m_f-1\}
```

The centered carrier phase is

```math
\theta(x)
:=
\operatorname{mod}\left(x+\dfrac{m_f}{2},m_f\right)
-\dfrac{m_f}{2}
```

The centered pulse function uses the GridKit
[`sigmoid`](../../../../../CommonMath.md#primitives)

```math
p(x,D) := \sigma(x+D)-\sigma(x-D)
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

```math
\begin{aligned}
s_a &\leftarrow \sum_{m=0}^{m_f-1}
  p\left(\theta\left(f_{\mathrm{c}}t-m-D_{a,m}\right),D_{a,m}\right) \\
s_b &\leftarrow \sum_{m=0}^{m_f-1}
  p\left(\theta\left(f_{\mathrm{c}}t-m-D_{b,m}\right),D_{b,m}\right) \\
s_c &\leftarrow \sum_{m=0}^{m_f-1}
  p\left(\theta\left(f_{\mathrm{c}}t-m-D_{c,m}\right),D_{c,m}\right)
\end{aligned}
```

## Initialization

None beyond the EMT initialization contract.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`s` | [-] | Three-phase switching function | $\mathbf{s} \in [0,1]^3$
