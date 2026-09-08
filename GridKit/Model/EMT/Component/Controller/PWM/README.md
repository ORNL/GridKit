# PWM Model

`PWM` produces a three-phase switching signal from a continuous modulation input.
Without an input, it generates sinusoidal PWM. The model adds no DAE variables
or residual rows.

## Block Diagram

![PWM model switching signal](../../../../../../docs/Figures/EMT/Controller/PWM/diagram.png)

Figure 1: Continuous PWM interface and centered sinusoidal switching signals for $M=0.8$, $f_{\mathrm{m}}=60\,\mathrm{Hz}$, and $f_{\mathrm{c}}=900\,\mathrm{Hz}$ at $\mu^{-1}=0.005\,\mathrm{ms}$ and $\mu^{-1}=1\,\mathrm{ms}$.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$M$ | [-] | `M` | Modulation index | Required without `m`, $M \in [0,1]$
$f_{\mathrm{m}}$ | [Hz] | `fm` | Modulation frequency | Required without `m`, positive
$f_{\mathrm{c}}$ | [Hz] | `fc` | Carrier frequency | Required, positive
$\alpha$ | [-] | `alignment` | Pulse alignment | Default $\frac{1}{2}$

### Parameter Validation

```math
\begin{aligned}
f_{\mathrm{c}} &> 0 \\
0 &\le \alpha \le 1
\end{aligned}
```

Without a modulation input, the sinusoidal parameters also satisfy

```math
\begin{aligned}
0 &\le M \le 1 \\
f_{\mathrm{c}} &> f_{\mathrm{m}} > 0.
\end{aligned}
```

### Derived Parameters

```math
\begin{aligned}
\omega_{\mathrm{m}} &= 2\pi f_{\mathrm{m}} \\
\omega_{\mathrm{c}} &= 2\pi f_{\mathrm{c}} \\
T_{\mathrm{c}} &= \dfrac{2\pi}{\omega_{\mathrm{c}}}
                   = \dfrac{1}{f_{\mathrm{c}}} \\
\boldsymbol{\phi}
&=
\begin{bmatrix}
\phi_a & \phi_b & \phi_c
\end{bmatrix}^{\mathsf{T}}
=
\begin{bmatrix}
0 & -\dfrac{2\pi}{3} & \dfrac{2\pi}{3}
\end{bmatrix}^{\mathsf{T}}
\end{aligned}
```

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{m}$ | `m` | Input | [-] | Three-phase modulation command | Optional, $\mathbf{m} \in [-1,1]^3$
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

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{m}$ | [-] | Modulation command | Differential-input configuration

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{m}$ | [-] | Modulation command | Algebraic-input configuration

## Model Equations

For phase $\ell\in\{a,b,c\}$, the instantaneous duty ratio is

```math
d_\ell(t)=\dfrac{1+m_\ell(t)}{2}.
```

Without a modulation input,

```math
m_\ell(t)=M\sin\left(\omega_{\mathrm{m}}t+\phi_\ell\right).
```

For duty argument $d\in[0,1]$, the periodic pulse edges and switching function are

```math
\begin{aligned}
a_k(d)&=\left[k+\alpha(1-d)\right]T_{\mathrm{c}} \\
b_k(d)&=\left[k+\alpha+(1-\alpha)d\right]T_{\mathrm{c}} \\
S_\mu(t,d)&=\sum_{k\in\mathbb{Z}}
\left[\sigma\left(t-a_k(d)\right)-\sigma\left(t-b_k(d)\right)\right].
\end{aligned}
```

Here $\sigma$ is the CommonMath
[`sigmoid`](../../../../../CommonMath.md#primitives) with sharpness $\mu>0$.
Every term uses the current duty argument; $k$ indexes periodic copies.

> [!NOTE]
> $\mu$ selects simulation resolution within the same continuous model:
> small values approach the instantaneous duty; large values resolve switching.
> The carrier-period mean is exactly $d$ for a fixed duty argument (Appendix A).

The isolated-edge width is

```math
\Delta t_{10\text{–}90}=\dfrac{2\ln 9}{\mu}.
```

Set solver `mu` before model construction. It also affects other CommonMath
primitives. The maximum integration step is $\mu^{-1}$ to resolve sharp edges;
monitor spacing does not determine integration accuracy.

### Internal Equations

#### Differential

None.

#### Algebraic

None.

### External Equations

```math
s_\ell(t)\leftarrow S_\mu\left(t,d_\ell(t)\right),
\qquad \ell\in\{a,b,c\}.
```

The outputs are algebraic expressions without owned DAE variables. Input values
and their derivatives are evaluated at the current solver iterate.

## Initialization

Evaluate the switching function from the initialized modulation input and time.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`s` | [-] | Three-phase switching function | $\mathbf{s} \in [0,1]^3$

See [case connections](../../../INPUT_FORMAT.md#case-connections) for vector signal wiring.

## Appendix A: Mean and smoothing limits

Let $a=\alpha(1-d)T_{\mathrm{c}}$ and
$b=[\alpha+(1-\alpha)d]T_{\mathrm{c}}$ for a fixed duty argument $d$.
Periodic summation gives

```math
\begin{aligned}
\dfrac{1}{T_{\mathrm{c}}}\int_0^{T_{\mathrm{c}}}S_\mu(t,d)\,\mathrm{d}t
&=\dfrac{1}{T_{\mathrm{c}}}\int_{-\infty}^{\infty}
\left[\sigma(t-a)-\sigma(t-b)\right]\,\mathrm{d}t \\
&=\dfrac{b-a}{T_{\mathrm{c}}}=d.
\end{aligned}
```

The sigmoid derivative $g_\mu(x)=\mathrm{d}\sigma(x)/\mathrm{d}x$ is
nonnegative with unit integral.
Thus $S_\mu$ is a periodic rectangular pulse convolved with $g_\mu$:
it is smooth and satisfies $0\le S_\mu\le1$.
As $\mu T_{\mathrm{c}}\to0$, the periodized kernel approaches
$T_{\mathrm{c}}^{-1}$; as $\mu T_{\mathrm{c}}\to\infty$, it concentrates
at each pulse edge. Consequently,

```math
\begin{aligned}
\lim_{\mu T_{\mathrm{c}}\to0}S_\mu(t,d)&=d, \\
\lim_{\mu T_{\mathrm{c}}\to\infty}S_\mu(t,d)
&=\sum_{k\in\mathbb{Z}}\mathbf{1}_{(a_k(d),b_k(d))}(t).
\end{aligned}
```

The switching limit holds away from the edges. For varying modulation,
$s_\ell(t)\to d_\ell(t)$ under broad smoothing. The fixed-duty mean identity
does not assert an exact carrier-period mean for a changing command.

## Appendix B: Pulse-sum evaluation

Outside a pulse interval, let $\delta$ be the distance to its nearest edge.
With $r=e^{-\mu dT_{\mathrm{c}}}$ and $z=e^{-\mu\delta}$, its contribution is

```math
\sigma(t-a_k)-\sigma(t-b_k)
=\dfrac{z(1-r)}{(1+z)(1+zr)}.
```

The nearest pulse is evaluated directly. Each outward replica updates
$z\leftarrow z e^{-\mu T_{\mathrm{c}}}$, avoiding repeated sigmoid evaluations.
The same pulse window and compensated summation are retained; this recurrence
is an algebraic rearrangement of the pulse sum.
