# PWM Model

`PWM` produces a three-phase switching signal from the $dq$ converter voltage
command. It limits the command to the voltage available from the DC supply,
returns the limited command to
[InnerCurrentControl](../InnerCurrentControl/README.md) for tracking
anti-windup, transforms the modulation command to phase coordinates, and
evaluates a smooth periodic switching function. Without inputs, it generates
sinusoidal PWM. Outputs are computed from the current input signals and time;
there are no internal DAE variables or residual rows. Enzyme generates the
output derivatives, which consumers compose through the signal Jacobian interface.

## Block Diagram

![PWM model block diagram](../../../../../../docs/Figures/EMT/Controller/PWM/diagram.png)

Figure 1: PWM model.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$M$ | [-] | `M` | Modulation index | Required without inputs, $M \in [0,1]$
$f_{\mathrm{m}}$ | [Hz] | `fm` | Modulation frequency | Required without inputs, positive
$f_{\mathrm{c}}$ | [Hz] | `fc` | Carrier frequency | Required, positive
$\alpha$ | [-] | `alignment` | Pulse alignment | Default $\frac{1}{2}$
$M^{\max}$ | [-] | `Mmax` | Sinusoidal modulation limit with voltage inputs | Default $1$, $0 < M^{\max} \le 1$

### Parameter Validation

```math
\begin{aligned}
f_{\mathrm{c}} &> 0 \\
0 &\le \alpha \le 1 \\
0 &< M^{\max} \le 1
\end{aligned}
```

Without inputs, the sinusoidal parameters also satisfy

```math
\begin{aligned}
0 &\le M \le 1 \\
f_{\mathrm{c}} &> f_{\mathrm{m}} > 0.
\end{aligned}
```

All parameters and the derived coefficient $a_u$ must be finite.

### Derived Parameters

```math
\begin{aligned}
\omega_{\mathrm{m}} &= 2\pi f_{\mathrm{m}} \\
\omega_{\mathrm{c}} &= 2\pi f_{\mathrm{c}} \\
T_{\mathrm{c}} &= \dfrac{2\pi}{\omega_{\mathrm{c}}}
                   = \dfrac{1}{f_{\mathrm{c}}} \\
a_u &= \dfrac{8}{3(M^{\max})^2} \\
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
$\mathbf{u}^{\mathrm{cmd}}$ | `u` | Input | [V] | Converter voltage command | Optional, $\mathbf{u}^{\mathrm{cmd}} \in \mathbb{R}^2$
$v_{\mathrm{dc}}$ | `vdc` | Input | [V] | DC voltage | With `u`, $v_{\mathrm{dc}} \ge 0$
$\theta$ | `theta` | Input | [rad] | Electrical reference angle | With `u`
$\mathbf{s}$ | `s` | Output | [-] | Three-phase switching function | $\mathbf{s} \in [0,1]^3$
$\mathbf{u}^{\mathrm{lim}}$ | `ulim` | Output | [V] | Limited voltage command | Requires `u`

The command uses $(d,q)$ order in the power-invariant
[Park](../../../Operators/Reference/Park/README.md) frame. Connect `u`,
`vdc`, and `theta` together, from `InnerCurrentControl.u`, an external DC
voltage, and `PLL.theta`. Return `ulim` to the current controller.

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

Connected voltage-command, DC, and angle variables may be differential.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{u}^{\mathrm{cmd}}$ | [V] | Converter voltage command | $\mathbf{u}^{\mathrm{cmd}} \in \mathbb{R}^2$
$v_{\mathrm{dc}}$ | [V] | DC voltage | $v_{\mathrm{dc}} \ge 0$
$\theta$ | [rad] | Electrical reference angle |

## Model Equations

With inputs, the limiter factor and modulation command are

```math
\begin{aligned}
\mathcal{L}_u(v_{\mathrm{dc}},\mathbf{u}^{\mathrm{cmd}}) &=
  \max\left(v_{\mathrm{dc}}^2,a_u\|\mathbf{u}^{\mathrm{cmd}}\|_2^2\right) \\
\mathbf{m}_{dq} &= \dfrac{2\mathbf{u}^{\mathrm{cmd}}}{\sqrt{\mathcal{L}_u(v_{\mathrm{dc}},\mathbf{u}^{\mathrm{cmd}})}}
\end{aligned}
```

The direction-preserving limit uses the CommonMath smooth
[`max`](../../../../../CommonMath.md#maximum) on squared voltages in
$\mathrm{V}^2$. In the power-invariant frame, a balanced phase peak of
$M^{\max}v_{\mathrm{dc}}/2$ corresponds to $\|\mathbf{u}^{\mathrm{cmd}}\|_2 =
\sqrt{3/8}\,M^{\max}v_{\mathrm{dc}}$, the configured linear modulation limit.
$M^{\max}=1$ gives the full sinusoidal PWM range without zero-sequence injection.[^spwm] The phase modulation command is the
inverse Park transform without zero sequence,

```math
\mathbf{m}_{abc} = \begin{bmatrix} m_a & m_b & m_c \end{bmatrix}^{\mathsf{T}}
= \mathbf{T}^{\mathsf{T}}(\theta)
\begin{bmatrix} m_d & m_q & 0 \end{bmatrix}^{\mathsf{T}},
```

so $|m_\ell| \le M^{\max}$. Without inputs,

```math
m_\ell(t)=M\sin\left(\omega_{\mathrm{m}}t+\phi_\ell\right).
```

For phase $\ell\in\{a,b,c\}$, the instantaneous duty ratio is

```math
d_\ell(t)=\dfrac{1+m_\ell(t)}{2}.
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
The command is continuous in time, with no carrier-synchronous sample-and-hold.
Single- or double-rate digital duty updates are distinct sampling models.[^sampling]

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

```math
\begin{aligned}
\mathbf{m}_{abc} &= \mathbf{T}^{\mathsf{T}}(\theta)
  \begin{bmatrix} m_d & m_q & 0 \end{bmatrix}^{\mathsf{T}} \\
\mathbf{u}^{\mathrm{lim}} &= \dfrac{v_{\mathrm{dc}}\mathbf{u}^{\mathrm{cmd}}}{\sqrt{\mathcal{L}_u(v_{\mathrm{dc}},\mathbf{u}^{\mathrm{cmd}})}} \\
s_\ell(t) &= S_\mu\left(t,d_\ell(t)\right),
\qquad \ell\in\{a,b,c\}
\end{aligned}
```

The modulation and limited-voltage equations apply only with voltage inputs, and
$\mathbf{u}^{\mathrm{lim}}=(v_{\mathrm{dc}}/2)\mathbf{m}_{dq}$.
In the hard-maximum limit, the command is unchanged inside the voltage limit
and keeps its direction at the limit magnitude outside it. CommonMath's smooth
maximum approaches this behavior continuously. Neither expression divides by
$v_{\mathrm{dc}}$: at zero DC voltage, the limited command is zero and the
modulation command remains finite. If the voltage command is also zero, both
are zero. Input values are evaluated at the current solver iterate.

### External Equations

None.

## Initialization

PWM owns no state. In [balanced initialization](../../../STATE.md#application),
`sa`, `sb`, `sc` are carrier means supplied by Converter. PWM transforms their
voltage to the PLL frame, checks the voltage limit, and requests `ud`, `uq`
from the controller. Otherwise, prescribed outputs must match the instantaneous
evaluated outputs. Limited-voltage outputs require voltage inputs.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`s` | [-] | Three-phase switching function | Expands to `sa`, `sb`, `sc`
`m` | [-] | Phase modulation command $\mathbf{m}_{abc}$ | Expands to `ma`, `mb`, `mc`
`ulim` | [V] | Limited voltage command | Expands to `ulimd`, `ulimq`; requires `u`

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

The numerical evaluation preserves the periodic pulse sum to roundoff.
For $\mu T_{\mathrm{c}} < 0.4$, the switching ripple is less than
$2.4\times10^{-20}$, so $S_\mu(t,d)=d$ to double precision.
Otherwise, reduce time to $[-T_{\mathrm{c}}/2,T_{\mathrm{c}}/2]$ and sum
$k=-R,\ldots,R$, where

```math
R=\begin{cases}
127, & 0.4\le\mu T_{\mathrm{c}}<4, \\
16, & 4\le\mu T_{\mathrm{c}}<16, \\
4, & \mu T_{\mathrm{c}}\ge16.
\end{cases}
```

The omitted sigmoid tails are below $10^{-21}$ for $d\in[0,1]$.

The broad-smoothing bound follows from the Fourier coefficients of the
periodized sigmoid kernel. With $x=2\pi^2/(\mu T_{\mathrm{c}})$,

```math
|S_\mu(t,d)-d|
\le \frac{4x}{\pi}\,
\frac{e^{-x}}{(1-e^{-x})(1-e^{-2x})}.
```

This bound selects a roundoff-equivalent evaluation; the resolved switching
function is evaluated directly from sigmoid pulses.

## Appendix C: Switching waveforms

![Centered sinusoidal switching signals](../../../../../../docs/Figures/EMT/Controller/PWM/waveforms.png)

Figure 2: Centered sinusoidal switching signals for $M=0.8$, $f_{\mathrm{m}}=60\,\mathrm{Hz}$, and $f_{\mathrm{c}}=900\,\mathrm{Hz}$ at $\mu^{-1}=0.005\,\mathrm{ms}$ and $\mu^{-1}=1\,\mathrm{ms}$.

[^spwm]: MathWorks, [*PWM Generator (Three-phase, Two-level)*](https://www.mathworks.com/help/sps/ref/pwmgeneratorthreephasetwolevel.html), "Overmodulation": $V_{\mathrm{peak}}=Mv_{\mathrm{dc}}/2$ for SPWM.
[^sampling]: Imperix, [*Carrier-based PWM*](https://imperix.com/doc/software/carrier-based-pwm), "PWM parameters update rate".
