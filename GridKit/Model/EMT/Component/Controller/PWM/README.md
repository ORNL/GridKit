# PWM Model

`PWM` produces a three-phase switching signal from a sampled modulation input.
Without an input, it generates sinusoidal PWM. The model adds no DAE variables
or residual rows.

## Block Diagram

![PWM model switching signal](../../../../../../docs/Figures/EMT/Controller/PWM/diagram.png)

Figure 1: Sampled PWM interface and centered sinusoidal switching signals for $M=0.8$, $f_{\mathrm{m}}=60\,\mathrm{Hz}$, and $f_{\mathrm{c}}=900\,\mathrm{Hz}$ at $\mu^{-1}=0.005\,\mathrm{ms}$ and $\mu^{-1}=1\,\mathrm{ms}$.

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
f_{\mathrm{c}} &> f_{\mathrm{m}} > 0 \\
\dfrac{f_{\mathrm{c}}}{f_{\mathrm{m}}} &\in 3\mathbb{N}.
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

The modulation input is sampled at accepted boundaries; it introduces no
continuous DAE dependency.

#### Differential

None.

#### Algebraic

None.

## Model Equations

For phase $\ell\in\{a,b,c\}$ and carrier interval $k\in\mathbb{Z}$,
regular sampling holds the modulation command for one carrier period:

```math
\begin{aligned}
t_k &= kT_{\mathrm{c}} \\
m_{\ell,k} &= m_\ell(t_{k-1}^-),
\qquad -1 \le m_{\ell,k} \le 1.
\end{aligned}
```

The command sampled at the accepted carrier boundary $t_{k-1}$ is applied over
interval $k$, a one-carrier computational delay. Trial residual evaluations
and interpolated monitor samples do not change the held command. Without an
input, the prescribed sinusoid retains its alignment-dependent sample:

```math
m_{\ell,k}
= M\sin\left(\omega_{\mathrm{m}}(k+\alpha)T_{\mathrm{c}}+\phi_\ell\right).
```

The full duty ratio and switching instants are

```math
\begin{aligned}
d_{\ell,k} &= \dfrac{1+m_{\ell,k}}{2} \\
t_{\ell,k}^{\mathrm{on}}
&= \left[k+\alpha(1-d_{\ell,k})\right]T_{\mathrm{c}} \\
t_{\ell,k}^{\mathrm{off}}
&= \left[k+\alpha+(1-\alpha)d_{\ell,k}\right]T_{\mathrm{c}}.
\end{aligned}
```

The switching function uses the GridKit
[`sigmoid`](../../../../../CommonMath.md#primitives) with shared sharpness
$\mu>0$.

The isolated-edge width and harmonic attenuation of a periodically repeated
pulse are

```math
\begin{aligned}
\Delta t_{10\text{–}90} &= \dfrac{2\ln 9}{\mu} \\
A(f,\mu) &= \dfrac{2\pi^2 f/\mu}{\sinh(2\pi^2 f/\mu)}
\end{aligned}
```

$\mu$ | $\Delta t_{10\text{–}90}$ | Interpretation at $f_{\mathrm{c}}=900\,\mathrm{Hz}$
----- | ------------------------- | ---------------------------------------------------------
$240$ | $18.3\,\mathrm{ms}$ | Broad smoothing; switching suppressed
$50000$ | $87.9\,\mathrm{\mu s}$ | Resolved smoothed switching with sufficiently fine steps
$200000$ | $22.0\,\mathrm{\mu s}$ | Sharper edges; finer steps required

Set solver `mu` before model construction. It also affects other CommonMath
primitives. Monitor spacing alone does not establish integration accuracy;
check switching harmonics against the sampled-edge prediction.


### Internal Equations

#### Differential

None.

#### Algebraic

None.

### External Equations

With a modulation input, each held command defines a periodic pulse train,

```math
p_{\ell,k}(t)
=
\sum_{r\in\mathbb{Z}}
\left[
  \sigma\left(t-t_{\ell,k}^{\mathrm{on}}-rT_{\mathrm{c}}\right)
  -\sigma\left(t-t_{\ell,k}^{\mathrm{off}}-rT_{\mathrm{c}}\right)
\right],
```

and the switching function crossfades the trains at the carrier boundaries
with the smoothed interval indicator:

```math
s_\ell(t)
\leftarrow
\sum_{k\in\mathbb{Z}}
\left[\sigma_{\mathrm{c}}(t-t_k)-\sigma_{\mathrm{c}}(t-t_{k+1})\right]
p_{\ell,k}(t),
\qquad
\sigma_{\mathrm{c}}(x)=\dfrac{1}{2}\left[1+\tanh\left(\dfrac{\mu_{\mathrm{c}}x}{2}\right)\right],
\qquad
\mu_{\mathrm{c}}=\max\left(\mu,\ \dfrac{2\ln(4/\varepsilon)}{T_{\mathrm{c}}}\right).
```

The crossfade rate equals $\mu$ unless that would spread a transition beyond
half a carrier period; the floor keeps the hold causal with the one-carrier
delay, so the weights over $t_k \le t < t_{k+1}$ involve only the committed
trains $k-1$, $k$, and $k+1$. The switching function is therefore smooth
everywhere and the integrator never restarts. Smoothing preserves the period
mean of each train:

```math
\dfrac{1}{T_{\mathrm{c}}}
\int_{t_k}^{t_{k+1}}p_{\ell,k}(t)\,\mathrm{d}t
=d_{\ell,k}.
```

Thus $A(nf_{\mathrm{c}},\mu)$ attenuates the carrier harmonics while retaining
the commanded mean. Each train approaches $d_{\ell,k}$ as $\mu$ decreases.
While consecutive commands differ, the crossfade shifts the interval mean of
$s_\ell$ from $d_{\ell,k}$ by at most $2\ln 2/(\mu_{\mathrm{c}}T_{\mathrm{c}})$;
equal neighbouring commands leave it exact.
The implementation evaluates a train in the time domain or as its Fourier
series, whichever needs fewer terms; both agree to rounding.

Without a modulation input, pulses retain their prescribed sinusoidal samples:

```math
s_\ell(t)
\leftarrow
\sum_{k\in\mathbb{Z}}
\left[
  \sigma\left(t-t_{\ell,k}^{\mathrm{on}}\right)
  -\sigma\left(t-t_{\ell,k}^{\mathrm{off}}\right)
\right],
\qquad
\ell\in\{a,b,c\}
```

For sampled input operation, the solver stops at each sampling instant to
commit the sample and continues without restarting. Its maximum step is
bounded by $\min(T_{\mathrm{c}}/20,\mu^{-1})$ to resolve carrier edges.

## Initialization

With a modulation input, the initialized input supplies the commands of the
previous, current, and next carrier intervals. Only these three commands are
retained; no modulation prehistory is required. Starting a new study resets
them. A restart within a carrier interval retains them.
Without an input, the sinusoidal switching sequence supplies its own prehistory.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`s` | [-] | Three-phase switching function | $\mathbf{s} \in [0,1]^3$

See [case connections](../../../INPUT_FORMAT.md#case-connections) for vector signal wiring.
