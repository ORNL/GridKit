# **Forced Oscillation Source (ForcedOscillation)**

ForcedOscillation is a smooth signal source/operator for forced-oscillation
studies. It can generate a standalone oscillatory signal or add the forcing to
an attached input signal.

## Block Diagram

```text
source mode:
  time ---------------------> [ ForcedOscillation ] ---> output signal

additive mode:
  input signal + time ------> [ ForcedOscillation ] ---> output signal
```

## Model Parameters

Symbol                    | Units    | Description                         | Typical Value | Note
--------------------------|----------|-------------------------------------|---------------|-----
$A$                       | [p.u.]   | Oscillation amplitude               | 0.0           | JSON key: `A`
$f$                       | [Hz]     | Initial oscillation frequency       | 0.0           | JSON key: `f`
$K_f$                     | [Hz/s]   | Linear frequency ramp               | 0.0           | JSON key: `Kf`
$\Phi$                    | [rad]    | Phase offset                        | 0.0           | JSON key: `Phi`
$B$                       | [p.u.]   | DC bias                             | 0.0           | JSON key: `Bias`
$K_{\mathrm{in}}$         | [-]      | Input gain                          | 1.0           | JSON key: `Kin`
$u_0$                     | [p.u.]   | Standalone input value              | 0.0           | JSON key: `u0`; used when `input` is not attached
$T_{\mathrm{on}}$         | [sec]    | Start time                          | 0.0           | JSON key: `Ton`
$T_{\mathrm{off}}$        | [sec]    | Stop time                           | -1.0          | JSON key: `Toff`; negative means no stop time
$T_r$                     | [sec]    | Raised-cosine rise time             | 0.0           | JSON key: `Tr`
$T_f$                     | [sec]    | Raised-cosine fall time             | 0.0           | JSON key: `Tf`; active only when $T_{\mathrm{off}}\ge0$
$K_d$                     | [1/sec]  | Exponential decay rate              | 0.0           | JSON key: `Kd`
$L^{\min}$                | [p.u.]   | Output lower limit                  | -1.0e6        | JSON key: `Lmin`
$L^{\max}$                | [p.u.]   | Output upper limit                  | 1.0e6         | JSON key: `Lmax`

JSON keys are exact enum names. Units belong in this table and are not part of
the JSON key names.

### Parameter Validation

Invalid ForcedOscillation parameter sets are rejected by the following checks:

```math
\begin{aligned}
  &A \ge 0 \\
  &f \ge 0,\quad K_f \ge 0 \\
  &T_r \ge 0,\quad T_f \ge 0,\quad K_d \ge 0 \\
  &T_{\mathrm{off}} \ge 0 \Rightarrow T_{\mathrm{off}} \ge T_{\mathrm{on}} \\
  &L^{\min} \le L^{\max}
\end{aligned}
```

The output port is required. If the optional input port is present, it must be
linked to an assigned signal node.

## Model Variables

### Internal Variables

#### Algebraic

Symbol | Units  | Description
-------|--------|------------
$y$    | [p.u.] | Limited output signal

### External Variables

#### Algebraic

Symbol | Units  | Description
-------|--------|------------
$u$    | [p.u.] | Optional input signal

## Model Equations

Let $\tau = \max(t - T_{\mathrm{on}}, 0)$ and define the active-window
indicator:

```math
\begin{aligned}
  a(t) &=
    \begin{cases}
      1, & t \ge T_{\mathrm{on}}
           \text{ and } (T_{\mathrm{off}} < 0
           \text{ or } t < T_{\mathrm{off}}) \\
      0, & \text{otherwise}
    \end{cases}
\end{aligned}
```

The raised-cosine rise and fall factors are:

```math
\begin{aligned}
  r(t) &=
    \begin{cases}
      \dfrac{1-\cos\!\left(\pi(t-T_{\mathrm{on}})/T_r\right)}{2},
        & T_r > 0 \text{ and } T_{\mathrm{on}} \le t < T_{\mathrm{on}} + T_r \\
      1, & \text{otherwise}
    \end{cases} \\
  d(t) &=
    \begin{cases}
      \dfrac{1-\cos\!\left(\pi(T_{\mathrm{off}}-t)/T_f\right)}{2},
        & T_{\mathrm{off}} \ge 0,\ T_f > 0,\ \text{and } T_{\mathrm{off}}-T_f < t < T_{\mathrm{off}} \\
      1, & \text{otherwise}
    \end{cases} \\
  e(t) &= a(t)r(t)d(t)
\end{aligned}
```

The chirped phase and oscillatory forcing are:

```math
\begin{aligned}
  \theta(t) &= \Phi + 2\pi\left(f\tau + \dfrac{1}{2}K_f\tau^2\right) \\
  F(t) &= A e(t)\exp(-K_d\tau)\sin\!\left(\theta(t)\right)
\end{aligned}
```

With no attached input signal, use $u=u_0$. The raw and limited outputs are:

```math
\begin{aligned}
  v_{\mathrm{raw}} &= K_{\mathrm{in}}u + B + F(t) \\
  v_{\mathrm{lim}} &= \text{clamp}\!\left(v_{\mathrm{raw}}, L^{\min}, L^{\max}\right)
\end{aligned}
```

The single algebraic residual is:

```math
\begin{aligned}
  0 &= y - v_{\mathrm{lim}}
\end{aligned}
```

The output limiter uses GridKit's smooth
[Clamp](../../../../CommonMath.md#derived-functions). The component does not
configure IDA step size. Studies should cap max step size and output cadence
finely enough to resolve the forcing. Discontinuous waveforms are not part of
Phase 1.

## Initialization

Let $t_0$ denote the initial simulation time. The component initializes the
algebraic output consistently:

```math
\begin{aligned}
  u_{\mathrm{init}} &=
    \begin{cases}
      u(t_0), & \text{if the input signal is attached} \\
      u_0, & \text{otherwise}
    \end{cases} \\
  y_0 &=
    \text{clamp}\!\left(
      K_{\mathrm{in}}u_{\mathrm{init}} + B + F(t_0),
      L^{\min},
      L^{\max}
    \right) \\
  \dot y_0 &= 0
\end{aligned}
```

## Model Outputs

Output   | Units  | Description
---------|--------|------------
`in`     | [p.u.] | Input value used by the model
`env`    | [-]    | Envelope value $e(t)$
`force`  | [p.u.] | Oscillatory forcing $F(t)$
`vraw`   | [p.u.] | Unlimited output value
`out`    | [p.u.] | Limited output signal
`active` | [-]    | Active-window indicator
