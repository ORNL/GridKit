# Signal sources

## Constant signal source (`ConstantSignalSource`)

`ConstantSignalSource` publishes a constant complex value through separate real
and imaginary output ports. Its values have no intrinsic units because the
consuming signal defines their physical meaning.

| Symbol | Units | JSON | Description | Default | Note |
| --- | --- | --- | --- | --- | --- |
| $S_r$ | unspecified | `Sr` | Real output value | 0.0 | Optional |
| $S_i$ | unspecified | `Si` | Imaginary output value | 0.0 | Optional |

The `sr` and `si` output ports publish $S_r$ and $S_i$, respectively.

## Forced-oscillation signal source (`ForcedOscillation`)

`ForcedOscillation` publishes a stateless, time-dependent periodic waveform
with an optional linear frequency ramp, raised-cosine activation window, and
exponential decay. A sine and smooth square, triangle, and sawtooth
approximations share the same phase and envelope model.

### Notes

- The output is exogenous: it is not a system DAE variable and its signal node
  carries `INVALID_INDEX`.
- The consuming port defines the physical units and base of the output and of
  the amplitude $A$.
- A negative $T_{\mathrm{off}}$ disables deactivation and the fall window.

### Block Diagram

```text
time ──> [ ForcedOscillation ] ──> output
```

### Model Parameters

| Symbol | Units | JSON | Description | Default | Note |
| --- | --- | --- | --- | --- | --- |
| $A$ | output-signal units | `A` | Oscillation amplitude | 0.0 | Optional |
| $f$ | Hz | `f` | Initial oscillation frequency | 0.0 | Optional |
| $K_f$ | Hz/s | `Kf` | Linear frequency ramp | 0.0 | Optional |
| $\Phi$ | rad | `Phi` | Phase offset | 0.0 | Optional |
| $T_{\mathrm{on}}$ | s | `Ton` | Activation time | 0.0 | Optional |
| $T_{\mathrm{off}}$ | s | `Toff` | Deactivation time | -1.0 | Optional; negative disables deactivation |
| $T_r$ | s | `Tr` | Raised-cosine rise time | 0.0 | Optional |
| $T_f$ | s | `Tf` | Raised-cosine fall time | 0.0 | Optional; used only for nonnegative $T_{\mathrm{off}}$ |
| $K_d$ | 1/s | `Kd` | Exponential decay rate | 0.0 | Optional |
| $W$ | integer | `waveform` | Carrier-waveform selector | 0 | Optional; 0 = sine, 1 = square, 2 = triangle, 3 = sawtooth |

Integer and real serialized values are accepted for every parameter. The
`waveform` value must represent an exact integer.

#### Parameter Validation

Every supplied parameter must be finite. The following domains and relationship
are enforced:

```math
\begin{aligned}
  A &\ge 0, & f &\ge 0, & K_f &\ge 0, \\
  T_r &\ge 0, & T_f &\ge 0, & K_d &\ge 0, \\
  W &\in \{0,1,2,3\}, \\
  T_{\mathrm{off}} &< 0
    \quad\text{or}\quad T_{\mathrm{off}} \ge T_{\mathrm{on}}.
\end{aligned}
```

#### Model Derived Parameters

None.

### Model Ports

The model has no bus or signal-input ports.

| Name | Port | Init | Description |
| --- | --- | --- | --- |
| Waveform output | `output` | Known | Required time-dependent exogenous signal |

### Model Variables

#### Internal Variables

##### Differential

None.

##### Algebraic

None. The published value is source-owned signal storage, not a DAE variable.

#### External Variables

##### Differential

None.

##### Algebraic

None.

### Model Equations

Let elapsed activation time be

```math
\tau(t) = \max\!\left(t-T_{\mathrm{on}},0\right).
```

The active-window indicator is

```math
a(t) =
\begin{cases}
1, & t \ge T_{\mathrm{on}}
     \text{ and }
     \left(T_{\mathrm{off}} < 0 \text{ or } t < T_{\mathrm{off}}\right), \\
0, & \text{otherwise}.
\end{cases}
```

The raised-cosine rise and fall factors are

```math
\begin{aligned}
r(t) &=
\begin{cases}
\dfrac{1-\cos\!\left(\pi(t-T_{\mathrm{on}})/T_r\right)}{2},
  & T_r>0 \text{ and } T_{\mathrm{on}} \le t<T_{\mathrm{on}}+T_r, \\
1, & \text{otherwise},
\end{cases} \\
d(t) &=
\begin{cases}
\dfrac{1-\cos\!\left(\pi(T_{\mathrm{off}}-t)/T_f\right)}{2},
  & T_{\mathrm{off}}\ge0,\quad T_f>0,
    \quad\text{and } T_{\mathrm{off}}-T_f<t<T_{\mathrm{off}}, \\
1, & \text{otherwise},
\end{cases} \\
e(t) &= a(t)r(t)d(t).
\end{aligned}
```

The non-sinusoidal carriers use fixed smoothing constants

```math
\kappa=3, \qquad \rho_{\triangle}=0.98, \qquad \rho_{\mathrm{saw}}=0.9.
```

The carrier selected by $W$ is

```math
q_W(\theta)=
\begin{cases}
\sin(\theta), & W=0, \\
\dfrac{\tanh(\kappa\sin(\theta))}{\tanh(\kappa)}, & W=1, \\
\dfrac{\arcsin(\rho_{\triangle}\sin(\theta))}
      {\arcsin(\rho_{\triangle})}, & W=2, \\
\dfrac{\operatorname{atan2}(\rho_{\mathrm{saw}}\sin(\theta),
                             1+\rho_{\mathrm{saw}}\cos(\theta))}
      {\arcsin(\rho_{\mathrm{saw}})}, & W=3.
\end{cases}
```

Because $0<\rho_{\triangle},\rho_{\mathrm{saw}}<1$, all four carriers are
smooth and periodic. Each crosses zero with positive slope at phase zero and
is normalized so that $|q_W|\le1$ with both unit extrema attained. Thus $A$ is
the peak carrier amplitude before the envelope and decay are applied. The
square and triangle extrema occur at $\theta=\pi/2$ and $3\pi/2$; the sawtooth
extrema occur where $\cos(\theta)=-\rho_{\mathrm{saw}}$. The phase and
published waveform are

```math
\begin{aligned}
\theta(t) &= \Phi + 2\pi\left(f\tau(t)+\frac{1}{2}K_f\tau(t)^2\right), \\
  s_{\mathrm{FO}}(t) &= A e(t)\exp\!\left(-K_d\tau(t)\right)q_W\!\left(\theta(t)\right).
\end{aligned}
```

#### Internal Equations

##### Differential

None.

##### Algebraic

None. `updateTime()` evaluates $s_{\mathrm{FO}}(t)$ directly.

#### External Equations

None.

### Initialization

#### Input Initialization

None.

#### Internal Initialization

None.

#### Output Initialization

At model initialization time $t_0$,

```math
s_{\mathrm{FO}} \leftarrow s_{\mathrm{FO}}(t_0).
```

Each subsequent `updateTime(t, alpha)` call republishes
$s_{\mathrm{FO}}(t)$.

### Monitorable Outputs

| Output | Units | Description | Note |
| --- | --- | --- | --- |
| `output` | output-signal units | Published waveform $s_{\mathrm{FO}}$ | Same value as the output signal |
| `envelope` | - | Raised-cosine activation envelope $e$ | Includes the active-window indicator |
| `active` | - | Active-window indicator $a$ | 0.0 or 1.0 |

### Testing

- `validation()` checks parameter types and domains, including waveform-selector
  values.
- `waveform()` checks stationary and chirped, decaying output.
- `carrierWaveforms()` checks the phase convention and amplitude of all four
  carriers.
- `activationWindowAndMonitors()` checks the rise/fall envelope and monitor values.
- `dependencyTracking()` checks the exogenous output and empty DAE structure.
