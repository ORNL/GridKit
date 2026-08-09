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

`ForcedOscillation` publishes a stateless, time-dependent sinusoid with an
optional linear frequency ramp, raised-cosine activation window, and
exponential decay.

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

Integer and real serialized values are accepted for every parameter.

#### Parameter Validation

Every supplied parameter must be finite. The following domains and relationship
are enforced:

```math
\begin{aligned}
  A &\ge 0, & f &\ge 0, & K_f &\ge 0, \\
  T_r &\ge 0, & T_f &\ge 0, & K_d &\ge 0, \\
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

The phase and published waveform are

```math
\begin{aligned}
\theta(t) &= \Phi + 2\pi\left(f\tau(t)+\frac{1}{2}K_f\tau(t)^2\right), \\
s_{\mathrm{FO}}(t) &= A e(t)\exp\!\left(-K_d\tau(t)\right)\sin\!\left(\theta(t)\right).
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

None.
