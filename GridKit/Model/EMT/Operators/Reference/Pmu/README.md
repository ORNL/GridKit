# Pmu Model

`Pmu` is a phasor measurement unit on a three-phase bus. It measures the RMS
voltage and the active and reactive power, each through a first-order lag of
half a cycle of the system frequency base $f^\mathrm{base}$, and the angle and
frequency of the voltage with a phase-locked loop.

## Block Diagram

None.

## Model Parameters

None.

### Parameter Validation

None.

### Model Derived Parameters

```math
T = \dfrac{1}{2 f^\mathrm{base}}
```

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{v}$ | `v` | Input | [V] | Bus voltage, phases a, b and c | $\mathbf{v} \in \mathbb{R}^3$
$\mathbf{i}$ | `i` | Input | [A] | Current, phases a, b and c | $\mathbf{i} \in \mathbb{R}^3$, zero when left out
$v_\mathrm{rms}$ | `vrms` | Output | [V] | Line-to-line RMS voltage |
$p$ | `p` | Output | [W] | Active power | Internal variable
$q$ | `q` | Output | [var] | Reactive power | Internal variable
$\theta$ | `theta` | Output | [rad] | Angle |
$\omega$ | `omega` | Output | [rad/s] | Frequency |

The powers are positive in the direction of $\mathbf{i}$.

## Submodels

Symbol | Description | Type | Order | JSON | Inputs | Outputs
------ | ----------- | ---- | ----- | ---- | ------ | -------
$\mathrm{pll}$ | Phase-locked loop | [Pll](../Pll/README.md) | 2 | `pll` | $\mathbb{R}^3$ | $\mathbb{R}$

### Submodel Validation

None.

### Submodel Wiring

The phase-locked loop reads $\mathbf{v}$. The unit reads its angle $\theta$
and integrator state $x$.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$v_\mathrm{ms}$ | $[\mathrm{V}^2]$ | Mean-square voltage |
$p$ | [W] | Active power |
$q$ | [var] | Reactive power |

#### Algebraic

None.

### External Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}$ | [V] | Bus voltage vector owned by EMT bus | $\mathbf{v} \in \mathbb{R}^3$
$\mathbf{i}$ | [A] | Current vector | $\mathbf{i} \in \mathbb{R}^3$
$\theta$ | [rad] | Angle owned by the phase-locked loop |
$x$ | [rad/s] | Integrator state owned by the phase-locked loop |

#### Algebraic

None.

## Model Equations

The mean over the phases is removed from the voltages:

```math
\tilde{v}_n = v_n - \dfrac{1}{3}(v_a + v_b + v_c)
```

### Internal Equations

#### Differential

```math
\begin{aligned}
0 &= -\dfrac{\mathrm{d}v_\mathrm{ms}}{\mathrm{d}t} + \dfrac{1}{T}(\tilde{v}_a^2 + \tilde{v}_b^2 + \tilde{v}_c^2 - v_\mathrm{ms}) \\
0 &= -\dfrac{\mathrm{d}p}{\mathrm{d}t} + \dfrac{1}{T}(\tilde{v}_a i_a + \tilde{v}_b i_b + \tilde{v}_c i_c - p) \\
0 &= -\dfrac{\mathrm{d}q}{\mathrm{d}t} + \dfrac{1}{T}(\dfrac{1}{\sqrt{3}}((v_b - v_c)i_a + (v_c - v_a)i_b + (v_a - v_b)i_c) - q)
\end{aligned}
```

#### Algebraic

None.

### External Equations

```math
\begin{aligned}
v_\mathrm{rms} &\leftarrow \sqrt{v_\mathrm{ms}} \\
\theta &\leftarrow \theta^\mathrm{pll} \\
\omega &\leftarrow \omega^\mathrm{pll}
\end{aligned}
```

On a balanced bus, $v_\mathrm{rms}$ is the line-to-line RMS voltage. Its
derivative is taken as zero on a dead bus.

## Initialization

None.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`vms` | $[\mathrm{V}^2]$ | Mean-square voltage |
`p` | [W] | Active power |
`q` | [var] | Reactive power |
`vrms` | [V] | Line-to-line RMS voltage | Output
`theta` | [rad] | Angle | Output
`omega` | [rad/s] | Frequency | Output
