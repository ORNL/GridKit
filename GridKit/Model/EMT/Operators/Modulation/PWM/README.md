# PWM Model

`PWM` maps duty ratios to smooth switching functions for $N$ channels.
Channels share carrier timing. Voltage limiting, coordinate transformations,
and duty generation belong to the upstream controller. Any limited-voltage
feedback for antiwindup remains there.

## Block Diagram

![PWM operator block diagram](../../../../../../docs/Figures/EMT/PWM/diagram.png)

Figure 1: PWM model; the input bubble denotes $1-\mathbf{d}$ componentwise,
and the stacked pulse branches sum over all carrier periods

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$N$ | [-] | `N` | Number of channels | Required, positive integer
$f_{\mathrm{c}}$ | [Hz] | `fc` | Carrier frequency | Required, positive
$\alpha$ | [-] | `alignment` | Pulse alignment | Default $\frac{1}{2}$

### Parameter Validation

All real parameters must be finite.

```math
\begin{aligned}
N &\in \mathbb{Z}_{>0} \\
f_{\mathrm{c}} &> 0 \\
0 &\le \alpha \le 1
\end{aligned}
```

### Model Derived Parameters

```math
T = \dfrac{1}{f_{\mathrm{c}}}
```

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{d}$ | `d` | Input | [-] | Duty ratios | $\mathbf{d} \in \mathbb{R}^N$
$\mathbf{s}$ | `s` | Output | [-] | Switching functions | $\mathbf{s} \in \mathbb{R}^N$

Each channel satisfies $0 \le d_n \le 1$ and $0 \le s_n \le 1$.
The switching functions can drive
[Bridge](../../../Component/Converter/Bridge/README.md).

## Submodels

None.

### Submodel Validation

None.

### Submodel Wiring

None.

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

None.

### External Variables

The tables apply to inputs that reference DAE variables owned by connected
models. Prescribed and computed signals introduce no variables.

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{d}$ | [-] | Duty ratios | Differential input, $\mathbf{d} \in \mathbb{R}^N$

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{d}$ | [-] | Duty ratios | Algebraic input, $\mathbf{d} \in \mathbb{R}^N$

## Model Equations

Smooth functions: [$\sigma$](../../../../../CommonMath.md#logistic-function).

Scalar offsets and smooth functions apply componentwise.
For each integer $k$, the pulse-edge vectors are

```math
\begin{aligned}
\mathbf{a}_k &= [k+\alpha(1-\mathbf{d})]T \\
\mathbf{b}_k &= \mathbf{a}_k+T\mathbf{d}
\end{aligned}
```

All channels use the same carrier periods starting at $t=kT$.
Electrical phase and phase sequence are carried by $\mathbf{d}$.
Every periodic copy uses the current duty vector. There is no sample-and-hold.
Alignment $\alpha=0$ fixes the rising edge at the period start,
$\alpha=1$ fixes the falling edge at the period end, and
$\alpha=\frac{1}{2}$ centers the pulse.

### Internal Equations

#### Differential

None.

#### Algebraic

None.

### External Equations

```math
\mathbf{s} \leftarrow \sum_{k\in\mathbb{Z}}
  \left[\sigma(t-\mathbf{a}_k)-\sigma(t-\mathbf{b}_k)\right]
```

## Initialization

None.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`s` | [-] | Switching functions | $\mathbf{s} \in \mathbb{R}^N$

## Development

For fixed duty, periodic summation preserves the carrier-period mean:

```math
\dfrac{1}{T}\int_0^T \mathbf{s}(t)\mathrm{d}t = \mathbf{d}
```

The shared [CommonMath smoothing scale](../../../../../CommonMath.md#primitives)
$\mu$ selects resolution within the same model. Small $\mu T$ approaches
instantaneous duty; large $\mu T$ resolves switching. The isolated-edge width is

```math
\Delta t_{10\text{–}90} = \dfrac{2\ln 9}{\mu}
```

![PWM switching functions](../../../../../../docs/Figures/EMT/PWM/waveforms.png)

Figure 2: Centered switching functions for sinusoidal duty inputs with modulation
index $\frac{4}{5}$, frequency $60\,\mathrm{Hz}$, and carrier frequency $900\,\mathrm{Hz}$
at $\mu^{-1}=5\,\mu\mathrm{s}$ and $\mu^{-1}=1\,\mathrm{ms}$
