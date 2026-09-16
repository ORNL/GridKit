# Transformer Model

`Transformer` represents a three-phase bank of two-winding transformers with
independent nonlinear $\pi$ cores. Internal quantities use transformer per-unit
bases; terminal voltages and injected currents use SI units.

## Block Diagram

![Transformer model circuit diagram](../../../../../docs/Figures/EMT/Transformer/diagram.png)

Figure 1: One phase of the per-unit winding circuit; terminal connection maps
are omitted. Resistors in the shunt branches are labeled by conductance;
nonlinear inductors carry $\beta_e g(\psi_e)$.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$N$ | [-] | `N` | Number of phases | Fixed at $3$
$S_\mathrm{b}$ | [VA] | `S` | Rated three-phase apparent power | Required, positive
$V_1$ | [V] | `V1` | Rated line-to-line RMS voltage of winding 1 | Required, positive
$V_2$ | [V] | `V2` | Rated line-to-line RMS voltage of winding 2 | Required, positive
$f_\mathrm{b}$ | [Hz] | `f` | Rated frequency | Required, positive
$\mathbf{P}_1$ | [-] | `P1` | Terminal 1 connection map | Default $\mathbf{I}_N$
$\mathbf{P}_2$ | [-] | `P2` | Terminal 2 connection map | Default $\mathbf{I}_N$
$\tau$ | [p.u.] | `tap` | Off-nominal ratio on winding 1 | Default $1$
$R$ | [p.u.] | `R` | Short-circuit resistance | Nonnegative
$X$ | [p.u.] | `X` | Short-circuit reactance | Required, positive
$I_0$ | [p.u.] | `I0` | No-load current at rated voltage | Required, positive
$P_0$ | [W] | `P0` | No-load loss at rated voltage | Nonnegative
$\psi_\mathrm{K}$ | [p.u.] | `knee` | Knee flux linkage | Default $1.2$
$L_\mathrm{sat}$ | [p.u.] | `Lsat` | Magnetizing saturation inductance | Positive, defaults to $2X$
$\beta$ | [-] | `split` | Winding 1 share of the magnetizing branch | Default $0.5$

### Parameter Validation

All parameters and derived coefficients must be finite.

```math
\begin{aligned}
N &= 3 \\
S_\mathrm{b}, V_1, V_2, f_\mathrm{b} &> 0 \\
\tau &> 0 \\
R &\ge 0 \\
X &> 0 \\
P_0 &\ge 0 \\
I_0 &> P_0 / S_\mathrm{b} \\
\psi_\mathrm{K}, L_\mathrm{sat} &> 0 \\
0 \le \beta &\le 1 \\
\mathbf{P}_1, \mathbf{P}_2 &\in \{-1, 0, 1\}^{N \times N}
\end{aligned}
```

Each column of $\mathbf{P}_e$ has one entry $+1$, at most one entry $-1$, and
zeros elsewhere. Winding-voltage magnitudes must agree across columns.
Ungrounded neutrals are not represented; $L_\mathrm{m}>L_\mathrm{sat}$.

### Derived Parameters

Use $\mathcal{N}=\{1,2,3\}$ in $(a,b,c)$ order and $e\in\{1,2\}$.
For $\boldsymbol{\gamma}=[0,-2\pi/3,2\pi/3]^\mathsf{T}$,
$\alpha_n=\exp(\mathrm{j}\gamma_n)$.

```math
\begin{aligned}
\omega_\mathrm{b} &= 2\pi f_\mathrm{b} \\
V_{\mathrm{w},e} &= \dfrac{V_e}{\sqrt{3}}\left|(\mathbf{P}_e^\mathsf{T}\boldsymbol{\alpha})_1\right| \\
V_{\mathrm{pk},e} &= \sqrt{2}\,V_{\mathrm{w},e} \\
I_{\mathrm{pk},e} &= \dfrac{\sqrt{2}\,S_\mathrm{b}}{3\,V_{\mathrm{w},e}} \\
R_1 = R_2 &= R/2 \\
G_\mathrm{c} &= P_0 / S_\mathrm{b} \\
I_\mathrm{m} &= \sqrt{I_0^2 - G_\mathrm{c}^2} \\
L_\mathrm{m} &= 1 / I_\mathrm{m} \\
\beta_1 = \beta, \qquad \beta_2 &= 1 - \beta
\end{aligned}
```

No-load data calibrate the shunt branches; terminal values also depend on
series impedance. Hysteresis and capacitances are neglected.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{v}_1$ | `v1` | Input | [V] | Terminal 1 bus voltage | $\mathbf{v}_1 \in \mathbb{R}^N$
$\mathbf{v}_2$ | `v2` | Input | [V] | Terminal 2 bus voltage | $\mathbf{v}_2 \in \mathbb{R}^N$
$\mathbf{i}_1$ | `i1` | Output | [A] | Current injection at terminal 1 | $\mathbf{i}_1 \in \mathbb{R}^N$
$\mathbf{i}_2$ | `i2` | Output | [A] | Current injection at terminal 2 | $\mathbf{i}_2 \in \mathbb{R}^N$

## Submodels

None.

### Submodel Validation

None.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}_{12}$ | [p.u.] | Series leakage current from node 1 to node 2 | $\mathbf{i}_{12} \in \mathbb{R}^N$
$\boldsymbol{\psi}_1$ | [p.u.] | Node 1 magnetizing flux linkage | $\boldsymbol{\psi}_1 \in \mathbb{R}^N$
$\boldsymbol{\psi}_2$ | [p.u.] | Node 2 magnetizing flux linkage | $\boldsymbol{\psi}_2 \in \mathbb{R}^N$

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{e}_1$ | [p.u.] | Node 1 voltage | $\mathbf{e}_1 \in \mathbb{R}^N$
$\mathbf{e}_2$ | [p.u.] | Node 2 voltage | $\mathbf{e}_2 \in \mathbb{R}^N$
$\mathbf{i}_\mathrm{m1}$ | [p.u.] | Node 1 magnetizing current | $\mathbf{i}_\mathrm{m1} \in \mathbb{R}^N$
$\mathbf{i}_\mathrm{m2}$ | [p.u.] | Node 2 magnetizing current | $\mathbf{i}_\mathrm{m2} \in \mathbb{R}^N$
$\mathbf{i}_\mathrm{w1}$ | [p.u.] | Winding 1 current | Positive into the winding, $\mathbf{i}_\mathrm{w1} \in \mathbb{R}^N$
$\mathbf{i}_\mathrm{w2}$ | [p.u.] | Winding 2 current | Positive into the winding, $\mathbf{i}_\mathrm{w2} \in \mathbb{R}^N$

### External Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}_1$ | [V] | Terminal 1 voltage owned by EMT bus | $\mathbf{v}_1 \in \mathbb{R}^N$
$\mathbf{v}_2$ | [V] | Terminal 2 voltage owned by EMT bus | $\mathbf{v}_2 \in \mathbb{R}^N$

#### Algebraic

Bus voltages are algebraic when their current-balance rows have no voltage derivative.

## Model Equations

The winding voltages are

```math
\mathbf{v}_{\mathrm{w},e} = \dfrac{1}{V_{\mathrm{pk},e}}\,\mathbf{P}_e^\mathsf{T}\mathbf{v}_e,
\quad e \in \{1,2\}.
```

The magnetizing characteristic applies elementwise, using the
[ramp](../../../../CommonMath.md#ramp) $\rho$:

```math
g(\psi) = \dfrac{\psi}{L_\mathrm{m}}
  + \left(\dfrac{1}{L_\mathrm{sat}} - \dfrac{1}{L_\mathrm{m}}\right)
    \left(\rho(\psi - \psi_\mathrm{K}) - \rho(-\psi - \psi_\mathrm{K})\right).
```

### Internal Equations

#### Differential

```math
\begin{aligned}
0 &= X\dfrac{\mathrm{d}\mathbf{i}_{12}}{\mathrm{d}t}
     + \omega_\mathrm{b}(\mathbf{e}_2 - \mathbf{e}_1) \\
0 &= \dfrac{\mathrm{d}\boldsymbol{\psi}_1}{\mathrm{d}t}
     - \omega_\mathrm{b}\mathbf{e}_1 \\
0 &= \dfrac{\mathrm{d}\boldsymbol{\psi}_2}{\mathrm{d}t}
     - \omega_\mathrm{b}\mathbf{e}_2
\end{aligned}
```

#### Algebraic

```math
\begin{aligned}
0 &= \mathbf{e}_1 - \mathbf{v}_{\mathrm{w},1} + R_1\,\mathbf{i}_\mathrm{w1} \\
0 &= \mathbf{e}_2 - \tau\,\mathbf{v}_{\mathrm{w},2} + \tau R_2\,\mathbf{i}_\mathrm{w2} \\
0 &= \mathbf{i}_\mathrm{m1} - \beta_1\,g(\boldsymbol{\psi}_1) - \beta_1 G_\mathrm{c}\,\mathbf{e}_1 \\
0 &= \mathbf{i}_\mathrm{m2} - \beta_2\,g(\boldsymbol{\psi}_2) - \beta_2 G_\mathrm{c}\,\mathbf{e}_2 \\
0 &= \mathbf{i}_\mathrm{w1} - \mathbf{i}_\mathrm{m1} - \mathbf{i}_{12} \\
0 &= \mathbf{i}_\mathrm{w2} + \tau\,\mathbf{i}_{12} - \tau\,\mathbf{i}_\mathrm{m2}
\end{aligned}
```

### External Equations

```math
\begin{aligned}
\mathbf{i}_1 &\leftarrow -I_{\mathrm{pk},1}\,\mathbf{P}_1\,\mathbf{i}_\mathrm{w1} \\
\mathbf{i}_2 &\leftarrow -I_{\mathrm{pk},2}\,\mathbf{P}_2\,\mathbf{i}_\mathrm{w2}
\end{aligned}
```

## Initialization

For sinusoidal terminal voltages at angular frequency $\omega>0$, balanced
positive-sequence samples determine their derivatives:

```math
\dfrac{\mathrm{d}\mathbf{v}_e}{\mathrm{d}t}
=\dfrac{\omega}{\sqrt{3}}
\begin{bmatrix}v_{e,c}-v_{e,b}\\v_{e,a}-v_{e,c}\\v_{e,b}-v_{e,a}\end{bmatrix}
```

Otherwise, supply the voltage derivatives. With $w=\omega/\omega_\mathrm{b}$,

```math
\begin{aligned}
\bar{\mathbf{v}}_{\mathrm{w},e}
  &= \mathbf{v}_{\mathrm{w},e}
     -\dfrac{\mathrm{j}}{\omega}\dfrac{\mathrm{d}\mathbf{v}_{\mathrm{w},e}}{\mathrm{d}t} \\
\bar{y}_e &= \beta_e\left(G_\mathrm{c}-\dfrac{\mathrm{j}}{wL_\mathrm{m}}\right) \\
\bar{z} &= \mathrm{j}wX
\end{aligned}
```

Solve the linear magnetizing circuit:

```math
\begin{bmatrix}
1+R_1\bar{y}_1+R_1/\bar{z} & -R_1/\bar{z} \\
-\tau^2R_2/\bar{z} & 1+\tau^2R_2\bar{y}_2+\tau^2R_2/\bar{z}
\end{bmatrix}
\begin{bmatrix}\bar{\mathbf{e}}_1\\\bar{\mathbf{e}}_2\end{bmatrix}
=
\begin{bmatrix}\bar{\mathbf{v}}_{\mathrm{w},1}\\\tau\bar{\mathbf{v}}_{\mathrm{w},2}\end{bmatrix}
```

```math
\begin{aligned}
\mathbf{i}_{12} &\leftarrow
  \operatorname{Re}\left(\dfrac{\bar{\mathbf{e}}_1-\bar{\mathbf{e}}_2}{\bar{z}}\right) \\
\boldsymbol{\psi}_e &\leftarrow
  \operatorname{Re}\left(-\dfrac{\mathrm{j}}{w}\bar{\mathbf{e}}_e\right),
  \qquad e\in\{1,2\}
\end{aligned}
```

Evaluate algebraic variables and derivatives from the nonlinear model equations.
This initializes a consistent state; saturated periodic steady state requires a
periodic solve. Without the sinusoidal assumption, instantaneous voltages do not fix the initial states.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`i12` | [p.u.] | Series leakage current from node 1 to node 2 | $\mathbf{i}_{12} \in \mathbb{R}^N$
`psi1` | [p.u.] | Node 1 magnetizing flux linkage | $\boldsymbol{\psi}_1 \in \mathbb{R}^N$
`psi2` | [p.u.] | Node 2 magnetizing flux linkage | $\boldsymbol{\psi}_2 \in \mathbb{R}^N$
`i1` | [A] | Terminal 1 current injection | $\mathbf{i}_1 \in \mathbb{R}^N$
`i2` | [A] | Terminal 2 current injection | $\mathbf{i}_2 \in \mathbb{R}^N$
