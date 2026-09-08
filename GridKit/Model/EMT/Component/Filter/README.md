# Filter Model

`Filter` represents a three-phase LCL filter between a converter and a terminal
Bus. It owns the converter-side current $\mathbf{i}$, capacitor voltage
$\mathbf{v}_{\mathrm{o}}$, and grid-side current $\mathbf{i}_{\mathrm{g}}$.
Both currents are positive from the converter toward the Bus.

## Block Diagram

![Filter model block diagram](../../../../../docs/Figures/EMT/Filter/diagram.png)

Figure 1: Filter model

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$\mathbf{R}_{\mathrm{s}}$ | [Ω] | `Rs` | Converter-side series resistance | $\mathbb{R}^{3\times3}$, default zero
$\mathbf{L}_{\mathrm{s}}$ | [H] | `Ls` | Converter-side series inductance | $\mathbb{R}^{3\times3}$, required
$\mathbf{C}$ | [F] | `C` | Shunt capacitance | $\mathbb{R}^{3\times3}$, required
$\mathbf{R}_{\mathrm{g}}$ | [Ω] | `Rg` | Grid-side series resistance | $\mathbb{R}^{3\times3}$, default zero
$\mathbf{L}_{\mathrm{g}}$ | [H] | `Lg` | Grid-side series inductance | $\mathbb{R}^{3\times3}$, required

### Parameter Validation

All matrices must be finite and symmetric. The resistance matrices must be
positive semidefinite; the inductance and capacitance matrices must be positive
definite. Off-diagonal entries represent coupling between phases.

### Derived Parameters

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{v}$ | `v` | Input | [V] | Terminal Bus voltage | $\mathbb{R}^3$
$\mathbf{e}$ | `e` | Input | [V] | Converter output voltage | $\mathbb{R}^3$
$\mathbf{i}$ | `i` | Output | [A] | Converter-side current | $\mathbb{R}^3$
$\mathbf{v}_{\mathrm{o}}$ | `vo` | Output | [V] | Capacitor voltage | $\mathbb{R}^3$
$\mathbf{i}_{\mathrm{g}}$ | `ig` | Output | [A] | Current injected into the terminal Bus | $\mathbb{R}^3$

Vector ports use phase order `a`, `b`, `c`. Scalar names append the phase letter,
for example `ea`, `voa`, and `iga`. The `bus` input shortcut binds `v` to the
terminal Bus. Output aliases are optional; the Bus receives the grid-side current
without an explicit signal alias.

Connect `Converter.vo` to `e` and return `i` to `Converter.i`. Controllers and
reference-frame operators read `i`, `vo`, and `ig` through signal ports. A PLL
can read `vo` and supply angle and frequency to the controllers.

## Submodels

None.

### Submodel Validation

None.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}$ | [A] | Converter-side inductor current | $\mathbb{R}^3$
$\mathbf{v}_{\mathrm{o}}$ | [V] | Capacitor voltage | $\mathbb{R}^3$
$\mathbf{i}_{\mathrm{g}}$ | [A] | Grid-side inductor current | $\mathbb{R}^3$

#### Algebraic

None.

### External Variables

#### Differential

The connected voltage signals may depend on differential variables owned by
other components. The Filter does not use their derivatives.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}$ | [V] | Terminal voltage owned by the Bus | $\mathbb{R}^3$, when algebraic
$\mathbf{e}$ | [V] | Converter output voltage | $\mathbb{R}^3$, may be computed from other signals

## Model Equations

### Internal Equations

#### Differential

```math
\begin{aligned}
0 &= \mathbf{L}_{\mathrm{s}}\dfrac{\mathrm{d}\mathbf{i}}{\mathrm{d}t}
   + \mathbf{R}_{\mathrm{s}}\mathbf{i} + \mathbf{v}_{\mathrm{o}}-\mathbf{e} \\
0 &= \mathbf{C}\dfrac{\mathrm{d}\mathbf{v}_{\mathrm{o}}}{\mathrm{d}t}
   + \mathbf{i}_{\mathrm{g}}-\mathbf{i} \\
0 &= \mathbf{L}_{\mathrm{g}}\dfrac{\mathrm{d}\mathbf{i}_{\mathrm{g}}}{\mathrm{d}t}
   + \mathbf{R}_{\mathrm{g}}\mathbf{i}_{\mathrm{g}} + \mathbf{v}-\mathbf{v}_{\mathrm{o}}
\end{aligned}
```

#### Algebraic

None.

### External Equations

None. The terminal Bus registers $\mathbf{i}_{\mathrm{g}}$ as a current
injection and owns its KCL equations.

## Initialization

The state file may prescribe `ia`, `ib`, `ic`, `voa`, `vob`, `voc`, `iga`,
`igb`, and `igc`; omitted values default to zero. The Filter publishes these
initial output values to downstream components, including a connected PLL.
Consistent initialization preserves all nine differential states and resolves
their derivatives under the EMT initialization contract.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`i` | [A] | Converter-side current | Expands to `ia`, `ib`, `ic`
`vo` | [V] | Capacitor voltage | Expands to `voa`, `vob`, `voc`
`ig` | [A] | Grid-side current | Expands to `iga`, `igb`, `igc`
