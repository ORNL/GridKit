# Converter Model

`Converter` maps a DC-link voltage and three-phase switching function to the
bridge voltage of a two-level voltage-source inverter. AC terminal currents
determine the current drawn from the DC link. The lossless operator adds no DAE
variables or residual rows.

## Block Diagram

![Converter model block diagram](../../../../../docs/Figures/EMT/Converter/diagram.png)

Figure 1: Converter model

## Model Parameters

None.

### Parameter Validation

None.

### Derived Parameters

The normalized phase incidence matrix and zero-sequence projector are

```math
\begin{aligned}
\mathbf{A} &=
\dfrac{1}{\sqrt{3}}\begin{bmatrix}
1 & -1 & 0 \\
0 & 1 & -1 \\
-1 & 0 & 1
\end{bmatrix} \\
\mathbf{P} &= \mathbf{A}^\mathsf{T}\mathbf{A}
=
\dfrac{1}{3}
\begin{bmatrix}
2 & -1 & -1 \\
-1 & 2 & -1 \\
-1 & -1 & 2
\end{bmatrix}
\end{aligned}
```

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{s}$ | `s` | Input | [-] | Switching function vector | $\mathbf{s} \in [0,1]^3$
$V_{\mathrm{dc}}$ | `vdc` | Input | [V] | DC-link voltage | $V_{\mathrm{dc}} \ge 0$
$\mathbf{i}$ | `i` | Input | [A] | AC terminal currents | Positive out of the bridge
$\mathbf{v}_{\mathrm{o}}$ | `vo` | Output | [V] | Bridge voltage vector | $\mathbf{v}_{\mathrm{o}} \in \mathbb{R}^3$
$I_{\mathrm{dc}}$ | `idc` | Output | [A] | DC-link current | Positive into the bridge

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

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{s}$ | [-] | Switching function vector | $\mathbf{s} \in [0,1]^3$
$V_{\mathrm{dc}}$ | [V] | DC-link voltage | $V_{\mathrm{dc}} \ge 0$
$\mathbf{i}$ | [A] | AC terminal currents | Positive out of the bridge

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

None.

### External Equations

```math
\begin{aligned}
\mathbf{v}_{\mathrm{o}} &\leftarrow V_\mathrm{dc}\mathbf{P}\mathbf{s} \\
I_\mathrm{dc} &\leftarrow (\mathbf{P}\mathbf{s})^\mathsf{T}\mathbf{i}
\end{aligned}
```

The current transformation preserves instantaneous power, including at zero
DC voltage, without division by $V_\mathrm{dc}$:

```math
V_\mathrm{dc} I_\mathrm{dc}
  = \mathbf{v}_\mathrm{o}^\mathsf{T}\mathbf{i}.
```

## Initialization

None beyond the EMT initialization contract.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`vo` | [V] | Bridge voltage | $\mathbf{v}_{\mathrm{o}} \in \mathbb{R}^3$
`idc` | [A] | DC-link current | Positive into the bridge

In case JSON, `mon: ["vo"]` expands to the scalar monitors `voa`, `vob`, `voc`.
See [case connections](../../INPUT_FORMAT.md#case-connections) for vector signal wiring.
