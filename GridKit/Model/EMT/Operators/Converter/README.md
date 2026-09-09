# Converter Model

`Converter` maps a DC voltage and three-phase switching function to the
bridge voltage of a two-level voltage-source inverter. The operator adds no
DAE variables or residual rows.

## Block Diagram

![Converter model block diagram](../../../../../docs/Figures/EMT/Converter/diagram.png)

Figure 1: Converter model

## Model Parameters

None.

### Parameter Validation

None.

### Derived Parameters

The normalized phase incidence matrix and zero-sequence removal projector are

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
$v_{\mathrm{dc}}$ | `vdc` | Input | [V] | DC voltage | $v_{\mathrm{dc}} \ge 0$
$\mathbf{e}$ | `e` | Output | [V] | Bridge voltage vector | $\mathbf{e} \in \mathbb{R}^3$

Vectors use $(a,b,c)$ order. The bridge voltage $\mathbf{e}$ is referred to
the AC neutral; $\mathbf{P}$ removes the common-mode pole voltage, so
$e_a+e_b+e_c=0$.

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

The connected DC voltage may be differential.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{s}$ | [-] | Switching function vector | $\mathbf{s} \in [0,1]^3$
$v_{\mathrm{dc}}$ | [V] | DC voltage | $v_{\mathrm{dc}} \ge 0$

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

None.

### External Equations

```math
\mathbf{e} \leftarrow v_{\mathrm{dc}}\mathbf{P}\mathbf{s}
```

## Initialization

[Balanced initialization](../../STATE.md#application) requires zero-sum bridge
voltage and positive DC voltage. It requests the carrier means
$\bar{s}_p=1/2+e_p/v_{\mathrm{dc}}$ from PWM. Otherwise, prescribed outputs must
match the evaluated bridge voltage.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`e` | [V] | Bridge voltage | $\mathbf{e} \in \mathbb{R}^3$

In case JSON, `mon: ["e"]` expands to the scalar monitors `ea`, `eb`, `ec`.
See [case connections](../../INPUT_FORMAT.md#case-connections) for vector signal wiring.
