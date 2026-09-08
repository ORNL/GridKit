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
$v_{\mathrm{dc}}$ | `vdc` | Input | [V] | DC-link voltage | $v_{\mathrm{dc}} \ge 0$
$\mathbf{i}$ | `i` | Input | [A] | AC terminal currents | $\mathbf{i} \in \mathbb{R}^3$, positive out of the bridge
$\mathbf{e}$ | `e` | Output | [V] | Bridge voltage vector | $\mathbf{e} \in \mathbb{R}^3$
$i_{\mathrm{dc}}$ | `idc` | Output | [A] | DC-link current | Positive into the bridge

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

The connected DC-link voltage and AC-current variables may be differential.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{s}$ | [-] | Switching function vector | $\mathbf{s} \in [0,1]^3$
$v_{\mathrm{dc}}$ | [V] | DC-link voltage | $v_{\mathrm{dc}} \ge 0$
$\mathbf{i}$ | [A] | AC terminal currents | $\mathbf{i} \in \mathbb{R}^3$, positive out of the bridge

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

None.

### External Equations

```math
\begin{aligned}
\mathbf{e} &\leftarrow v_{\mathrm{dc}}\mathbf{P}\mathbf{s} \\
i_{\mathrm{dc}} &\leftarrow (\mathbf{P}\mathbf{s})^\mathsf{T}\mathbf{i}
\end{aligned}
```

The current transformation preserves instantaneous power, including at zero
DC voltage, without division by $v_{\mathrm{dc}}$:

```math
v_{\mathrm{dc}} i_{\mathrm{dc}}
  = \mathbf{e}^\mathsf{T}\mathbf{i}.
```

Positive $i_{\mathrm{dc}}$ flows physically into the bridge, but its signal is
an output: the switching function and AC currents determine the DC current
drawn. [DCLink](../../Component/Controller/DCLink/README.md) receives this current
and supplies $v_{\mathrm{dc}}$. Its source-current input $i_{\mathrm{src}}$ is
separate from the converter current.

## Initialization

None beyond the EMT initialization contract.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`e` | [V] | Bridge voltage | $\mathbf{e} \in \mathbb{R}^3$
`idc` | [A] | DC-link current | Positive into the bridge

In case JSON, `mon: ["e"]` expands to the scalar monitors `ea`, `eb`, `ec`.
See [case connections](../../INPUT_FORMAT.md#case-connections) for vector signal wiring.
