# Switch Model

`Switch` represents an ideal three-phase EMT switch between two buses. Series
current $\mathbf{i}_{12}$ is directed from terminal 1 to terminal 2. The
Boolean command `open` operates all phases: `true` is open and
`false` is closed. The switch contains no energy storage; switching
transients arise from the connected EMT network.

## Block Diagram

![Switch model block diagram](../../../../../docs/Figures/EMT/Switch/diagram.png)

Figure 1: Switch model

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$I$ | [A] | `I` | Nominal phase RMS current | Optional, positive; absolute-tolerance scale
$N$ | [-] | `N` | Number of phases | Fixed at $3$
$\mathrm{open}$ | [-] | `open` | Ganged switch command | Default `false` (closed)

### Parameter Validation

```math
N = 3
```

### Derived Parameters

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{v}_1$ | `v1` | Input | [V] | Terminal 1 bus voltage | $\mathbf{v}_1 \in \mathbb{R}^N$
$\mathbf{v}_2$ | `v2` | Input | [V] | Terminal 2 bus voltage | $\mathbf{v}_2 \in \mathbb{R}^N$
$\mathbf{i}_{12}$ | `i12` | Output | [A] | Series current from terminal 1 to terminal 2 | $\mathbf{i}_{12} \in \mathbb{R}^N$

## Submodels

None.

### Submodel Validation

None.

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}_{12}$ | [A] | Series current from terminal 1 to terminal 2 | $\mathbf{i}_{12} \in \mathbb{R}^N$

### External Variables

#### Differential

When a connected equation depends on the bus-voltage derivative:

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}_1$ | [V] | Terminal 1 voltage owned by EMT bus | $\mathbf{v}_1 \in \mathbb{R}^N$
$\mathbf{v}_2$ | [V] | Terminal 2 voltage owned by EMT bus | $\mathbf{v}_2 \in \mathbb{R}^N$

#### Algebraic

Otherwise, the bus-voltage variables above are algebraic.

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

```math
\mathbf{0} =
\begin{cases}
\mathbf{i}_{12}, & \mathrm{open}=\mathrm{true} \\
\mathbf{v}_2-\mathbf{v}_1, & \mathrm{open}=\mathrm{false}
\end{cases}
```

### Terminal Currents

```math
\begin{aligned}
\mathbf{i}_1 &= -\mathbf{i}_{12} \\
\mathbf{i}_2 &= \mathbf{i}_{12}.
\end{aligned}
```

The bus registers these current signals and owns their residual and Jacobian
contributions to KCL.

## Initialization

The state-file `open` value overrides the parameter at $t_0$. Initial
`i12a`, `i12b`, and `i12c` values default to zero. The assembled consistent
initialization enforces the switch equations. Events update the command
through `setOpen`; it is not a signal input.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`open` | [-] | Switch command | Boolean; `true` open, `false` closed
`i12` | [A] | Series current from terminal 1 to terminal 2 | $\mathbf{i}_{12} \in \mathbb{R}^N$
