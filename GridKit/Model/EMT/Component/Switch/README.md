# Switch Model

`Switch` represents an ideal $N$-phase EMT switch between two buses. Series
current $\mathbf{i}_{12}$ is directed from terminal 1 to terminal 2. The
Boolean `open` operates all phases: `true` is open and `false` is closed. The
switch contains no energy storage; switching transients arise from the
connected EMT network.

## Block Diagram

![Switch model block diagram](../../../../../docs/Figures/EMT/Switch/diagram.png)

Figure 1: Switch model

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$N$ | [-] | `N` | Number of phases | Required, positive integer
$\mathrm{open}$ | [-] | `open` | Initial switch command | Optional Boolean, `true` open, `false` closed (default)

### Parameter Validation

```math
N \in \mathbb{Z}_{>0}
```

### Model Derived Parameters

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{v}_1$, $\mathbf{i}_1$ | `bus1` | Bus | [V], [A] | Terminal 1 bus voltage and current injection | $\mathbf{v}_1, \mathbf{i}_1 \in \mathbb{R}^N$
$\mathbf{v}_2$, $\mathbf{i}_2$ | `bus2` | Bus | [V], [A] | Terminal 2 bus voltage and current injection | $\mathbf{v}_2, \mathbf{i}_2 \in \mathbb{R}^N$

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

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}_{12}$ | [A] | Series current from terminal 1 to terminal 2 | $\mathbf{i}_{12} \in \mathbb{R}^N$

### External Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}_1$ | [V] | Terminal 1 voltage owned by EMT bus | $\mathbf{v}_1 \in \mathbb{R}^N$
$\mathbf{v}_2$ | [V] | Terminal 2 voltage owned by EMT bus | $\mathbf{v}_2 \in \mathbb{R}^N$

#### Algebraic

None.

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

```math
\begin{cases}
\mathbf{i}_{12} = \mathbf{0}, & \text{open}\\
\mathbf{v}_2-\mathbf{v}_1 = \mathbf{0}, & \text{closed}
\end{cases}
```

### External Equations

```math
\begin{aligned}
\mathbf{i}_1 &\leftarrow -\mathbf{i}_{12} \\
\mathbf{i}_2 &\leftarrow \mathbf{i}_{12}
\end{aligned}
```

## Initialization

The `open` parameter sets the command at $t_0$, before the algebraic
equations are enforced. The switch has no additional model-specific
initialization.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`i12` | [A] | Series current from terminal 1 to terminal 2 | $\mathbf{i}_{12} \in \mathbb{R}^N$
