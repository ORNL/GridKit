# Bridge Model

`Bridge` represents an ideal two-level converter with $N$ legs and two DC
terminals. Switching functions set the pole voltages and DC current injections.
The bridge is lossless and contains no energy storage.

## Block Diagram

![Bridge model diagram](../../../../../../docs/Figures/EMT/Bridge/diagram.png)

Figure 1: Representative three-phase bridge with binary contact positions;
the equations apply to $N$ legs and define continuous interpolation

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$N$ | [-] | `N` | Number of legs | Required, positive integer

### Parameter Validation

```math
N \in \mathbb{Z}_{>0}
```

### Model Derived Parameters

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{s}$ | `s` | Input | [-] | Switching functions | $\mathbf{s} \in \mathbb{R}^N$
$\mathbf{v}$, $\mathbf{i}$ | `ac` | Bus | [V], [A] | AC pole voltages and current injections | $\mathbf{v} \in \mathbb{R}^N$ and $\mathbf{i} \in \mathbb{R}^N$
$v_+$, $i_+$ | `vpos` | Bus | [V], [A] | Positive DC terminal voltage and current injection | Scalar
$v_-$, $i_-$ | `vneg` | Bus | [V], [A] | Negative DC terminal voltage and current injection | Scalar

All terminal currents are positive from the bridge into the connected buses.
Each switching function satisfies $0 \le s_n \le 1$ and may come from
[PWM](../../../Operators/Modulation/PWM/README.md).

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
$\mathbf{i}$ | [A] | AC current injections | $\mathbf{i} \in \mathbb{R}^N$

### External Variables

Classification follows the assembled DAE. The tables apply only to connected
DAE variables; prescribed and computed switching signals introduce no variables.

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}$ | [V] | AC pole voltages | Differential voltage, $\mathbf{v} \in \mathbb{R}^N$
$v_+$ | [V] | Positive DC terminal voltage | Differential voltage
$v_-$ | [V] | Negative DC terminal voltage | Differential voltage
$\mathbf{s}$ | [-] | Switching functions | Differential input, $\mathbf{s} \in \mathbb{R}^N$

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}$ | [V] | AC pole voltages | Algebraic voltage, $\mathbf{v} \in \mathbb{R}^N$
$v_+$ | [V] | Positive DC terminal voltage | Algebraic voltage
$v_-$ | [V] | Negative DC terminal voltage | Algebraic voltage
$\mathbf{s}$ | [-] | Switching functions | Algebraic input, $\mathbf{s} \in \mathbb{R}^N$

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

For each leg,

```math
0 = -v_n + (1-s_n)v_- + s_n v_+
```

### External Equations

The owned currents $\mathbf{i}$ contribute directly to the AC bus balances.
The DC current injections are

```math
\begin{aligned}
i_+ &\leftarrow -\sum_{n=1}^N s_n i_n \\
i_- &\leftarrow -\sum_{n=1}^N (1-s_n)i_n
\end{aligned}
```

## Initialization

None.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`i` | [A] | AC current injections | $\mathbf{i} \in \mathbb{R}^N$
`ipos` | [A] | Positive DC current injection | Scalar
`ineg` | [A] | Negative DC current injection | Scalar
`e` | [V] | Pole voltages relative to their mean | $\mathbf{e} \in \mathbb{R}^N$

## Development

Current and instantaneous power are conserved:

```math
\begin{aligned}
\sum_{n=1}^N i_n + i_+ + i_- &= 0 \\
\sum_{n=1}^N v_n i_n + v_+ i_+ + v_- i_- &= 0
\end{aligned}
```

Binary switching selects one DC rail per leg. Fractional switching functions
give a lossless interpolation; they do not represent partial semiconductor
conduction. Losses, dead time, and diode blocking are outside this model.

The derived voltage monitor is

```math
e_n \leftarrow v_n - \dfrac{1}{N}\sum_{k=1}^N v_k
```

It preserves line-to-line voltages without defining a load neutral or imposing
zero-sum AC current. Neutral and grounding connections belong to the network.
