# LoadZ Model

`LoadZ` represents an $N$-phase impedance load. Current $\mathbf{i}$ is injected
from the load into the EMT bus.

> [!NOTE]
> The initial end-to-end implementation will support three-phase systems only
> to establish a proof of concept. The formulation below remains $N$-phase.

## Block Diagram

None.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$N$ | [-] | `N` | Number of phases | Required, positive integer

### Parameter Validation

```math
N > 0
```

### Derived Parameters

None.

## Submodels

Submodel | Type | Order | Inputs | Parameters | Outputs
-------- | ---- | ----- | ------ | ---------- | -------
Impedance $f^{\mathbf{z}}$ | [`VectorFit`](../../../Operators/Rational/VectorFit/README.md) | $Q_{\mathbf{z}}$ | $\mathbf{i} \in \mathbb{R}^N$ | `Z` | $f^{\mathbf{z}}(\mathbf{i}) \in \mathbb{R}^N$

### Submodel Validation

```math
\begin{aligned}
f^{\mathbf{z}} &: \mathbb{R}^N \rightarrow \mathbb{R}^N \\
\mathbf{E}^{\mathbf{z}} &\in \mathbb{R}^{N \times N} \\
\mathrm{rank}\left(\mathbf{E}^{\mathbf{z}}\right) &= N
\end{aligned}
```

$\mathbf{E}^{\mathbf{z}}$ is the linear coefficient of the impedance fit. The
full-rank condition is required by this formulation's differential-current
classification.
Static or singular fits require a corresponding algebraic-current formulation.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}$ | [A] | Current injection from load into EMT bus | $\mathbf{i} \in \mathbb{R}^N$

#### Algebraic

None.

### External Variables

External variables are owned by the EMT bus.

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}$ | [V] | Port voltage vector owned by EMT bus | $\mathbf{v} \in \mathbb{R}^N$

#### Algebraic

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{i}^{\mathrm{inj}}$ | `i` | Output | [A] | Current injection at load port | $\mathbf{i}^{\mathrm{inj}} \in \mathbb{R}^N$

## Model Equations

### Differential Equations

```math
0 = f^{\mathbf{z}}(\mathbf{i}) + \mathbf{v}
```

### Algebraic Equations

None.

### Wiring

```math
\mathbf{i}^{\mathrm{inj}} = \mathbf{i}
```

## Initialization

### Input Initialization

```math
\begin{aligned}
\mathbf{v}
  &\leftarrow \text{initialized bus voltage}
\end{aligned}
```

### Internal Initialization

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$\mathbf{i}_0$ | [A] | `init.i0` | Initial current-injection vector | Optional, defaults to $\mathbf{0}$

The current and a provisional derivative seed the impedance submodel before the
system-level consistent-initial-condition calculation:

```math
\begin{aligned}
\mathbf{i}
  &\leftarrow \mathbf{i}_0 \\
\dfrac{\mathrm{d}\mathbf{i}}{\mathrm{d}t}
  &\leftarrow \mathbf{0}
\end{aligned}
```

The assembled load and network residuals then replace the provisional seed with
consistent differential derivatives and algebraic outputs while enforcing

```math
f^{\mathbf{z}}(\mathbf{i}) + \mathbf{v} \leftarrow \mathbf{0}
```

A network operating point should provide a consistent $\mathbf{i}_0$ when a
zero-current seed is not appropriate.

### Output Initialization

```math
\mathbf{i}^{\mathrm{inj}} \leftarrow \mathbf{i}
```

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`i` | [A] | Load current injection | $\mathbf{i} \in \mathbb{R}^N$
