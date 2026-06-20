# Bus Model

`Bus` represents a three-phase bus in instantaneous abc coordinates. The
bus voltages are differential variables, and the model equations enforce
three-phase current balance at the bus.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$N$ | [-] | `N` | Number of connected-device ports | Required, positive integer

### Parameter Validation

```math
N > 0
```

### Model Derived Parameters

None.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}$ | [V] | Bus voltage vector | $\mathbf{v} = [v_a, v_b, v_c]^T \in \mathbb{R}^3$

#### Algebraic

None.

### External Variables

#### Differential

None.

#### Algebraic

None.

## Model Ports

Symbol | Type | Units | Description | Note
------ | ---- | ----- | ----------- | ----
$\mathbf{i}^{\mathrm{inj}}_e$ | Input | [A] | Current injection from connected device $e$ | $e=1,\ldots,N$

## Model Equations

### Differential Equations

```math
\begin{aligned}
0 &= \sum_{e=1}^{N} \mathbf{i}^{\mathrm{inj}}_e
\end{aligned}
```

Each $\mathbf{i}^{\mathrm{inj}}_e$ may depend on the bus voltage and bus voltage derivative.

### Algebraic Equations

None.

## Initialization

For a balanced three-phase initialization derived from phasor voltage
$V = |V| \angle \phi$ and nominal angular frequency $\omega_0 = 2 \pi f_0$:

```math
\begin{aligned}
\mathbf{v}(0) = \sqrt{2}|V|
\begin{bmatrix}
  \cos(\phi) \\
  \cos(\phi - \tfrac{2\pi}{3}) \\
  \cos(\phi + \tfrac{2\pi}{3})
\end{bmatrix} \\
\dot{\mathbf{v}}(0) = -\sqrt{2}|V|\omega_0
\begin{bmatrix}
  \sin(\phi) \\
  \sin(\phi - \tfrac{2\pi}{3}) \\
  \sin(\phi + \tfrac{2\pi}{3})
\end{bmatrix}
\end{aligned}
```

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`v` | [V] | Bus voltage | $\mathbf{v} \in \mathbb{R}^3$
