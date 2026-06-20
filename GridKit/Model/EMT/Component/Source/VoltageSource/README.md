# VoltageSource Model

`VoltageSource` represents a sinusoidal EMT voltage source in instantaneous
phase coordinates. The source voltage vector is connected to the EMT bus
through terminal conductance $\mathbf{g}_{\mathrm{s}}$.

## Model Parameters

For phase count $N$:

Symbol                         | Units   | JSON     | Description                       | Note
-------------------------------|---------|----------|-----------------------------------|-----
$\mathbf{E}$                   | [V]     | `E`      | Source voltage magnitudes         | $\mathbf{E}\in\mathbb{R}^N$, RMS
$\boldsymbol{\phi}$            | [rad]   | `phi`    | Source phase offsets              | $\boldsymbol{\phi}\in\mathbb{R}^N$
$\omega$                       | [rad/s] | `omega`  | Source angular frequency          |
$\mathbf{g}_{\mathrm{s}}$      | [S]     | `G`      | Terminal conductance              | $\mathbf{g}_{\mathrm{s}}\in\mathbb{R}^N$

### Parameter Validation

```math
\begin{aligned}
\mathbf{E} &\ge \mathbf{0} \\
\omega &> 0 \\
\mathbf{g}_{\mathrm{s}} &> \mathbf{0}
\end{aligned}
```

### Model Derived Parameters

The phase count $N$ is the length of $\mathbf{E}$.

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

None.

### External Variables

External variables are owned by the EMT bus.

#### Differential

Symbol           | Units | Description                                  | Note
-----------------|-------|----------------------------------------------|---------------------------------
$\mathbf{v}$     | [V]   | Bus voltage vector                          | $\mathbf{v} \in \mathbb{R}^N$

#### Algebraic

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{i}^{\mathrm{inj}}$ | `i` | Output | [A] | Current injection at source port | $\mathbf{i}^{\mathrm{inj}} \in \mathbb{R}^N$

## Model Equations

### Differential Equations

None.

### Algebraic Equations

None.

### Wiring

```math
i^{\mathrm{inj}}_n
=
g_{\mathrm{s},n}
\left(
\sqrt{2}E_n\cos\left(\omega t + \phi_n\right) - v_n
\right)
```

## Initialization

No internal state is initialized. Initial values satisfy the wiring equations:

```math
i^{\mathrm{inj}}_{n,0}
=
g_{\mathrm{s},n}
\left(
\sqrt{2}E_n\cos\left(\phi_n\right) - v_{n,0}
\right)
```

## Monitors

None.
