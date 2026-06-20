# LoadRL Model

`LoadRL` represents a three-phase RL load in instantaneous abc coordinates.
The load owns the three-phase differential current vector $\mathbf{i}$,
which is the current injection from the load into the EMT bus.

## Model Parameters

Symbol  | Units      | JSON | Description             | Note
--------|------------|------|-------------------------|-----
$R_a$   | [$\Omega$] | `Ra` | Load resistance, phase a |
$R_b$   | [$\Omega$] | `Rb` | Load resistance, phase b |
$R_c$   | [$\Omega$] | `Rc` | Load resistance, phase c |
$L_a$   | [H]        | `La` | Load inductance, phase a |
$L_b$   | [H]        | `Lb` | Load inductance, phase b |
$L_c$   | [H]        | `Lc` | Load inductance, phase c |

### Parameter Validation

```math
\begin{aligned}
R_a, R_b, R_c &\ge 0 \\
L_a, L_b, L_c &> 0
\end{aligned}
```

### Model Derived Parameters

```math
\begin{aligned}
  \mathbf{R} &= \operatorname{diag}(R_a, R_b, R_c) \\
  \mathbf{L} &= \operatorname{diag}(L_a, L_b, L_c)
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

Symbol              | Units | Description                                    | Note
--------------------|-------|------------------------------------------------|---------------------------------
$\mathbf{i}$        | [A]   | Current injection from load into EMT bus       | $\mathbf{i} = [i_a, i_b, i_c]^T \in \mathbb{R}^3$

#### Algebraic

None.

### External Variables

External variables are owned by the EMT bus.

#### Differential

Symbol           | Units | Description                                  | Note
-----------------|-------|----------------------------------------------|---------------------------------
$\mathbf{v}$     | [V]   | Port voltage vector                         | $\mathbf{v} = [v_a, v_b, v_c]^T \in \mathbb{R}^3$

#### Algebraic

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{i}^{\mathrm{inj}}$ | `i` | Output | [A] | Current injection at load port | $\mathbf{i}^{\mathrm{inj}} \in \mathbb{R}^3$

## Model Equations

### Differential Equations

```math
0 = \mathbf{R}\,\mathbf{i} + \mathbf{L}\dot{\mathbf{i}} + \mathbf{v}
```

### Algebraic Equations

None.

### Wiring

```math
\mathbf{i}^{\mathrm{inj}} = \mathbf{i}
```

## Initialization

For a balanced three-phase initialization derived from phasor current injection
$I^{\mathrm{inj}} = |I^{\mathrm{inj}}| \angle \theta$:

```math
\mathbf{i}_0 = \sqrt{2}\,|I^{\mathrm{inj}}|
\begin{bmatrix}
  \cos(\theta) \\
  \cos(\theta - \tfrac{2\pi}{3}) \\
  \cos(\theta + \tfrac{2\pi}{3})
\end{bmatrix}
```

Only $\mathbf{i}_0$ is initialized. The solver computes
$\dot{\mathbf{i}}_0$ from the differential residual.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`i` | [A] | Load current injection | $\mathbf{i} \in \mathbb{R}^3$
