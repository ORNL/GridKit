# LoadRL Model

`LoadRL` represents a three-phase RL load in instantaneous abc coordinates.
The load owns the three-phase differential current vector $\mathbf{i}$,
which is directed from the load into the bus.

## Model Parameters

Symbol  | Units      | Description                | Note
--------|------------|----------------------------|--------------------------------
$R_a$   | [$\Omega$] | Load resistance, phase a    |
$R_b$   | [$\Omega$] | Load resistance, phase b    |
$R_c$   | [$\Omega$] | Load resistance, phase c    |
$L_a$   | [H]        | Load inductance, phase a    |
$L_b$   | [H]        | Load inductance, phase b    |
$L_c$   | [H]        | Load inductance, phase c    |

## Model Derived Parameters

``` math
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
$\mathbf{i}$        | [A]   | Load current vector, directed from load into bus | $\mathbf{i} = [i_a, i_b, i_c]^T \in \mathbb{R}^3$

#### Algebraic

None.

### External Variables

External variables enter component model equations but are owned by
other components. The EMT bus at the load terminal owns the voltage
variable and provides the equation needed to have a balanced system
of equations.

#### Differential

Symbol           | Units | Description                                  | Note
-----------------|-------|----------------------------------------------|---------------------------------
$\mathbf{v}$     | [V]   | Terminal voltage vector, owned by EMT bus    | $\mathbf{v} = [v_a, v_b, v_c]^T \in \mathbb{R}^3$

#### Algebraic

None.

## Model Equations

### Differential Equations

``` math
0 = \mathbf{R}\,\mathbf{i} + \mathbf{L}\dot{\mathbf{i}} + \mathbf{v}
```

### Algebraic Equations

None.

### Bus Residual Contributions

The RL load contributes to the KCL residual at its terminal bus. The
expression is accumulated into the owning bus residual.

``` math
\mathbf{i}^\text{inj} := \mathbf{i}
```

## Initialization

The initialization assumes a balanced three-phase system. Given the power
flow phasor load current $I = |I| \angle \theta$, the initial load
current is:

``` math
\mathbf{i}(0) = \sqrt{2}\,|I|
\begin{bmatrix}
  \cos(\theta) \\
  \cos(\theta - \tfrac{2\pi}{3}) \\
  \cos(\theta + \tfrac{2\pi}{3})
\end{bmatrix}
```

The initial derivative is then given by the RL load equation for DAE
consistency:

``` math
\dot{\mathbf{i}}(0) = -\mathbf{L}^{-1}\left(\mathbf{v}(0) + \mathbf{R}\,\mathbf{i}(0)\right)
```

## Model Outputs

Candidate monitorable outputs include the load current components $i_a$, $i_b$, and $i_c$ (into the bus).
