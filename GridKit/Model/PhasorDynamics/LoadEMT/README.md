# LoadEMT Model

`LoadEMT` represents a three-phase RL load in instantaneous abc coordinates.
The load owns three differential current variables, one for each phase. Load
current $\mathbf{i}$ is directed from the load into the bus.

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

Symbol              | Units | Description              | Note
--------------------|-------|--------------------------|---------------------------------
$i_a$               | [A]   | Load current, phase a     | directed from load into bus
$i_b$               | [A]   | Load current, phase b     | directed from load into bus
$i_c$               | [A]   | Load current, phase c     | directed from load into bus

#### Algebraic

None.

### External Variables

External variables enter component model equations but are owned by
other components. The EMT bus at the load terminal owns the voltage
variable and provides the equation needed to have a balanced system
of equations.

#### Differential

Symbol  | Units | Description             | Note
--------|-------|-------------------------|------------------
$v_a$   | [V]   | Terminal voltage, phase a | owned by EMT bus
$v_b$   | [V]   | Terminal voltage, phase b | owned by EMT bus
$v_c$   | [V]   | Terminal voltage, phase c | owned by EMT bus

#### Algebraic

None.

## Model Equations

### Differential Equations

``` math
\dot{\mathbf{i}} = -\mathbf{L}^{-1}\left(\mathbf{v} + \mathbf{R}\,\mathbf{i}\right)
```

(or, alternatively)

``` math
0 = \mathbf{v} + \mathbf{R}\,\mathbf{i} + \mathbf{L}\dot{\mathbf{i}}
```

### Algebraic Equations

None.

### Bus Residual Contributions

The RL load contributes to the KCL residual at its terminal bus. The
expression is accumulated into the owning bus residual.

``` math
\Delta \mathbf{i} = \mathbf{i}
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
