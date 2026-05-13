# BranchEMT Model

`BranchEMT` represents a lumped-parameter EMT transmission line. The nominal
pi model is obtained by spatially discretizing the telegrapher equations over
a segment of length $\Delta x$, with a half shunt placed at each terminal.
Series current $\mathbf{i}_s$ is directed from bus 1 to bus 2. The positive
flow direction is into buses. All parameters are $3 \times 3$ matrices
capturing self and mutual coupling between phases.

## Model Parameters

Symbol           | Units          | Description                              | Note
-----------------|----------------|------------------------------------------|---------------------------------
$\mathbf{R}'$    | [$\Omega$/m]   | Series resistance matrix per unit length | $\mathbb{R}^{3 \times 3}$
$\mathbf{L}'$    | [H/m]          | Series inductance matrix per unit length | $\mathbb{R}^{3 \times 3}$
$\mathbf{G}'$    | [S/m]          | Shunt conductance matrix per unit length | $\mathbb{R}^{3 \times 3}$
$\mathbf{C}'$    | [F/m]          | Shunt capacitance matrix per unit length | $\mathbb{R}^{3 \times 3}$
$\Delta x$       | [m]            | Line segment length                      | $\mathbb{R}$

### Model Derived Parameters

``` math
\begin{aligned}
  \mathbf{R} &= \mathbf{R}'\Delta x   & \mathbf{G} &= \mathbf{G}'\Delta x \\
  \mathbf{L} &= \mathbf{L}'\Delta x   & \mathbf{C} &= \mathbf{C}'\Delta x
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

Symbol           | Units  | Description           | Note
-----------------|--------|-----------------------|---------------------------------
$\mathbf{i}_s$   | [A]    | Series branch current | $\mathbb{R}^3$, directed bus 1 to bus 2

#### Algebraic

None.

### External Variables

External variables enter component model equations but are owned by
other components. The EMT bus at each terminal owns the voltage
variable and provides the equation needed to have a balanced system
of equations.

#### Differential

Symbol           | Units  | Description              | Note
-----------------|--------|--------------------------|------------------
$\mathbf{v}_1$   | [V]    | Terminal voltage at bus 1 | $\mathbb{R}^3$, owned by EMT bus
$\mathbf{v}_2$   | [V]    | Terminal voltage at bus 2 | $\mathbb{R}^3$, owned by EMT bus

#### Algebraic

None.

## Model Equations

### Differential Equations

``` math
\dot{\mathbf{i}}_s = \mathbf{L}^{-1}\left((\mathbf{v}_1 - \mathbf{v}_2) - \mathbf{R}\,\mathbf{i}_s\right)
```

(or, alternatively)

``` math
0 = (\mathbf{v}_1 - \mathbf{v}_2) - \mathbf{R}\,\mathbf{i}_s - \mathbf{L}\dot{\mathbf{i}}_s
```

### Algebraic Equations

None.

### Bus Residual Contributions

The lumped line contributes to the KCL residual at each terminal bus.
Each expression is accumulated into the owning bus residual.

``` math
\Delta \mathbf{i}_1 = -\mathbf{i}_s - \dfrac{\mathbf{C}}{2}\,\dot{\mathbf{v}}_1 - \dfrac{\mathbf{G}}{2}\,\mathbf{v}_1
```

``` math
\Delta \mathbf{i}_2 = +\mathbf{i}_s - \dfrac{\mathbf{C}}{2}\,\dot{\mathbf{v}}_2 - \dfrac{\mathbf{G}}{2}\,\mathbf{v}_2
```

## Initialization

The initialization assumes a balanced three-phase system. Given bus
voltages $\mathbf{v}_1(0)$, $\mathbf{v}_2(0)$ and their time
derivatives $\dot{\mathbf{v}}_1(0)$, $\dot{\mathbf{v}}_2(0)$ from
the EMT bus, and the power flow phasor series current
$I_s = |I_s| \angle \theta$, the initial series current is:

``` math
\mathbf{i}_s(0) = \sqrt{2}\,|I_s|
\begin{bmatrix}
  \cos(\theta) \\
  \cos(\theta - \tfrac{2\pi}{3}) \\
  \cos(\theta + \tfrac{2\pi}{3})
\end{bmatrix}
```

The initial derivative is then given by the series branch equation for
DAE consistency:

``` math
\dot{\mathbf{i}}_s(0) = \mathbf{L}^{-1}\left((\mathbf{v}_1(0) - \mathbf{v}_2(0)) - \mathbf{R}\,\mathbf{i}_s(0)\right)
```
