# BranchLumpedConstant Model

`BranchLumpedConstant` represents a lumped-parameter EMT transmission line.
The nominal $\pi$-model is obtained by spatially discretizing the telegrapher equations over
a segment of length $\Delta x$, with a half shunt placed at each terminal.
Series current $\mathbf{i}$ is directed from bus 1 to bus 2. Bus residual
current injections are positive into buses. All electrical parameter matrices
are $3 \times 3$ and capture self and mutual coupling between phases.

<div align="center">
   <img align="center" src="../../../../../docs/Figures/EMT/lumped_constant_diagram.svg">

  Figure 1: Lumped constant EMT branch model
</div>

## Model Parameters

Symbol           | Units          | Description                              | Note
-----------------|----------------|------------------------------------------|---------------------------------
$\mathbf{R}'$    | [$\Omega$/m]   | Series resistance matrix per unit length | $\mathbb{R}^{3 \times 3}$
$\mathbf{L}'$    | [H/m]          | Series inductance matrix per unit length | $\mathbb{R}^{3 \times 3}$
$\mathbf{G}'$    | [S/m]          | Shunt conductance matrix per unit length | $\mathbb{R}^{3 \times 3}$
$\mathbf{C}'$    | [F/m]          | Shunt capacitance matrix per unit length | $\mathbb{R}^{3 \times 3}$
$\Delta x$       | [m]            | Line segment length                      | $\mathbb{R}$

## Model Derived Parameters

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
$\mathbf{i}$   | [A]    | Series branch current, directed bus 1 to bus 2 | $\mathbf{i} = [i_a, i_b, i_c]^T \in \mathbb{R}^3$

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
$\mathbf{v}_1$   | [V]    | Terminal voltage at bus 1, owned by bus 1 | $\mathbf{v}_1 = [v_{1,a}, v_{1,b}, v_{1,c}]^T \in \mathbb{R}^3$
$\mathbf{v}_2$   | [V]    | Terminal voltage at bus 2, owned by bus 2 | $\mathbf{v}_2 = [v_{2,a}, v_{2,b}, v_{2,c}]^T \in \mathbb{R}^3$

#### Algebraic

None.

## Model Equations

### Differential Equations


``` math
0 = \mathbf{R}\,\mathbf{i} + \mathbf{L}\dot{\mathbf{i}} + \mathbf{v}_2 - \mathbf{v}_1
```

### Algebraic Equations

None.

### Bus Residual Contributions

The lumped line contributes to the KCL residual at each terminal bus.
Each expression is accumulated into the owning bus residual.

``` math
\mathbf{i}^\text{inj}_1 := - \dfrac{\mathbf{G}}{2}\,\mathbf{v}_1 - \dfrac{\mathbf{C}}{2}\,\dot{\mathbf{v}}_1 - \mathbf{i}
```

``` math
\mathbf{i}^\text{inj}_2 := - \dfrac{\mathbf{G}}{2}\,\mathbf{v}_2 - \dfrac{\mathbf{C}}{2}\,\dot{\mathbf{v}}_2 + \mathbf{i}
```

## Initialization

The initialization assumes a balanced three-phase system. Given bus
voltages $\mathbf{v}_1(0)$, $\mathbf{v}_2(0)$ and their time
derivatives $\dot{\mathbf{v}}_1(0)$, $\dot{\mathbf{v}}_2(0)$ from
the EMT bus, and the power flow phasor series current
$I = |I| \angle \theta$, the initial series current is:

``` math
\mathbf{i}(0) = \sqrt{2}\,|I|
\begin{bmatrix}
  \cos(\theta) \\
  \cos(\theta - \tfrac{2\pi}{3}) \\
  \cos(\theta + \tfrac{2\pi}{3})
\end{bmatrix}
```

The initial derivative is then given by the series branch equation for
DAE consistency:

``` math
\dot{\mathbf{i}}(0) = \mathbf{L}^{-1}\left(\mathbf{v}_1(0) - \mathbf{v}_2(0) - \mathbf{R}\,\mathbf{i}(0)\right)
```

## Model Outputs

Candidate monitorable outputs include the series branch current components
$i_a$, $i_b$, and $i_c$.

Terminal current injection expressions are documented above as
$\mathbf{i}^\text{inj}_1$ and $\mathbf{i}^\text{inj}_2$.
