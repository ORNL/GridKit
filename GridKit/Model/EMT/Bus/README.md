# Bus Model

`Bus` represents a three-phase bus in instantaneous abc coordinates. The
bus voltages are differential variables, and the model equations enforce
three-phase current balance at the bus.

## Model Parameters

None.

## Model Derived Parameters

None.

## Model Variables

### Internal Variables

#### Differential

Symbol   | Units | Description        | Note
---------|-------|--------------------|--------------------------------
$\mathbf{v}$ | [V] | Bus voltage vector | $\mathbf{v} = [v_a, v_b, v_c]^T \in \mathbb{R}^3$

#### Algebraic

None.

### External Variables

#### Differential

None.

#### Algebraic

None.

## Model Equations

### Differential Equations

An explicit representation for $\dot{\mathbf{v}}$ is not used because
the effective shunt admittances depend on connected components and are not
known at the bus level. The implicit DAE solver operates directly on
the accumulated KCL residual:

``` math
\begin{aligned}
0 &= \sum_{e \in \mathcal{E}} \mathbf{i}^\text{inj}_e
\end{aligned}
```

where $\mathbf{i}^\text{inj}_e$ is the vector of phase-current injections
of connected component $e$ into the bus, which are a function of the bus voltage and bus voltage derivative.

### Algebraic Equations

None.

## Initialization

For a balanced three-phase initialization derived from the phasor voltage
$V = |V| \angle \phi$ and nominal angular frequency $\omega_0 = 2 \pi f_0$,

``` math
\mathbf{v}(0) = \sqrt{2}\,|V|
\begin{bmatrix}
  \cos(\phi) \\
  \cos(\phi - \tfrac{2\pi}{3}) \\
  \cos(\phi + \tfrac{2\pi}{3})
\end{bmatrix}
```

and

``` math
\dot{\mathbf{v}}(0) = -\sqrt{2}\,|V|\,\omega_0
\begin{bmatrix}
  \sin(\phi) \\
  \sin(\phi - \tfrac{2\pi}{3}) \\
  \sin(\phi + \tfrac{2\pi}{3})
\end{bmatrix}
```

## Model Outputs

Phase voltages $v_a$, $v_b$, and $v_c$ are monitorable model outputs.

Phase-voltage derivatives $\dot{v}_a$, $\dot{v}_b$, and $\dot{v}_c$ are also available as monitorable outputs.
