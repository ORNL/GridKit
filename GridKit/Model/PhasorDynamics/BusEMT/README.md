# EMT Bus Model

`BusEMT` represents a three-phase bus in instantaneous abc coordinates. The
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
$v_a$    | [V]   | Bus voltage, phase a |
$v_b$    | [V]   | Bus voltage, phase b |
$v_c$    | [V]   | Bus voltage, phase c |

#### Algebraic

None.

### External Variables

#### Differential

None.

#### Algebraic

None.

## Model Equations

### Differential Equations

An explicit $\dot{\mathbf{v}} = \ldots$ form is not used because the
effective shunt admittances depends on connected components and is not
known at the bus level. The implicit DAE solver operates directly on
the accumulated KCL residual.

``` math
\begin{aligned}
0 &= \sum_{e \in \mathcal{E}} \Delta i_{a}^{e} \\
0 &= \sum_{e \in \mathcal{E}} \Delta i_{b}^{e} \\
0 &= \sum_{e \in \mathcal{E}} \Delta i_{c}^{e}
\end{aligned}
```

where $\Delta i_{a}^{e}$, $\Delta i_{b}^{e}$, and $\Delta i_{c}^{e}$ are the phase-current contributions
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
