# VoltageSource Model

`VoltageSource` represents a three-phase voltage source in instantaneous abc
coordinates. The source waveform is configurable by phase magnitude and phase
offset for each phase and is otherwise constant. Each source terminal is
connected to the EMT bus through a phase resistance.

## Model Parameters

Symbol        | Units      | Description                         | Note
--------------|------------|-------------------------------------|--------------------------------
$E_a$         | [V]        | Source voltage magnitude, phase a   | RMS
$E_b$         | [V]        | Source voltage magnitude, phase b   | RMS
$E_c$         | [V]        | Source voltage magnitude, phase c   | RMS
$\phi_a$      | [rad]      | Source phase offset, phase a        |
$\phi_b$      | [rad]      | Source phase offset, phase b        |
$\phi_c$      | [rad]      | Source phase offset, phase c        |
$\omega_0$    | [rad/s]    | Source angular frequency            |
$R_a$         | [$\Omega$] | Terminal resistance, phase a        |
$R_b$         | [$\Omega$] | Terminal resistance, phase b        |
$R_c$         | [$\Omega$] | Terminal resistance, phase c        |

## Model Derived Parameters

None.

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

None.

### External Variables

External variables enter component model equations but are owned by
other components. The EMT bus at the source terminal owns the voltage
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

None.

### Algebraic Equations

None.

### Bus Residual Contributions

The source contributes current to the KCL residual at its terminal bus.
The injection vector is accumulated into the owning bus residual. Given source
angular frequency $\omega_0$, the source waveform is:

``` math
\begin{aligned}
e_a(t) &= \sqrt{2}\,E_a\cos(\omega_0 t + \phi_a) \\
e_b(t) &= \sqrt{2}\,E_b\cos(\omega_0 t + \phi_b) \\
e_c(t) &= \sqrt{2}\,E_c\cos(\omega_0 t + \phi_c)
\end{aligned}
```

The current contribution is positive into the bus:

``` math
\mathbf{i}^\text{inj} :=
\begin{bmatrix}
  \dfrac{e_a(t)-v_a}{R_a} \\
  \dfrac{e_b(t)-v_b}{R_b} \\
  \dfrac{e_c(t)-v_c}{R_c}
\end{bmatrix}
```

## Initialization

No internal state is initialized. At $t = 0$, the source waveform is:

``` math
\begin{aligned}
e_a(0) &= \sqrt{2}\,E_a\cos(\phi_a) \\
e_b(0) &= \sqrt{2}\,E_b\cos(\phi_b) \\
e_c(0) &= \sqrt{2}\,E_c\cos(\phi_c)
\end{aligned}
```

## Model Outputs

Candidate monitorable outputs include the source waveform components
$e_a(t)$, $e_b(t)$, and $e_c(t)$.

The terminal current injection expression is documented above as
$\mathbf{i}^\text{inj}$.
