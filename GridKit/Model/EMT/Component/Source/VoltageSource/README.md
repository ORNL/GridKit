# VoltageSource Model

`VoltageSource` represents a three-phase voltage source in instantaneous abc
coordinates. The source waveform is configurable by phase magnitude and phase
offset for each phase and is otherwise constant. Each source port is
connected to the EMT bus through a phase resistance.

## Model Parameters

Symbol        | Units      | JSON     | Description                         | Note
--------------|------------|----------|-------------------------------------|-----
$E_a$         | [V]        | `Ea`     | Source voltage magnitude, phase a   | RMS
$E_b$         | [V]        | `Eb`     | Source voltage magnitude, phase b   | RMS
$E_c$         | [V]        | `Ec`     | Source voltage magnitude, phase c   | RMS
$\phi_a$      | [rad]      | `phia`   | Source phase offset, phase a        |
$\phi_b$      | [rad]      | `phib`   | Source phase offset, phase b        |
$\phi_c$      | [rad]      | `phic`   | Source phase offset, phase c        |
$\omega_0$    | [rad/s]    | `omega0` | Source angular frequency            |
$R_a$         | [$\Omega$] | `Ra`     | Terminal resistance, phase a        |
$R_b$         | [$\Omega$] | `Rb`     | Terminal resistance, phase b        |
$R_c$         | [$\Omega$] | `Rc`     | Terminal resistance, phase c        |

### Parameter Validation

```math
\begin{aligned}
E_a, E_b, E_c &\ge 0 \\
\omega_0 &> 0 \\
R_a, R_b, R_c &> 0
\end{aligned}
```

### Model Derived Parameters

None.

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
$\mathbf{v}$     | [V]   | Port voltage vector                         | $\mathbf{v} = [v_a, v_b, v_c]^T \in \mathbb{R}^3$

#### Algebraic

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{i}^{\mathrm{inj}}$ | `i` | Output | [A] | Current injection at source port | $\mathbf{i}^{\mathrm{inj}} \in \mathbb{R}^3$

## Model Equations

### Differential Equations

None.

### Algebraic Equations

None.

### Wiring

Given source angular frequency $\omega_0$, the source waveform is:

```math
\begin{aligned}
e_a(t) &= \sqrt{2}\,E_a\cos(\omega_0 t + \phi_a) \\
e_b(t) &= \sqrt{2}\,E_b\cos(\omega_0 t + \phi_b) \\
e_c(t) &= \sqrt{2}\,E_c\cos(\omega_0 t + \phi_c)
\end{aligned}
```

```math
\mathbf{i}^{\mathrm{inj}} =
\begin{bmatrix}
  \dfrac{e_a(t)-v_a}{R_a} \\
  \dfrac{e_b(t)-v_b}{R_b} \\
  \dfrac{e_c(t)-v_c}{R_c}
\end{bmatrix}
```

## Initialization

No internal state is initialized. At $t = 0$, the source waveform is:

```math
\begin{aligned}
e_a(0) &= \sqrt{2}\,E_a\cos(\phi_a) \\
e_b(0) &= \sqrt{2}\,E_b\cos(\phi_b) \\
e_c(0) &= \sqrt{2}\,E_c\cos(\phi_c)
\end{aligned}
```

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`e` | [V] | Source waveform | $\mathbf{e} \in \mathbb{R}^3$
