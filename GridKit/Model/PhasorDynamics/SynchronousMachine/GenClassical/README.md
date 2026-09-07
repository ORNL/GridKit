# **Classical Generator Model (GENCLS)**

GENCLS is the 2nd order synchronous machine model: a constant transient EMF
behind the transient reactance. See the
[General Synchronous Machine Model](../README.md) for general synchronous
machine information.

## Notes

- No field dynamics: the `efd` input is the EMF behind $X'_d$ and is held
  constant when unconnected.
- No saliency and no saturation.

## Model Parameters

Symbol            | Units  | JSON  | Description                      | Typical Value | Note
------------------|--------|-------|----------------------------------|---------------|------
$P_0$             | [p.u.] | `p0`  | Initial active power injection   | 1.0           | System base; required initialization source
$Q_0$             | [p.u.] | `q0`  | Initial reactive power injection | 0.0           | System base; required initialization source
$S^\mathrm{base}$ | [MVA]  | `mva` | GENCLS component power base      | 100.0         |
$H$               | [sec]  | `H`   | Rotor inertia                    | 3.0           |
$D$               | [p.u.] | `D`   | Damping coefficient              | 0.0           |
$R_a$             | [p.u.] | `Ra`  | Armature resistance              | 0.0           | Component base
$X'_d$            | [p.u.] | `Xdp` | Direct-axis transient reactance  | 0.2           | Component base

### Parameter Validation

None.

### Model Derived Parameters

```math
\begin{aligned}
  G
    &= \dfrac{R_a}{R_a^2+(X'_d)^2} \\
  B
    &= -\dfrac{X'_d}{R_a^2+(X'_d)^2} \\
  k_\mathrm{base}
    &= \dfrac{S^\mathrm{sys}}{S^\mathrm{base}}
\end{aligned}
```

Multiplying by $k_\mathrm{base}$ converts system base to component base.
$S^\mathrm{sys}$ and $S^\mathrm{base}$ are stored in VA; $f^\mathrm{sys}$ is
the system frequency base in Hz.

## Model Ports

Name    | Port   | Init    | Description
--------|--------|---------|------------
`bus`   | Bus    | Known   | Terminal bus voltage
`pmech` | Input  | Unknown | Mechanical-power input
`efd`   | Input  | Unknown | Field-voltage input
`speed` | Output | Known   | Speed-deviation output

## Model Variables

### Internal Variables

#### Differential

Symbol   | Units  | Description         | Note
---------|--------|---------------------|------
$\delta$ | [rad]  | Rotor angle         |
$\omega$ | [p.u.] | Speed deviation     | Exported through `speed` when assigned

#### Algebraic

Symbol         | Units  | Description                           | Note
---------------|--------|---------------------------------------|------
$T_\mathrm{e}$ | [p.u.] | Electrical torque                     | Component base
$I_\mathrm{r}$ | [p.u.] | Terminal current, real component      | Component base
$I_\mathrm{i}$ | [p.u.] | Terminal current, imaginary component | Component base

### External Variables

#### Differential

None.

#### Algebraic

Symbol          | Units  | Init    | Description                           | Note
----------------|--------|---------|---------------------------------------|------
$V_\mathrm{r}$  | [p.u.] | Known   | Terminal voltage, real component      | Bus input
$V_\mathrm{i}$  | [p.u.] | Known   | Terminal voltage, imaginary component | Bus input
$P_\mathrm{m}$  | [p.u.] | Unknown | Mechanical power                      | Optional signal port `pmech`; system base
$E_\mathrm{fd}$ | [p.u.] | Unknown | Field voltage                         | Optional signal port `efd`; component base

## Model Equations

### Internal Equations

#### Differential

```math
\begin{aligned}
  0 &= -\dot\delta + 2\pi f^\mathrm{sys}\omega \\
  0 &= -\dot\omega + \dfrac{1}{2H}\left(\dfrac{k_\mathrm{base}P_\mathrm{m}-D\omega}{1+\omega}-T_\mathrm{e}\right)
\end{aligned}
```

#### Algebraic

```math
\begin{aligned}
  0 &= -T_\mathrm{e} + GE_\mathrm{fd}^2 - E_\mathrm{fd}[(GV_\mathrm{r}-BV_\mathrm{i})\cos(\delta)+(BV_\mathrm{r}+GV_\mathrm{i})\sin(\delta)] \\
  0 &= -I_\mathrm{r} + E_\mathrm{fd}(G\cos(\delta)-B\sin(\delta)) - GV_\mathrm{r} + BV_\mathrm{i} \\
  0 &= -I_\mathrm{i} + E_\mathrm{fd}(B\cos(\delta)+G\sin(\delta)) - BV_\mathrm{r} - GV_\mathrm{i}
\end{aligned}
```

### External Equations

```math
\begin{aligned}
  \Delta I^\mathrm{bus}_\mathrm{r} &\mathrel{+}= \dfrac{I_\mathrm{r}}{k_\mathrm{base}} \\
  \Delta I^\mathrm{bus}_\mathrm{i} &\mathrel{+}= \dfrac{I_\mathrm{i}}{k_\mathrm{base}}
\end{aligned}
```

## Initialization

### Input Initialization

```math
\begin{aligned}
  V_\mathrm{r}, V_\mathrm{i}
    &\leftarrow \text{terminal-bus voltage} \\
  P_0, Q_0
    &\leftarrow \text{power-flow injection on system base}
\end{aligned}
```

### Internal Initialization

All internal derivatives are initialized to zero:

```math
\begin{aligned}
  I_\mathrm{r}
    &\leftarrow k_\mathrm{base}\dfrac{V_\mathrm{r}P_0+V_\mathrm{i}Q_0}{V_\mathrm{r}^2+V_\mathrm{i}^2} \\
  I_\mathrm{i}
    &\leftarrow k_\mathrm{base}\dfrac{V_\mathrm{i}P_0-V_\mathrm{r}Q_0}{V_\mathrm{r}^2+V_\mathrm{i}^2} \\
  \delta
    &\leftarrow \text{arg}[V_\mathrm{r}+jV_\mathrm{i}+(R_a+jX'_d)(I_\mathrm{r}+jI_\mathrm{i})] \\
  \omega
    &\leftarrow 0 \\
  T_\mathrm{e}
    &\leftarrow V_\mathrm{r}I_\mathrm{r}+V_\mathrm{i}I_\mathrm{i}+R_a(I_\mathrm{r}^2+I_\mathrm{i}^2)
\end{aligned}
```

### Output Initialization

```math
\begin{aligned}
  P_\mathrm{m}
    &\leftarrow \dfrac{T_\mathrm{e}}{k_\mathrm{base}} \\
  E_\mathrm{fd}
    &\leftarrow |V_\mathrm{r}+jV_\mathrm{i}+(R_a+jX'_d)(I_\mathrm{r}+jI_\mathrm{i})|
\end{aligned}
```

## Monitors

Monitor | Units  | Description                           | Note
--------|--------|---------------------------------------|------
`ir`    | [p.u.] | Terminal current, real component      | System base; oriented leaving the machine
`ii`    | [p.u.] | Terminal current, imaginary component | System base; oriented leaving the machine
`p`     | [p.u.] | Active power                          | System base; $V_\mathrm{r}I_\mathrm{r}+V_\mathrm{i}I_\mathrm{i}$
`q`     | [p.u.] | Reactive power                        | System base; $V_\mathrm{i}I_\mathrm{r}-V_\mathrm{r}I_\mathrm{i}$
`delta` | [rad]  | Rotor angle                           |
`omega` | [p.u.] | Speed deviation                       | $\omega=0$ at synchronous speed
`speed` | [p.u.] | Per-unit machine speed                | $1+\omega$
