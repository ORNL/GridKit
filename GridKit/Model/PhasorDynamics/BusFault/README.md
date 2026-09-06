# BusFault

Represents an impedance fault at a bus. This device can exist in two states, on or off, controlled by the user. Following a state change, generally the solver needs to be reset as this is a discrete event.

## Model Parameters

Symbol | Units     | JSON     | Description          | Typical Value | Note
-------|-----------|----------|----------------------|---------------|-----------------------------------------------------------------------------------
$R$    | [p.u.]    | `R`      | Fault resistance     |               |
$X$    | [p.u.]    | `X`      | Fault reactance      |               |
$U$    | [boolean] | `state0` | Initial fault status |               | JSON boolean; `true` puts the fault on. Changed at run time through `setStatus()`.

### Parameter Validation

None.

### Model Derived Parameters
```math
\begin{aligned}
  G   &=\dfrac{R}{R^2+ X^2} \\
  B   &= -\dfrac{X}{R^2 + X^2}\\
\end{aligned}
```

## Model Ports

Name             | Port  | Init  | Description
-----------------|-------|-------|------------
`bus`            | Bus   | Known | Required bus where the fault is applied
`control_signal` | Input | N/A   | Accepted by the parser but not read by the model; fault status is set through `state0` and `setStatus()`

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol      | Units   | Description                           | Note
------------|---------|---------------------------------------| ------
$I_r$       | [p.u.]  | Terminal current, real component      | Read by bus
$I_i$       | [p.u.]  | Terminal current, imaginary component |  Read by bus

### External Variables

#### Differential

None.

#### Algebraic

Symbol      | Units   | Description                           | Note
------------|---------|---------------------------------------| ------
$V_r$       | [p.u.]  | Terminal voltage, real component      | owned by bus object
$V_i$       | [p.u.]  | Terminal voltage, imaginary component | owned by bus object

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

```math
\begin{aligned}
0 &= -I_r + U (-G V_r + B V_i) \\
0 &= -I_i + U (-B V_r - G V_i)
\end{aligned}
```

### External Equations

The fault currents are added to the connected bus residuals:

```math
\begin{aligned}
I_r^{\mathrm{bus}} &\leftarrow I_r^{\mathrm{bus}} + I_r \\
I_i^{\mathrm{bus}} &\leftarrow I_i^{\mathrm{bus}} + I_i
\end{aligned}
```

## Initialization

The initial fault status and bus voltage determine the terminal currents:

```math
\begin{aligned}
U &\leftarrow \text{initial fault status} \\
I_r &\leftarrow U(-GV_r+BV_i) \\
I_i &\leftarrow U(-BV_r-GV_i)
\end{aligned}
```

The derivative vector entries initialize to zero.

## Monitors

Monitor | Units    | Description                               | Note
--------|----------|-------------------------------------------|-----
`state` | [binary] | Fault status                              | `1` when on; `0` when off
`ir`    | [p.u.]   | Fault-current real component              | Added to the connected-bus residual
`ii`    | [p.u.]   | Fault-current imaginary component         | Added to the connected-bus residual
