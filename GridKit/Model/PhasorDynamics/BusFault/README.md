# BusFault

Represents an impedance fault at a bus. The `active` input is fixed while a
solve is running and may be changed while the solve is stopped. The solver is
then reinitialized at the event time.

## Model Parameters

Symbol   | Units      | Description                     | Note
---------|------------|---------------------------------|-------
$R$      | [p.u.]     | Fault resistance                | 
$X$      | [p.u.]     | Fault reactance                 | 

## Model Inputs

Input    | Units | Description                                      | Default
---------|-------|--------------------------------------------------|--------
`active` | [-]   | Fault status; zero is inactive and nonzero active | 0

## Model Derived Parameters
``` math
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

Here $U$ is one when `active` is nonzero and zero otherwise.

## Initialization

When active, the terminal current is initialized from the bus voltage and
fault admittance. When inactive, both current components are initialized to
zero.

## Monitors

Name     | Description
---------|-----------------------------------
`active` | Fault status
`ir`     | Terminal current, real component
`ii`     | Terminal current, imaginary component
