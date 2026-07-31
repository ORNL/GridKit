# LoadZ

Static constant-impedance load model. `LoadZ` owns terminal current states and
adds their current contribution to the connected bus residual.

## Model Parameters

Symbol | Units  | JSON | Description     | Typical Value | Note
-------|--------|------|-----------------|---------------|-----
$R$    | [p.u.] | `R`  | Load resistance |               |
$X$    | [p.u.] | `X`  | Load reactance  |               |

### Parameter Validation

None.

### Model Derived Parameters

```math
\begin{aligned}
G &= \dfrac{R}{R^2 + X^2} \\
B &= -\dfrac{X}{R^2 + X^2}
\end{aligned}
```

## Model Inputs

Input    | Units | Description                                      | Default
---------|-------|--------------------------------------------------|--------
`online` | [-]   | Connection status; zero is disconnected          | 1

Any nonzero `online` value connects the load. The input is read when the
network contribution and monitors are evaluated. Disconnecting the load does
not remove its internal algebraic equations. Change the input only while the
solve is stopped.

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units  | Description                           | Note
-------|--------|---------------------------------------|------
$I_r$  | [p.u.] | Terminal current, real component      | Added to connected bus residual
$I_i$  | [p.u.] | Terminal current, imaginary component | Added to connected bus residual

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units  | Description                           | Note
-------|--------|---------------------------------------|------
$V_r$  | [p.u.] | Terminal voltage, real component      | Owned by connected bus
$V_i$  | [p.u.] | Terminal voltage, imaginary component | Owned by connected bus

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

```math
\begin{aligned}
0 &= I_r + G V_r - B V_i \\
0 &= I_i + B V_r + G V_i
\end{aligned}
```

### External Equations

```math
\begin{aligned}
I_r^{\mathrm{bus}} &\leftarrow I_r^{\mathrm{bus}} + I_r \\
I_i^{\mathrm{bus}} &\leftarrow I_i^{\mathrm{bus}} + I_i
\end{aligned}
```

## Initialization

The initial bus voltage determines the terminal currents:

```math
\begin{aligned}
I_r &\leftarrow -G V_r + B V_i \\
I_i &\leftarrow -B V_r - G V_i
\end{aligned}
```

The derivative vector entries initialize to zero.

## Monitors

Monitor | Units  | Description                                  | Note
--------|--------|----------------------------------------------|------
`p`  | [p.u.] | Active power at the connected bus terminal   | Positive for injection into the connected bus
`q`  | [p.u.] | Reactive power at the connected bus terminal | Positive for injection into the connected bus
